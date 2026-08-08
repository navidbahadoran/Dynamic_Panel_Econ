"""Deterministic, parallel, resume-safe Monte Carlo orchestration."""

from __future__ import annotations

import hashlib
import json
import os
import subprocess
import time
from concurrent.futures import ProcessPoolExecutor
from copy import deepcopy
from dataclasses import asdict
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd
from threadpoolctl import threadpool_limits

from .calibration import CalibrationResult, calibrate_cell, calibrate_rank_stress_cell
from .config import config_hash, write_resolved_config
from .dgp import (
    DGPParameters,
    generate_panel,
    generate_rank_stress_panel,
    group_gap_pilot,
)
from .estimation import (
    observe_fixed_rank_fits,
    observe_nuclear_fits,
)
from .inference import (
    infer_corrected_target,
    infer_target,
    prepare_riesz_system,
    prepare_split_fits,
)
from .lowrank import numerical_rank
from .mc_accounting import (
    FIT_DIAGNOSTIC_COLUMNS,
    INFERENCE_DIAGNOSTIC_COLUMNS,
    REPLICATION_COLUMNS,
    build_replication_records,
    canonical_status,
    method_name,
    semantic_replication_id,
)
from .rank_selection import (
    RankPilotFailure,
    RankSelectionFailure,
    fit_fixed_rank_multistart,
    select_ranks,
)
from .seeds import seed_sequence
from .targets import target_direction, target_regularity_diagnostics, target_value

FAILURE_CODES = (
    "success",
    "calibration_failure",
    "nonfinite_input",
    "full_fit_not_converged",
    "first_order_residual_high",
    "coefficient_bound_active",
    "rank_at_cap",
    "rank_pilot_failure",
    "rank_selection_failure",
    "target_unsupported_selected_rank",
    "split_target_unsupported_selected_rank",
    "tangent_gram_eigensolver_failure",
    "tangent_gram_nearly_singular",
    "split_tangent_gram_eigensolver_failure",
    "split_tangent_gram_nearly_singular",
    "riesz_not_converged",
    "riesz_target_instability",
    "split_fit_not_converged",
    "split_rank_loss",
    "split_target_support_loss",
    "split_riesz_target_instability",
    "spatial_variance_failure",
    "nonpositive_variance",
    "unexpected_exception",
)

Task = tuple[int, int, int, int, tuple[int, ...] | None]


def _json_default(value: Any) -> Any:
    if isinstance(value, np.generic):
        return value.item()
    if isinstance(value, np.ndarray):
        return value.tolist()
    if isinstance(value, tuple):
        return list(value)
    raise TypeError(f"Object of type {type(value).__name__} is not JSON serializable")


def _git_commit() -> str:
    try:
        return subprocess.check_output(
            ["git", "-c", "safe.directory=D:/Programming/Dynamic_Panel_Econ", "rev-parse", "HEAD"],
            text=True,
            stderr=subprocess.DEVNULL,
        ).strip()
    except (OSError, subprocess.SubprocessError):
        return "unknown"


def _params(config: dict[str, Any]) -> DGPParameters:
    return DGPParameters.from_mapping(config["dgp"])


def resolve_group_gap(config: dict[str, Any]) -> tuple[dict[str, Any], list[dict[str, Any]]]:
    """Run the prespecified DGP-4 pilot and freeze one global loading gap."""

    resolved = deepcopy(config)
    dgp_cfg = resolved["dgp"]
    params = _params(resolved)
    reports, chosen_gaps = [], []
    for n, t in resolved["run"]["cells"]:
        _, report = group_gap_pilot(
            int(n),
            int(t),
            int(resolved["run"]["master_seed"]),
            params,
            [float(value) for value in dgp_cfg["prespecified_gap_grid"]],
            float(dgp_cfg["min_abs_postscale_group_diff_a"]),
            bool(dgp_cfg["auto_adjust_group_gap"]),
            draws=int(dgp_cfg.get("group_gap_pilot_draws", 20)),
        )
        report.update({"N": int(n), "T": int(t)})
        reports.append(report)
        chosen_gaps.append(float(report["chosen_gap"]))
    chosen = max(chosen_gaps)
    center = 0.5 * (params.mu_lambda_a_1 + params.mu_lambda_a_2)
    dgp_cfg["mu_lambda_a_1"] = center - chosen / 2.0
    dgp_cfg["mu_lambda_a_2"] = center + chosen / 2.0
    dgp_cfg["resolved_group_gap"] = chosen
    return resolved, reports


def calibrate_design(
    config: dict[str, Any],
    *,
    failures: dict[tuple[int, int, int, tuple[int, ...] | None], str] | None = None,
) -> dict[tuple[int, int, int, tuple[int, ...] | None], CalibrationResult]:
    params = _params(config)
    calibration_seed = config["dgp"].get("calibration_seed")
    calibration_seed = (
        int(config["run"]["master_seed"])
        if calibration_seed is None
        else int(calibration_seed)
    )
    results = {}
    for dgp in config["run"]["dgps"]:
        for n, t in config["run"]["cells"]:
            for true_rank in _task_designs(config):
                common = {
                    "params": params,
                    "pi_h": float(config["dgp"]["pi_h"]),
                    "target_r2": float(config["dgp"]["target_r2"]),
                    "draws": int(config["dgp"]["calibration_draws"]),
                    "tolerance": float(config["dgp"]["calibration_tolerance"]),
                }
                key = (int(dgp), int(n), int(t), true_rank)
                try:
                    if true_rank is None:
                        result = calibrate_cell(
                            int(dgp), int(n), int(t), calibration_seed, **common
                        )
                    else:
                        result = calibrate_rank_stress_cell(
                            int(dgp), int(n), int(t), true_rank, calibration_seed,
                            component_strengths=tuple(config["rank_stress"]["component_strengths"]),
                            **common,
                        )
                    results[key] = result
                except Exception as exc:
                    if failures is None:
                        raise
                    failures[key] = f"{type(exc).__name__}: {exc}"
    return results


def _selection_options(config: dict[str, Any]) -> dict[str, Any]:
    estimation = config["estimation"]
    keys = (
        "coefficient_bound",
        "max_sweeps",
        "objective_rtol",
        "stationarity_tol",
        "lstsq_rcond",
        "nuclear_gamma",
        "nuclear_epsilon",
        "nuclear_max_iter",
        "nuclear_tol",
        "dykstra_max_iter",
        "dykstra_tol",
        "eta_for_penalty",
        "spatial_dimension",
        "ic_multiplier",
        "threshold_multiplier",
        "start_objective_stability_tol",
        "rank_adaptive_improvement_tol",
        "rank_adaptive_removal_tol",
        "rank_adaptive_max_steps",
        "rank_adaptive_max_routes",
        "cap_pilot_start_envelope_fraction",
        "dense_nuclear_gamma",
        "threshold_sensitivity_multipliers",
        "ic_sensitivity_multipliers",
        "larger_rank_caps",
        "compute_rank_sensitivities",
        "compute_dense_grid_sensitivity",
        "compute_larger_cap_sensitivity",
    )
    return {key: estimation[key] for key in keys}


def _fit_options(config: dict[str, Any]) -> dict[str, Any]:
    estimation = config["estimation"]
    return {
        key: estimation[key]
        for key in (
            "coefficient_bound",
            "max_sweeps",
            "objective_rtol",
            "stationarity_tol",
            "lstsq_rcond",
        )
    }


def _task_parts(task: Task) -> tuple[int, int, int, int, tuple[int, ...]]:
    dgp, n, t, replication, stress_rank = task
    return dgp, n, t, replication, tuple(stress_rank or (1, 1, 1))


def _failure_record(
    task: Task,
    config: dict[str, Any],
    status: str,
    detail: str,
    extra: dict[str, Any] | None = None,
) -> dict[str, Any]:
    dgp, n, t, replication, true_rank = _task_parts(task)
    method = method_name(config["run"]["rank_mode"])
    record = {
        "record_type": "failure",
        "dgp": dgp,
        "N": n,
        "T": t,
        "replication": replication,
        "target": "__replication__",
        "true_rank_vector": json.dumps(true_rank),
        "status": status,
        "primary_status": canonical_status(status),
        "method": method,
        "run_id": config_hash(config),
        "semantic_replication_id": semantic_replication_id(
            dgp, n, t, replication, true_rank
        ),
        "supplied_rank_vector": (
            json.dumps(true_rank if task[4] is not None else config["estimation"]["fixed_ranks"])
            if method == "fixed_rank"
            else None
        ),
        "failure_detail": detail,
        "config_hash": config_hash(config),
        "git_commit": _git_commit(),
    }
    if extra:
        record.update(extra)
    return record


def _rank_record(
    task: Task,
    config: dict[str, Any],
    selection: Any,
    rank_runtime: float,
    status: str,
    panel_diagnostics: dict[str, Any] | None = None,
) -> dict[str, Any]:
    dgp, n, t, replication, true_rank = _task_parts(task)
    diagnostics = selection.diagnostics
    selected = tuple(selection.selected.ranks)
    candidates = {tuple(rank) for rank in diagnostics["candidate_rank_vectors"]}
    true_record = next(
        (
            item
            for item in diagnostics["baseline_candidate_records"]
            if tuple(item["ranks"]) == true_rank
        ),
        None,
    )
    true_sources = list(true_record["sources"]) if true_record is not None else []
    sensitivity = diagnostics.get("sensitivities", {})
    selected_starts = selection.selected.start_diagnostics
    cap_pilot = diagnostics["cap_pilot"]
    record: dict[str, Any] = {
        "record_type": "rank",
        "method": "selected_rank",
        "run_id": config_hash(config),
        "primary_status": canonical_status(status),
        "semantic_replication_id": semantic_replication_id(
            dgp, n, t, replication, true_rank
        ),
        "dgp": dgp,
        "N": n,
        "T": t,
        "replication": replication,
        "status": status,
        "true_rank_vector": json.dumps(true_rank),
        "selected_rank_vector": json.dumps(selected),
        "rank_cap_thresholded_vector": json.dumps(diagnostics["rank_cap_thresholded_vector"]),
        "cap_pilot_rank": json.dumps(diagnostics["rank_cap_thresholded_vector"]),
        "true_rank_in_candidates": true_rank in candidates,
        "candidate_coverage": true_rank in candidates,
        "rank_selection_diagnostics": (
            json.dumps(diagnostics, sort_keys=True, default=_json_default)
            if config["run"].get("save_candidate_details", True)
            else None
        ),
        "true_rank_sources": json.dumps(true_sources),
        "true_rank_from_nuclear_path": any(source.startswith("nuclear_path_") for source in true_sources),
        "true_rank_from_rank_cap_pilot": "rank_cap_pilot" in true_sources,
        "true_rank_from_neighbor_completion": any("neighbor" in source for source in true_sources),
        "true_rank_from_multiple_sources": len(true_sources) > 1,
        "true_rank_ic": true_record["ic"] if true_record is not None else np.nan,
        "selected_ic": selection.selected.ic,
        "selected_minus_true_ic": (
            selection.selected.ic - true_record["ic"] if true_record is not None else np.nan
        ),
        "exact_rank_recovery": selected == true_rank,
        "zero_rank_recovery": all(
            selected[index] == 0 for index, rank in enumerate(true_rank) if rank == 0
        ),
        "candidate_count_initial": diagnostics["candidate_count_initial"],
        "candidate_count_final": diagnostics["candidate_count_final"],
        "smallest_ic": diagnostics["smallest_ic"],
        "second_smallest_ic": diagnostics["second_smallest_ic"],
        "ic_gap": diagnostics["ic_gap"],
        "rank_at_cap": diagnostics["selected_rank_at_cap"],
        "selected_converged": selection.selected.fit.converged,
        "selected_stationarity_residual": selection.selected.fit.stationarity_residual,
        "selected_max_envelope_ratio": selection.selected.fit.max_envelope_ratio,
        "selected_start_objectives": json.dumps(
            selected_starts["start_objectives"], default=_json_default
        ),
        "selected_start_stationarity_residuals": json.dumps(
            selected_starts["start_stationarity_residuals"], default=_json_default
        ),
        "selected_best_two_objective_gap": selected_starts["best_two_objective_gap"],
        "selected_objective_stability_pass": selected_starts["objective_stability_pass"],
        "selected_third_start_used": selected_starts["third_start_used"],
        "cap_pilot_converged": diagnostics["cap_pilot_converged"],
        "cap_pilot_stationarity_residual": diagnostics["cap_pilot_stationarity_residual"],
        "cap_pilot_max_envelope_ratio": diagnostics["cap_pilot_max_envelope_ratio"],
        "cap_pilot_numerical_rank_before_thresholding": json.dumps(
            cap_pilot["numerical_rank_before_thresholding"]
        ),
        "cap_pilot_outer_start_rank_vectors": json.dumps(
            cap_pilot["outer_start_rank_vectors"]
        ),
        "cap_pilot_outer_final_rank_vectors": json.dumps(
            cap_pilot["outer_final_rank_vectors"]
        ),
        "cap_pilot_outer_objectives": json.dumps(
            cap_pilot["outer_objectives"], default=_json_default
        ),
        "cap_pilot_best_two_objective_gap": cap_pilot["best_two_objective_gap"],
        "cap_pilot_objective_stability_pass": cap_pilot["objective_stability_pass"],
        "cap_pilot_attempted_route_count": cap_pilot.get("attempted_route_count", 0),
        "cap_pilot_valid_route_count": cap_pilot.get("valid_route_count", 0),
        "cap_pilot_stable_route_count": cap_pilot.get("stable_route_count", 0),
        "cap_pilot_stable_final_numerical_ranks_agree": cap_pilot.get(
            "stable_final_numerical_ranks_agree", False
        ),
        "cap_pilot_stable_final_thresholded_ranks_agree": cap_pilot.get(
            "stable_final_thresholded_ranks_agree", False
        ),
        "cap_pilot_start_attempts": json.dumps(
            cap_pilot["outer_start_attempts"], sort_keys=True, default=_json_default
        ),
        "cap_pilot_confirmation_attempts": json.dumps(
            cap_pilot.get("basin_confirmation_attempts", []),
            sort_keys=True,
            default=_json_default,
        ),
        "cap_pilot_original_best_objective": cap_pilot.get("original_best_objective"),
        "cap_pilot_original_second_best_objective": cap_pilot.get(
            "original_second_best_objective"
        ),
        "cap_pilot_original_stability_gap": cap_pilot.get("original_stability_gap"),
        "cap_pilot_confirmation_best_objective": cap_pilot.get(
            "confirmation_best_objective"
        ),
        "cap_pilot_number_confirmation_valid": cap_pilot.get(
            "number_confirmation_valid", 0
        ),
        "cap_pilot_number_confirmation_matching_best": cap_pilot.get(
            "number_confirmation_matching_best", 0
        ),
        "cap_pilot_final_acceptance_basis": cap_pilot.get(
            "final_pilot_acceptance_basis"
        ),
        "rank_runtime_seconds": rank_runtime,
        "rank_diagnostics_json": (
            json.dumps(diagnostics, sort_keys=True, default=_json_default)
            if config["run"].get("save_candidate_details", True)
            else None
        ),
        "config_hash": config_hash(config),
        "git_commit": _git_commit(),
    }
    for multiplier, result in sensitivity.get("ic_multipliers", {}).items():
        record[f"ic_multiplier_{multiplier}_selected_rank"] = json.dumps(
            result["selected_rank"]
        )
        record[f"ic_multiplier_{multiplier}_true_rank_in_candidates"] = (
            true_rank
            in {tuple(rank) for rank in result.get("candidate_rank_vectors", [])}
        )
    for multiplier, result in sensitivity.get("threshold_multipliers", {}).items():
        record[f"threshold_multiplier_{multiplier}_selected_rank"] = json.dumps(
            result["selected_rank"]
        )
        record[f"threshold_multiplier_{multiplier}_true_rank_in_candidates"] = (
            true_rank
            in {tuple(rank) for rank in result.get("candidate_rank_vectors", [])}
        )
    for index, block in enumerate(("A", "B", "H")):
        record[f"{block}_underselected"] = selected[index] < true_rank[index]
        record[f"{block}_overselected"] = selected[index] > true_rank[index]
        record[f"{block}_true_rank"] = true_rank[index]
        record[f"{block}_selected_rank"] = selected[index]
    if panel_diagnostics:
        record.update(
            {key: value for key, value in panel_diagnostics.items() if np.isscalar(value)}
        )
    return record


def _fixed_rank_record(
    task: Task,
    config: dict[str, Any],
    fit: Any,
    runtime: float,
    panel_diagnostics: dict[str, Any],
) -> dict[str, Any]:
    dgp, n, t, replication, true_rank = _task_parts(task)
    supplied = tuple(int(value) for value in fit.ranks)
    record = {
        "record_type": "rank",
        "method": "fixed_rank",
        "run_id": config_hash(config),
        "primary_status": "success",
        "dgp": dgp,
        "N": n,
        "T": t,
        "replication": replication,
        "status": "success",
        "true_rank_vector": json.dumps(true_rank),
        "supplied_rank_vector": json.dumps(supplied),
        "selected_rank_vector": None,
        "cap_pilot_rank": None,
        "candidate_coverage": None,
        "rank_selection_diagnostics": None,
        "rank_runtime_seconds": runtime,
        "full_objective": fit.objective,
        "full_iterations": fit.iterations,
        "selected_converged": fit.converged,
        "selected_stationarity_residual": fit.stationarity_residual,
        "selected_max_envelope_ratio": fit.max_envelope_ratio,
        "semantic_replication_id": semantic_replication_id(
            dgp, n, t, replication, true_rank
        ),
        "config_hash": config_hash(config),
        "git_commit": _git_commit(),
    }
    record.update(
        {key: value for key, value in panel_diagnostics.items() if np.isscalar(value)}
    )
    return record


def _fit_diagnostic_record(
    task: Task,
    config: dict[str, Any],
    fit: Any,
    *,
    fit_type: str,
    runtime_seconds: float | None = None,
    initialization_route: str | None = None,
    start_number: int | None = None,
    candidate_source: str | None = None,
    ic: float | None = None,
    ic_valid: bool | None = None,
) -> dict[str, Any]:
    dgp, n, t, replication, true_rank = _task_parts(task)
    singular = [np.linalg.svd(matrix, compute_uv=False) for matrix in fit.theta.matrices()]
    positive = [
        (float(values[0]), float(values[rank - 1]))
        for values, rank in zip(singular, fit.ranks, strict=True)
        if rank > 0 and len(values) >= rank
    ]
    sigma_1 = min((item[0] for item in positive), default=np.nan)
    sigma_r = min((item[1] for item in positive), default=np.nan)
    objective_initial = fit.objective_history[0] if fit.objective_history else np.nan
    relative = (
        abs(objective_initial - fit.objective) / max(1.0, abs(objective_initial))
        if np.isfinite(objective_initial)
        else np.nan
    )
    return {
        "record_type": "fit_diagnostic",
        "run_id": config_hash(config),
        "semantic_replication_id": semantic_replication_id(
            dgp, n, t, replication, true_rank
        ),
        "dgp": dgp,
        "N": n,
        "T": t,
        "replication": replication,
        "method": method_name(config["run"]["rank_mode"]),
        "fit_type": fit_type,
        "target": None,
        "requested_rank": json.dumps(fit.ranks),
        "numerical_rank": json.dumps(
            tuple(numerical_rank(matrix) for matrix in fit.theta.matrices())
        ),
        "initialization_route": (
            initialization_route
            if initialization_route is not None
            else fit.diagnostics.get("diagnostic_context")
        ),
        "start_number": start_number,
        "objective_initial": objective_initial,
        "objective_final": fit.objective,
        "relative_objective_change": relative,
        "iterations": fit.iterations,
        "convergence_flag": fit.converged,
        "iteration_cap_hit": fit.iterations >= int(config["estimation"]["max_sweeps"]),
        "stationarity_residual": fit.stationarity_residual,
        "stationarity_pass": fit.stationarity_residual
        <= float(config["estimation"]["stationarity_tol"]),
        "coefficient_envelope": fit.max_envelope_ratio
        * float(config["estimation"]["coefficient_bound"]),
        "coefficient_envelope_ratio": fit.max_envelope_ratio,
        "coefficient_bound_hit": fit.max_envelope_ratio >= 1.0,
        "sigma_1": sigma_1,
        "sigma_r": sigma_r,
        "sigma_r_over_sigma_1": sigma_r / max(sigma_1, np.finfo(float).tiny)
        if positive
        else np.nan,
        "best_start_objective": None,
        "second_start_objective": None,
        "objective_stability_gap": None,
        "objective_stability_pass": None,
        "runtime_seconds": (
            runtime_seconds
            if runtime_seconds is not None
            else fit.diagnostics.get("runtime_seconds")
        ),
        "exception_type": None,
        "exception_message": None,
        "nuclear_path_index": None,
        "lambda": None,
        "thresholded_rank": None,
        "candidate_source": candidate_source,
        "IC": ic,
        "IC_valid": ic_valid,
        "diagnostic_context": fit.diagnostics.get("diagnostic_context"),
        "initial_coefficient_envelope": fit.diagnostics.get(
            "initial_coefficient_envelope"
        ),
        "final_coefficient_envelope": fit.diagnostics.get(
            "final_coefficient_envelope"
        ),
        "coefficient_envelope_history": json.dumps(
            fit.diagnostics.get("coefficient_envelope_history"),
            default=_json_default,
        ),
    }


def _nuclear_fit_diagnostic_record(
    task: Task, config: dict[str, Any], fit: Any, path_index: int
) -> dict[str, Any]:
    dgp, n, t, replication, true_rank = _task_parts(task)
    numerical = tuple(numerical_rank(matrix) for matrix in fit.theta.matrices())
    initial = fit.objective_history[0] if fit.objective_history else np.nan
    return {
        "record_type": "fit_diagnostic", "run_id": config_hash(config),
        "semantic_replication_id": semantic_replication_id(
            dgp, n, t, replication, true_rank
        ),
        "dgp": dgp, "N": n, "T": t, "replication": replication,
        "method": "selected_rank", "fit_type": "nuclear_path", "target": None,
        "requested_rank": None, "numerical_rank": json.dumps(numerical),
        "initialization_route": "warm_nuclear_path" if path_index else "zero",
        "start_number": path_index + 1, "objective_initial": initial,
        "objective_final": fit.objective,
        "relative_objective_change": abs(initial - fit.objective) / max(1.0, abs(initial)),
        "iterations": fit.iterations, "convergence_flag": fit.converged,
        "iteration_cap_hit": fit.iterations >= int(config["estimation"]["nuclear_max_iter"]),
        "stationarity_residual": None, "stationarity_pass": None,
        "coefficient_envelope": max(
            float(np.max(np.abs(matrix))) for matrix in fit.theta.matrices()
        ),
        "coefficient_envelope_ratio": max(
            float(np.max(np.abs(matrix))) for matrix in fit.theta.matrices()
        ) / float(config["estimation"]["coefficient_bound"]),
        "coefficient_bound_hit": max(
            float(np.max(np.abs(matrix))) for matrix in fit.theta.matrices()
        ) >= float(config["estimation"]["coefficient_bound"]),
        "sigma_1": None, "sigma_r": None, "sigma_r_over_sigma_1": None,
        "best_start_objective": None, "second_start_objective": None,
        "objective_stability_gap": None, "objective_stability_pass": None,
        "runtime_seconds": fit.runtime_seconds,
        "exception_type": None, "exception_message": None,
        "nuclear_path_index": path_index, "lambda": fit.penalty,
        "thresholded_rank": None, "candidate_source": "nuclear_screening",
        "IC": None, "IC_valid": None,
    }


def _selected_fit_diagnostic_records(
    task: Task, config: dict[str, Any], selection: Any
) -> list[dict[str, Any]]:
    rows: list[dict[str, Any]] = []
    for index, nuclear in enumerate(selection.nuclear_fits):
        rows.append(
            {
                **{key: None for key in (
                    "requested_rank", "numerical_rank", "initialization_route", "start_number",
                    "objective_initial", "relative_objective_change", "stationarity_residual",
                    "stationarity_pass", "coefficient_envelope", "coefficient_envelope_ratio",
                    "coefficient_bound_hit", "sigma_1", "sigma_r", "sigma_r_over_sigma_1",
                    "best_start_objective", "second_start_objective", "objective_stability_gap",
                    "objective_stability_pass", "runtime_seconds", "exception_type",
                    "exception_message", "thresholded_rank", "candidate_source", "IC", "IC_valid",
                )},
                "record_type": "fit_diagnostic",
                "run_id": config_hash(config),
                "dgp": task[0], "N": task[1], "T": task[2], "replication": task[3],
                "method": "selected_rank", "fit_type": "nuclear_path", "target": None,
                "objective_final": nuclear.objective, "iterations": nuclear.iterations,
                "convergence_flag": nuclear.converged, "iteration_cap_hit": False,
                "nuclear_path_index": index, "lambda": nuclear.penalty,
            }
        )
    for route in selection.diagnostics.get("cap_pilot", {}).get(
        "outer_start_attempts", []
    ):
        path = route.get("path", [])
        for step_number, step in enumerate(path, 1):
            rows.append(
                {
                    "record_type": "fit_diagnostic",
                    "run_id": config_hash(config),
                    "dgp": task[0], "N": task[1], "T": task[2],
                    "replication": task[3], "method": "selected_rank",
                    "fit_type": "rank_cap_pilot", "target": None,
                    "requested_rank": json.dumps(step.get("ranks")),
                    "numerical_rank": json.dumps(route.get("final_numerical_rank_vector"))
                    if step_number == len(path) else None,
                    "initialization_route": route.get("route_source"),
                    "start_number": route.get("route_number"),
                    "objective_initial": path[0].get("objective") if path else None,
                    "objective_final": step.get("objective"),
                    "relative_objective_change": None,
                    "iterations": step_number, "convergence_flag": step.get("valid"),
                    "iteration_cap_hit": False,
                    "stationarity_residual": route.get("final_stationarity_residual")
                    if step_number == len(path) else None,
                    "stationarity_pass": route.get("final_valid")
                    if step_number == len(path) else None,
                    "coefficient_envelope": None,
                    "coefficient_envelope_ratio": route.get("final_max_envelope_ratio")
                    if step_number == len(path) else None,
                    "coefficient_bound_hit": (
                        route.get("final_max_envelope_ratio", 0.0) >= 1.0
                        if step_number == len(path) else None
                    ),
                    "sigma_1": None, "sigma_r": None, "sigma_r_over_sigma_1": None,
                    "best_start_objective": None, "second_start_objective": None,
                    "objective_stability_gap": route.get("start_fit_diagnostics", {}).get("best_two_objective_gap"),
                    "objective_stability_pass": route.get("start_fit_diagnostics", {}).get("objective_stability_pass"),
                    "runtime_seconds": None, "exception_type": None,
                    "exception_message": None,
                    "nuclear_path_index": route.get("nuclear_path_index"),
                    "lambda": route.get("nuclear_penalty"),
                    "thresholded_rank": json.dumps(route.get("final_thresholded_rank_vector"))
                    if step_number == len(path) else None,
                    "candidate_source": step.get("move"), "IC": None, "IC_valid": None,
                }
            )
    for candidate in selection.candidates.values():
        starts = candidate.start_diagnostics
        objectives = list(starts.get("start_objectives", []))
        residuals = list(starts.get("start_stationarity_residuals", []))
        valid_starts = list(starts.get("start_valid", []))
        sorted_objectives = sorted(float(value) for value in objectives if np.isfinite(value))
        for start_index, objective in enumerate(objectives):
            row = _fit_diagnostic_record(
                task,
                config,
                candidate.fit,
                fit_type="candidate_post_refit",
                initialization_route=(
                    "deterministic" if start_index < 2 else "third_randomized"
                ),
                start_number=start_index + 1,
                candidate_source="|".join(candidate.sources),
                ic=candidate.ic,
                ic_valid=candidate.valid,
            )
            row["objective_final"] = objective
            row["stationarity_residual"] = (
                residuals[start_index] if start_index < len(residuals) else None
            )
            row["stationarity_pass"] = (
                valid_starts[start_index] if start_index < len(valid_starts) else None
            )
            row["best_start_objective"] = (
                sorted_objectives[0] if sorted_objectives else None
            )
            row["second_start_objective"] = (
                sorted_objectives[1] if len(sorted_objectives) > 1 else None
            )
            row["objective_stability_gap"] = starts.get("best_two_objective_gap")
            row["objective_stability_pass"] = starts.get("objective_stability_pass")
            rows.append(row)
    return rows


def classify_inference_status(result: Any, config: dict[str, Any]) -> str:
    """Apply the paper's primary and split-inference failure precedence."""

    if getattr(result, "failure_code", None) is not None:
        return str(result.failure_code)
    if not result.riesz.converged:
        return "riesz_not_converged"
    floor = float(config["inference"]["riesz_target_rayleigh_floor"])
    if (
        not np.isfinite(result.riesz.target_rayleigh_quotient)
        or result.riesz.target_rayleigh_quotient < floor
    ):
        return "riesz_target_instability"
    if result.variance <= 0.0 or not np.isfinite(result.variance):
        return "nonpositive_variance"
    if not result.corrected:
        return "success"
    split_fits = result.diagnostics["split_fits"]
    if len(split_fits) != 4 or not all(
        item["converged"]
        and item["stationarity_residual"] <= config["estimation"]["stationarity_tol"]
        and item["max_envelope_ratio"] < 1.0
        for item in split_fits
    ):
        return "split_fit_not_converged"
    if not all(item["rank_supported"] for item in split_fits):
        return "split_rank_loss"
    if not all(item["target_supported"] for item in split_fits):
        return "split_target_support_loss"
    if not all(
        item["riesz_converged"] and item["riesz_target_stable"] for item in split_fits
    ):
        return "split_riesz_target_instability"
    return "success"


def _use_spatial_variance(config: dict[str, Any], dgp: int) -> bool:
    variance_type = config["inference"].get("variance_type", "auto")
    if variance_type == "spatial":
        return True
    if variance_type == "diagonal":
        return False
    return dgp >= 2


def _dgp_realization_hash(panel: Any, calibration: dict[str, Any]) -> str:
    """Hash the complete realized estimation sample and frozen DGP calibration."""

    digest = hashlib.sha256()
    arrays = (
        ("y", panel.y),
        *((f"regressor_{index}", value) for index, value in enumerate(panel.design.regressors())),
        *((f"theta0_{index}", value) for index, value in enumerate(panel.theta0.matrices())),
        ("u", panel.u),
        ("u_tilde", panel.u_tilde),
        ("u_tilde_lag", panel.u_tilde_lag),
        ("groups", panel.groups),
    )
    for name, value in arrays:
        array = np.ascontiguousarray(value)
        digest.update(name.encode("ascii"))
        digest.update(str(array.dtype).encode("ascii"))
        digest.update(np.asarray(array.shape, dtype=np.int64).tobytes())
        digest.update(array.tobytes())
    for name in ("c_h", "c_xi"):
        digest.update(name.encode("ascii"))
        digest.update(np.float64(calibration[name]).tobytes())
    return digest.hexdigest()


def run_replication(
    task: Task,
    config: dict[str, Any],
    calibration: dict[str, Any],
) -> list[dict[str, Any]]:
    dgp, n, t, replication, true_rank = _task_parts(task)
    started = time.perf_counter()
    try:
        params = _params(config)
        dgp_seed = seed_sequence(
            config["run"]["master_seed"], "production", dgp, n, t, replication, true_rank, "dgp"
        )
        if task[4] is None:
            panel = generate_panel(
                dgp,
                n,
                t,
                dgp_seed,
                c_h=float(calibration["c_h"]),
                c_xi=float(calibration["c_xi"]),
                params=params,
                coefficient_bound=float(config["estimation"]["coefficient_bound"]),
                simulation_interior_margin=float(
                    config["estimation"]["simulation_interior_margin"]
                ),
            )
        else:
            panel = generate_rank_stress_panel(
                dgp,
                n,
                t,
                true_rank,
                dgp_seed,
                component_strengths=tuple(config["rank_stress"]["component_strengths"]),
                c_h=float(calibration["c_h"]),
                c_xi=float(calibration["c_xi"]),
                params=params,
                coefficient_bound=float(config["estimation"]["coefficient_bound"]),
                simulation_interior_margin=float(
                    config["estimation"]["simulation_interior_margin"]
                ),
            )
        arrays = [panel.y, *panel.design.regressors(), *panel.theta0.matrices(), panel.u]
        if not all(np.all(np.isfinite(array)) for array in arrays):
            return [_failure_record(task, config, "nonfinite_input", "generated array is nonfinite")]
        panel.diagnostics["dgp_realization_hash"] = _dgp_realization_hash(panel, calibration)

        selection_started = time.perf_counter()
        selection = None
        if config["run"]["rank_mode"] == "fixed":
            supplied_ranks = (
                true_rank
                if task[4] is not None
                else tuple(int(value) for value in config["estimation"]["fixed_ranks"])
            )
            fit, fixed_rank_verification = fit_fixed_rank_multistart(
                panel.y,
                panel.design,
                supplied_ranks,
                seed=seed_sequence(
                    config["run"]["master_seed"],
                    dgp,
                    n,
                    t,
                    replication,
                    true_rank,
                    "fixed_rank_starts",
                ),
                fit_options=_fit_options(config),
                stationarity_tolerance=float(config["estimation"]["stationarity_tol"]),
                start_objective_stability_tol=float(
                    config["estimation"]["start_objective_stability_tol"]
                ),
                start_envelope_fraction=float(
                    config["estimation"]["cap_pilot_start_envelope_fraction"]
                ),
            )
            rank_runtime = time.perf_counter() - selection_started
            rank_status = "success"
            rank_row = _fixed_rank_record(
                task, config, fit, rank_runtime, panel.diagnostics
            )
            rank_row["fixed_rank_multistart_diagnostics"] = json.dumps(
                fixed_rank_verification, sort_keys=True, default=_json_default
            )
        else:
            try:
                selection = select_ranks(
                    panel.y,
                    panel.design,
                    tuple(int(value) for value in config["estimation"]["rank_caps"]),
                    seed=seed_sequence(
                        config["run"]["master_seed"], dgp, n, t, replication, true_rank, "rank_starts"
                    ),
                    **_selection_options(config),
                )
            except RankPilotFailure as exc:
                pilot = exc.diagnostics
                return [
                    _failure_record(
                        task,
                        config,
                        "rank_pilot_failure",
                        str(exc),
                        {
                            "cap_pilot_attempted_route_count": pilot.get(
                                "attempted_route_count", 0
                            ),
                            "cap_pilot_valid_route_count": pilot.get("valid_route_count", 0),
                            "cap_pilot_stable_route_count": pilot.get(
                                "stable_route_count", 0
                            ),
                            "cap_pilot_objective_stability_pass": pilot.get(
                                "objective_stability_pass", False
                            ),
                            "cap_pilot_start_attempts": json.dumps(
                                pilot.get("outer_start_attempts", []),
                                sort_keys=True,
                                default=_json_default,
                            ),
                            "cap_pilot_confirmation_attempts": json.dumps(
                                pilot.get("basin_confirmation_attempts", []),
                                sort_keys=True,
                                default=_json_default,
                            ),
                            "cap_pilot_original_best_objective": pilot.get(
                                "original_best_objective"
                            ),
                            "cap_pilot_original_second_best_objective": pilot.get(
                                "original_second_best_objective"
                            ),
                            "cap_pilot_original_stability_gap": pilot.get(
                                "original_stability_gap"
                            ),
                            "cap_pilot_confirmation_best_objective": pilot.get(
                                "confirmation_best_objective"
                            ),
                            "cap_pilot_number_confirmation_valid": pilot.get(
                                "number_confirmation_valid", 0
                            ),
                            "cap_pilot_number_confirmation_matching_best": pilot.get(
                                "number_confirmation_matching_best", 0
                            ),
                            "cap_pilot_final_acceptance_basis": pilot.get(
                                "final_pilot_acceptance_basis", "failure"
                            ),
                            "replication_runtime_seconds": time.perf_counter() - started,
                            **{
                                key: value
                                for key, value in panel.diagnostics.items()
                                if np.isscalar(value)
                            },
                        },
                    )
                ]
            except RankSelectionFailure as exc:
                detail = str(exc)
                status = (
                    "candidate_numerically_unresolved"
                    if "valid candidate post-refit" in detail
                    or "stable" in detail.lower()
                    else "rank_selection_failure"
                )
                return [
                    _failure_record(
                        task,
                        config,
                        status,
                        detail,
                        {"dgp_realization_hash": panel.diagnostics["dgp_realization_hash"]},
                    )
                ]
            rank_runtime = time.perf_counter() - selection_started
            rank_status = "rank_at_cap" if selection.diagnostics["selected_rank_at_cap"] else "success"
            rank_row = _rank_record(
                task,
                config,
                selection,
                rank_runtime,
                rank_status,
                panel.diagnostics,
            )
            fit = selection.selected.fit
        rank_row["replication_runtime_seconds"] = time.perf_counter() - started
        if rank_status != "success":
            return [
                rank_row,
                _failure_record(
                    task,
                    config,
                    rank_status,
                    "selected rank equals imposed cap",
                    {"dgp_realization_hash": panel.diagnostics["dgp_realization_hash"]},
                ),
            ]

        if not fit.converged:
            return [
                rank_row,
                _failure_record(
                    task,
                    config,
                    "full_fit_failure",
                    f"rank={fit.ranks}",
                    {"dgp_realization_hash": panel.diagnostics["dgp_realization_hash"]},
                ),
            ]
        if fit.stationarity_residual > float(config["estimation"]["stationarity_tol"]):
            return [
                rank_row,
                _failure_record(
                    task,
                    config,
                    "full_fit_failure",
                    str(fit.stationarity_residual),
                    {"dgp_realization_hash": panel.diagnostics["dgp_realization_hash"]},
                ),
            ]
        if fit.max_envelope_ratio >= 1.0:
            return [
                rank_row,
                _failure_record(
                    task,
                    config,
                    "coefficient_bound_hit",
                    str(fit.max_envelope_ratio),
                    {"dgp_realization_hash": panel.diagnostics["dgp_realization_hash"]},
                ),
            ]
        if tuple(numerical_rank(matrix) for matrix in fit.theta.matrices()) != fit.ranks:
            return [
                rank_row,
                _failure_record(
                    task,
                    config,
                    "full_fit_failure",
                    "fixed-rank numerical rank was not preserved",
                    {"dgp_realization_hash": panel.diagnostics["dgp_realization_hash"]},
                ),
            ]
        if (
            config["run"]["rank_mode"] == "fixed"
            and not fixed_rank_verification["objective_stability_pass"]
        ):
            return [
                rank_row,
                _failure_record(
                    task,
                    config,
                    "full_fit_failure",
                    "fixed-rank best objective was not independently reproduced",
                    {"dgp_realization_hash": panel.diagnostics["dgp_realization_hash"]},
                ),
            ]

        target_specs = [
            target_direction(str(name), panel.theta0, panel.groups, dgp=dgp)
            for name in config["inference"]["targets"]
        ]
        compute_gram = bool(config["inference"]["compute_tangent_gram"])
        full_system = (
            prepare_riesz_system(
                fit.theta,
                panel.design,
                fit.ranks,
                compute_tangent_gram=compute_gram,
                tangent_gram_tolerance=float(config["inference"]["tangent_gram_tol"]),
                tangent_gram_max_iter=int(config["inference"]["tangent_gram_max_iter"]),
            )
            if target_specs
            else None
        )
        split_started = time.perf_counter()
        split_bundle = (
            prepare_split_fits(
                fit,
                panel.y,
                panel.design,
                panel.groups,
                time_seed=seed_sequence(
                    config["run"]["master_seed"], dgp, n, t, replication, true_rank, "time_split"
                ),
                unit_seed=seed_sequence(
                    config["run"]["master_seed"], dgp, n, t, replication, true_rank, "unit_split"
                ),
                split_relative_rank_floor=float(
                    config["inference"]["split_relative_rank_floor"]
                ),
                compute_tangent_gram=compute_gram,
                tangent_gram_tolerance=float(config["inference"]["tangent_gram_tol"]),
                tangent_gram_max_iter=int(config["inference"]["tangent_gram_max_iter"]),
                fit_options=_fit_options(config),
            )
            if any(spec.broad for spec in target_specs)
            else None
        )
        split_fit_runtime = time.perf_counter() - split_started if split_bundle else 0.0
        split_coefficient_fit_count = (
            split_bundle.coefficient_fit_count if split_bundle is not None else 0
        )
        rank_row["split_coefficient_fit_count"] = split_coefficient_fit_count
        rank_row["split_fit_runtime_seconds"] = split_fit_runtime
        fit_metadata_rows = (
            [] if selection is None else _selected_fit_diagnostic_records(task, config, selection)
        )
        if split_bundle is not None:
            fit_metadata_rows.extend(
                _fit_diagnostic_record(
                    task,
                    config,
                    split_record.fit,
                    fit_type=f"{split_record.kind}_split_{split_record.part}",
                    initialization_route="full_fit_restriction",
                )
                for split_record in split_bundle.records
            )
        rows: list[dict[str, Any]] = [rank_row, *fit_metadata_rows]
        for spec in target_specs:
            name = spec.name
            truth = target_value(spec.direction, panel.theta0)
            regularity = target_regularity_diagnostics(spec, panel.theta0)
            inference_started = time.perf_counter()
            if spec.broad:
                if split_bundle is None or full_system is None:
                    raise RuntimeError("broad target requires prepared full and split systems")
                result = infer_corrected_target(
                    spec.direction,
                    fit,
                    full_system,
                    panel.y,
                    panel.design,
                    split_bundle,
                    spatial=_use_spatial_variance(config, dgp),
                    c_sp=float(config["inference"]["spatial_c"]),
                    riesz_tolerance=float(config["inference"]["riesz_tol"]),
                    riesz_max_iter=int(config["inference"]["riesz_max_iter"]),
                    riesz_solver=str(config["inference"]["riesz_solver"]),
                    target_rayleigh_floor=float(
                        config["inference"]["riesz_target_rayleigh_floor"]
                    ),
                    target_support_tolerance=float(
                        config["inference"]["target_support_tolerance"]
                    ),
                    tangent_gram_min_eigenvalue_floor=float(
                        config["inference"]["tangent_gram_min_eigenvalue_floor"]
                    ),
                )
            else:
                if full_system is None:
                    raise RuntimeError("target inference requires a prepared full Riesz system")
                result = infer_target(
                    spec.direction,
                    fit,
                    panel.y,
                    panel.design,
                    spatial=_use_spatial_variance(config, dgp),
                    c_sp=float(config["inference"]["spatial_c"]),
                    riesz_tolerance=float(config["inference"]["riesz_tol"]),
                    riesz_max_iter=int(config["inference"]["riesz_max_iter"]),
                    riesz_solver=str(config["inference"]["riesz_solver"]),
                    target_support_tolerance=float(
                        config["inference"]["target_support_tolerance"]
                    ),
                    tangent_gram_min_eigenvalue_floor=float(
                        config["inference"]["tangent_gram_min_eigenvalue_floor"]
                    ),
                    riesz_system=full_system,
                )
            status = classify_inference_status(result, config)

            plugin_estimate = float(result.diagnostics.get("plugin_estimate", result.estimate))
            plugin_se = float(result.diagnostics.get("plugin_standard_error", result.standard_error))
            plugin_variance = float(result.diagnostics.get("plugin_variance", result.variance))
            corrected_estimate = float(result.estimate)
            corrected_se = float(result.standard_error)
            corrected_variance = float(result.variance)
            primary_se = corrected_se if spec.broad else plugin_se
            primary_estimate = corrected_estimate if spec.broad else plugin_estimate
            if status == "success" and not np.isfinite(primary_estimate):
                status = "nonfinite_estimate"
            elif status == "success" and not np.isfinite(primary_se):
                status = "nonfinite_standard_error"
            elif status == "success" and (
                not np.isfinite(row_variance := (corrected_variance if spec.broad else plugin_variance))
                or row_variance <= 0.0
            ):
                status = "invalid_variance"
            valid_interval = bool(
                status == "success" and np.isfinite(primary_se) and primary_se > 0.0
            )
            z_true = abs(primary_estimate - truth) / primary_se if valid_interval else np.nan
            z_zero = abs(primary_estimate) / primary_se if valid_interval else np.nan
            best_neighbor_rank = (
                selection.diagnostics.get("best_neighbor_rank")
                if selection is not None
                else None
            )
            neighbor_change = None
            if best_neighbor_rank is not None and selection is not None:
                neighbor = selection.candidates.get(tuple(best_neighbor_rank))
                if neighbor is not None and neighbor.valid:
                    neighbor_change = target_value(spec.direction, neighbor.fit.theta) - target_value(
                        spec.direction, fit.theta
                    )
            row: dict[str, Any] = {
                "record_type": "target",
                "run_id": config_hash(config),
                "method": method_name(config["run"]["rank_mode"]),
                "primary_status": canonical_status(status),
                "warning_flags": json.dumps([]),
                "semantic_replication_id": semantic_replication_id(
                    dgp, n, t, replication, true_rank
                ),
                "dgp": dgp,
                "N": n,
                "T": t,
                "replication": replication,
                "target": name,
                "status": status,
                "true_rank_vector": json.dumps(true_rank),
                "truth": truth,
                "estimate": primary_estimate,
                "standard_error": primary_se,
                "variance": corrected_variance if spec.broad else plugin_variance,
                "plugin_estimate": plugin_estimate,
                "plugin_standard_error": plugin_se,
                "plugin_variance": plugin_variance,
                "plugin_error": plugin_estimate - truth,
                "corrected_estimate": corrected_estimate,
                "corrected_standard_error": corrected_se,
                "corrected_variance": corrected_variance,
                "corrected_error": corrected_estimate - truth,
                "phi_full": float(result.diagnostics.get("phi_full", plugin_estimate)),
                "phi_time_sum": result.diagnostics.get("phi_time_sum"),
                "phi_unit_sum": result.diagnostics.get("phi_unit_sum"),
                "split_fit_count": len(result.diagnostics.get("split_fits", [])),
                "split_coefficient_fit_count": split_coefficient_fit_count,
                "split_fit_runtime_seconds": split_fit_runtime,
                "split_diagnostics_json": json.dumps(
                    result.diagnostics.get("split_fits", []),
                    sort_keys=True,
                    default=_json_default,
                ),
                "time_split_assignments_json": json.dumps(result.diagnostics.get("time_parts")),
                "unit_split_assignments_json": json.dumps(result.diagnostics.get("unit_parts")),
                "centered_reject_5pct": bool(
                    valid_interval and z_true > 1.959963984540054
                ),
                "reject_zero_5pct": bool(valid_interval and z_zero > 1.959963984540054),
                "covered_95pct": bool(valid_interval and z_true <= 1.959963984540054),
                "corrected": result.corrected,
                "supplied_rank_vector": (
                    json.dumps(fit.ranks) if selection is None else None
                ),
                "selected_rank_vector": (
                    json.dumps(fit.ranks) if selection is not None else None
                ),
                "selected_rank": json.dumps(fit.ranks),
                "selected_ic": selection.selected.ic if selection is not None else None,
                "rank_at_cap": (
                    selection.diagnostics["selected_rank_at_cap"]
                    if selection is not None
                    else None
                ),
                "candidate_count": (
                    selection.diagnostics["candidate_count_final"]
                    if selection is not None
                    else None
                ),
                "candidate_coverage": (
                    true_rank
                    in {
                        tuple(rank)
                        for rank in selection.diagnostics["candidate_rank_vectors"]
                    }
                    if selection is not None
                    else None
                ),
                "cap_pilot_rank": (
                    json.dumps(selection.diagnostics["rank_cap_thresholded_vector"])
                    if selection is not None
                    else None
                ),
                "rank_selection_diagnostics": (
                    json.dumps(
                        selection.diagnostics,
                        sort_keys=True,
                        default=_json_default,
                    )
                    if selection is not None
                    and config["run"].get("save_candidate_details", True)
                    else None
                ),
                "ic_gap": selection.diagnostics["ic_gap"] if selection is not None else None,
                "best_neighbor_rank": json.dumps(best_neighbor_rank),
                "best_neighbor_target_change": neighbor_change,
                "full_objective": fit.objective,
                "full_iterations": fit.iterations,
                "stationarity_residual": fit.stationarity_residual,
                "max_envelope_ratio": fit.max_envelope_ratio,
                "riesz_iterations": result.riesz.iterations,
                "riesz_equation_residual": result.riesz.equation_residual,
                "riesz_target_rayleigh_quotient": result.riesz.target_rayleigh_quotient,
                "riesz_target_tangent_norm": result.riesz.target_tangent_norm,
                "tangent_gram_coordinate_dimension": result.diagnostics.get(
                    "tangent_gram_coordinate_dimension"
                ),
                "tangent_gram_smallest_eigenvalue": result.diagnostics.get(
                    "tangent_gram_smallest_eigenvalue"
                ),
                "tangent_gram_largest_eigenvalue": result.diagnostics.get(
                    "tangent_gram_largest_eigenvalue"
                ),
                "tangent_gram_condition_number": result.diagnostics.get(
                    "tangent_gram_condition_number"
                ),
                "tangent_gram_eigensolver_converged": result.diagnostics.get(
                    "tangent_gram_eigensolver_converged"
                ),
                "tangent_gram_eigensolver_status": result.diagnostics.get(
                    "tangent_gram_eigensolver_status"
                ),
                "weighted_residual_identity": result.riesz.weighted_residual_identity,
                "coefficient_bound_active": fit.max_envelope_ratio >= 1.0,
                "spatial_cutoff": result.diagnostics["spatial_cutoff"],
                "spatial_c": float(config["inference"]["spatial_c"]),
                "rank_runtime_seconds": rank_runtime,
                "inference_runtime_seconds": time.perf_counter() - inference_started,
                "replication_runtime_seconds": time.perf_counter() - started,
                "seed_dgp": str(dgp_seed.entropy),
                "config_hash": config_hash(config),
                "git_commit": _git_commit(),
            }
            if config["run"].get("experiment") == "power":
                block = str(config["run"].get("power_block", "A"))
                default_truth = panel.truths.get(
                    f"{block}_G2_minus_G1_time_average_true",
                    panel.truths.get(f"{block}_G2_minus_G1_fixed_time_true", truth),
                )
                row.update(
                    {
                        "nominal_delta": float(config["run"]["nominal_delta"]),
                        "power_block": block,
                        "null_or_alternative": (
                            "size/null"
                            if float(config["run"]["nominal_delta"]) == 0.0
                            else "alternative"
                        ),
                        "realized_true_contrast": float(default_truth),
                    }
                )
            row.update({key: value for key, value in panel.diagnostics.items() if np.isscalar(value)})
            row.update(regularity)
            row.update(panel.truths)
            rows.append(row)
            split_statuses = {
                f"{item.get('kind')}_split_{item.get('part')}_status": (
                    "success"
                    if item.get("converged", False)
                    and item.get("rank_supported", False)
                    and item.get("target_supported", True)
                    else "failure"
                )
                for item in result.diagnostics.get("split_fits", [])
            }
            rows.append(
                {
                    "record_type": "inference_diagnostic",
                    "run_id": config_hash(config),
                    "semantic_replication_id": semantic_replication_id(
                        dgp, n, t, replication, true_rank
                    ),
                    "dgp": dgp, "N": n, "T": t, "replication": replication,
                    "method": method_name(config["run"]["rank_mode"]),
                    "target": name, "primary_status": canonical_status(status),
                    "target_tangent_norm": result.riesz.target_tangent_norm,
                    "target_supported": result.diagnostics.get("target_supported", True),
                    "riesz_iterations": result.riesz.iterations,
                    "riesz_residual": result.riesz.equation_residual,
                    "riesz_converged": result.riesz.converged,
                    "riesz_target_rayleigh_quotient": result.riesz.target_rayleigh_quotient,
                    "tangent_gram_min_eigenvalue": result.diagnostics.get("tangent_gram_smallest_eigenvalue"),
                    "tangent_gram_max_eigenvalue": result.diagnostics.get("tangent_gram_largest_eigenvalue"),
                    "tangent_gram_condition_number": result.diagnostics.get("tangent_gram_condition_number"),
                    "variance_estimate": row["variance"], "standard_error": primary_se,
                    "interval_length": 2.0 * 1.959963984540054 * primary_se
                    if np.isfinite(primary_se) else np.nan,
                    "estimate_finite": np.isfinite(primary_estimate),
                    "standard_error_finite": np.isfinite(primary_se),
                    "phi_full": row["phi_full"], "phi_time_sum": row["phi_time_sum"],
                    "phi_unit_sum": row["phi_unit_sum"], "phi_corrected": corrected_estimate,
                    "true_target_projection_ratio": regularity.get("true_target_projection_ratio"),
                    "true_entry_unit_leverage_scaled": regularity.get("true_entry_unit_leverage_scaled"),
                    "true_entry_time_leverage_scaled": regularity.get("true_entry_time_leverage_scaled"),
                    "target_applicability": regularity.get("target_applicability"),
                    "headline_theorem_target": regularity.get("headline_theorem_target"),
                    **split_statuses,
                }
            )
        return rows
    except Exception as exc:  # Every requested replication leaves an auditable record.
        return [_failure_record(task, config, "unexpected_exception", f"{type(exc).__name__}: {exc}")]


def _worker(payload: tuple[Task, dict[str, Any], dict[str, Any]]) -> list[dict[str, Any]]:
    task, config, calibration = payload
    observed_fits = []
    observed_nuclear_fits = []
    with threadpool_limits(limits=int(config["run"]["blas_threads"])):
        with observe_fixed_rank_fits(observed_fits.append):
            with observe_nuclear_fits(observed_nuclear_fits.append):
                rows = run_replication(task, config, calibration)
    fit_metadata = [row for row in rows if row["record_type"] == "fit_diagnostic"]
    rows = [row for row in rows if row["record_type"] != "fit_diagnostic"]
    split_count = max(
        (int(row.get("split_coefficient_fit_count", 0) or 0) for row in rows),
        default=0,
    )
    split_labels = ("time_split_0", "time_split_1", "unit_split_0", "unit_split_1")
    exact_fit_rows = []
    for index, fit in enumerate(observed_fits):
        split_offset = len(observed_fits) - split_count
        if index >= split_offset and split_count == 4:
            fit_type = split_labels[index - split_offset]
        elif config["run"]["rank_mode"] == "fixed" and index == 0:
            fit_type = "full_fixed_rank"
        else:
            fit_type = "coefficient_fit"
        exact = _fit_diagnostic_record(
                task,
                config,
                fit,
                fit_type=fit_type,
                start_number=index + 1,
            )
        matching = next(
            (
                row
                for row in fit_metadata
                if row.get("fit_type") != "nuclear_path"
                and row.get("requested_rank") == exact["requested_rank"]
                and row.get("objective_final") is not None
                and np.isclose(
                    float(row["objective_final"]),
                    float(exact["objective_final"]),
                    rtol=1e-10,
                    atol=1e-12,
                )
            ),
            None,
        )
        if matching is not None:
            for key in (
                "initialization_route", "candidate_source", "IC", "IC_valid",
                "best_start_objective", "second_start_objective",
                "objective_stability_gap", "objective_stability_pass",
            ):
                if matching.get(key) is not None:
                    exact[key] = matching[key]
            fit_metadata.remove(matching)
        exact_fit_rows.append(exact)
    nuclear_rows = [
        _nuclear_fit_diagnostic_record(task, config, fit, index)
        for index, fit in enumerate(observed_nuclear_fits)
    ]
    rows.extend([*nuclear_rows, *exact_fit_rows])
    dgp_realization_hash = next(
        (
            row.get("dgp_realization_hash")
            for row in rows
            if row.get("dgp_realization_hash") is not None
        ),
        None,
    )
    for row in rows:
        row.setdefault("dgp_realization_hash", dgp_realization_hash)
    dgp, n, t, replication, true_rank = _task_parts(task)
    failures = [row for row in rows if row["record_type"] == "failure"]
    primary_status = (
        failures[0]["primary_status"] if failures else "success"
    )
    method = method_name(config["run"]["rank_mode"])
    rows.append(
        {
            "record_type": "replication",
            "run_id": config_hash(config),
            "dgp": dgp,
            "N": n,
            "T": t,
            "replication": replication,
            "method": method,
            "target": "__replication__",
            "primary_status": primary_status,
            "status": primary_status,
            "completed_dgp_replication": primary_status != "calibration_failure",
            "point_estimate_valid": False,
            "inference_valid": False,
            "retained_for_bias_rmse": False,
            "retained_for_coverage": False,
            "retained_for_rejection": False,
            "supplied_rank_vector": (
                json.dumps(true_rank if task[4] is not None else config["estimation"]["fixed_ranks"])
                if method == "fixed_rank"
                else None
            ),
            "selected_rank_vector": next(
                (row.get("selected_rank_vector") for row in rows if row["record_type"] == "rank"),
                None,
            ),
            "cap_pilot_rank": next(
                (row.get("cap_pilot_rank") for row in rows if row["record_type"] == "rank"),
                None,
            ),
            "candidate_coverage": next(
                (row.get("candidate_coverage") for row in rows if row["record_type"] == "rank"),
                None,
            ),
            "rank_selection_diagnostics": next(
                (row.get("rank_selection_diagnostics") for row in rows if row["record_type"] == "rank"),
                None,
            ),
            "expected_fit_count": sum(row["record_type"] == "fit_diagnostic" for row in rows),
            "semantic_replication_id": semantic_replication_id(
                dgp, n, t, replication, true_rank
            ),
            "dgp_realization_hash": dgp_realization_hash,
            "replication_runtime_seconds": max(
                (float(row.get("replication_runtime_seconds", 0.0) or 0.0) for row in rows),
                default=0.0,
            ),
            "exception_type": (
                failures[0].get("failure_detail", "").split(":", 1)[0]
                if primary_status == "software_exception"
                else None
            ),
            "exception_message": failures[0].get("failure_detail") if failures else None,
        }
    )
    return rows


def _calibration_failure_task_rows(
    task: Task, config: dict[str, Any], detail: str
) -> list[dict[str, Any]]:
    failure = _failure_record(task, config, "calibration_failure", detail)
    dgp, n, t, replication, true_rank = _task_parts(task)
    method = method_name(config["run"]["rank_mode"])
    attempt = {
        "record_type": "replication", "run_id": config_hash(config),
        "dgp": dgp, "N": n, "T": t, "replication": replication,
        "method": method, "target": "__replication__",
        "primary_status": "calibration_failure", "status": "calibration_failure",
        "completed_dgp_replication": False, "point_estimate_valid": False,
        "inference_valid": False, "retained_for_bias_rmse": False,
        "retained_for_coverage": False, "retained_for_rejection": False,
        "supplied_rank_vector": (
            json.dumps(true_rank if task[4] is not None else config["estimation"]["fixed_ranks"])
            if method == "fixed_rank" else None
        ),
        "selected_rank_vector": None, "cap_pilot_rank": None,
        "candidate_coverage": None, "rank_selection_diagnostics": None,
        "expected_fit_count": 0,
        "semantic_replication_id": semantic_replication_id(dgp, n, t, replication, true_rank),
        "replication_runtime_seconds": 0.0,
        "exception_type": detail.split(":", 1)[0], "exception_message": detail,
    }
    return [failure, attempt]


def _atomic_parquet(rows: list[dict[str, Any]], destination: Path) -> None:
    destination.parent.mkdir(parents=True, exist_ok=True)
    temporary = destination.with_suffix(destination.suffix + ".partial")
    frame = pd.DataFrame(rows) if rows else pd.DataFrame({"record_type": pd.Series(dtype="object")})
    frame.to_parquet(temporary, index=False)
    os.replace(temporary, destination)


def _consolidate_chunks(root: Path, subdirectory: str, stem: str, columns: tuple[str, ...]) -> None:
    files = sorted((root / subdirectory).glob("*.parquet"))
    frames = [pd.read_parquet(path) for path in files]
    usable = [frame for frame in frames if len(frame.columns)]
    combined = pd.concat(usable, ignore_index=True) if usable else pd.DataFrame(columns=columns)
    combined.to_parquet(root / f"{stem}.parquet", index=False)
    combined.to_csv(root / f"{stem}.csv", index=False)


def _task_designs(config: dict[str, Any]) -> list[tuple[int, ...] | None]:
    if "rank_stress" not in config:
        return [None]
    return [tuple(int(rank) for rank in ranks) for ranks in config["rank_stress"]["true_rank_vectors"]]


def run_monte_carlo(
    config: dict[str, Any],
    *,
    resume: bool = False,
    overwrite: bool = False,
    n_jobs: int | None = None,
    cli_argv: list[str] | None = None,
) -> tuple[dict[str, Any], Path]:
    group_gap_failure = None
    try:
        resolved, gap_reports = resolve_group_gap(config)
    except Exception as exc:
        resolved = deepcopy(config)
        gap_reports = []
        group_gap_failure = f"{type(exc).__name__}: {exc}"
    if resolved["dgp"].get("calibration_seed") is None:
        resolved["dgp"]["calibration_seed"] = int(resolved["run"]["master_seed"])
    digest = config_hash(resolved)
    root = Path(resolved["run"]["output_root"]) / resolved["run"]["name"] / digest
    if root.exists() and not resume and not overwrite:
        raise FileExistsError(f"output exists; use --resume or --overwrite: {root}")
    root.mkdir(parents=True, exist_ok=True)
    from .cli import resolved_config_text

    write_resolved_config(resolved, root / "resolved_config.json")
    (root / "resolved_config.toml").write_text(
        resolved_config_text(resolved), encoding="utf-8"
    )
    (root / "command.txt").write_text(
        " ".join(cli_argv or []) + "\n", encoding="utf-8"
    )
    (root / "group_gap_pilot.json").write_text(
        json.dumps(gap_reports, indent=2) + "\n", encoding="utf-8"
    )
    calibration_failures: dict[
        tuple[int, int, int, tuple[int, ...] | None], str
    ] = {}
    if group_gap_failure is None:
        calibrations = calibrate_design(resolved, failures=calibration_failures)
    else:
        calibrations = {}
        for dgp in resolved["run"]["dgps"]:
            for n, t in resolved["run"]["cells"]:
                for true_rank in _task_designs(resolved):
                    calibration_failures[(int(dgp), int(n), int(t), true_rank)] = (
                        group_gap_failure
                    )
    calibration_rows = []
    for key, result in calibrations.items():
        row = asdict(result)
        row.update(row.pop("diagnostics"))
        row["true_rank_vector"] = json.dumps(key[3])
        calibration_rows.append(row)
    for key, detail in calibration_failures.items():
        calibration_rows.append(
            {
                "dgp": key[0], "n": key[1], "t": key[2],
                "true_rank_vector": json.dumps(key[3]),
                "primary_status": "calibration_failure", "failure_detail": detail,
            }
        )
    coefficient_bound = float(resolved["estimation"]["coefficient_bound"])
    required_margin = float(resolved["estimation"]["simulation_interior_margin"])
    envelopes = [
        float(row["theoretical_coefficient_envelope"])
        for row in calibration_rows
        if row.get("theoretical_coefficient_envelope") is not None
    ]
    maximum_envelope = max(envelopes) if envelopes else None
    if maximum_envelope is not None and maximum_envelope > coefficient_bound - required_margin:
        raise ValueError(
            "deterministic coefficient envelope violates simulation interior condition: "
            f"max envelope={maximum_envelope:.6f}, B={coefficient_bound:.6f}, "
            f"required margin={required_margin:.6f}"
        )
    pd.DataFrame(calibration_rows).to_parquet(root / "calibration.parquet", index=False)

    replications = int(resolved["run"]["replications"])
    chunk_size = int(resolved["run"]["chunk_size"])
    # A command-line worker override is an execution detail, not part of the
    # resolved statistical design or its content-addressed run identity.
    jobs = int(resolved["run"]["n_jobs"] if n_jobs is None else n_jobs)
    for true_rank in _task_designs(resolved):
        rank_suffix = "" if true_rank is None else "_true" + "-".join(map(str, true_rank))
        for dgp in resolved["run"]["dgps"]:
            for n, t in resolved["run"]["cells"]:
                calibration_key = (int(dgp), int(n), int(t), true_rank)
                calibration = (
                    asdict(calibrations[calibration_key])
                    if calibration_key in calibrations
                    else None
                )
                for begin in range(0, replications, chunk_size):
                    end = min(replications, begin + chunk_size)
                    stem = f"dgp{dgp}_N{n}_T{t}{rank_suffix}_r{begin:05d}-{end - 1:05d}.parquet"
                    target_destination = root / "raw" / stem
                    rank_destination = root / "rank" / stem
                    fit_destination = root / "fit" / stem
                    inference_destination = root / "inference" / stem
                    replication_destination = root / "replications" / stem
                    destinations = (
                        target_destination,
                        rank_destination,
                        fit_destination,
                        inference_destination,
                        replication_destination,
                    )
                    if all(path.exists() for path in destinations) and resume:
                        continue
                    tasks: list[Task] = [
                        (int(dgp), int(n), int(t), replication, true_rank)
                        for replication in range(begin, end)
                    ]
                    payloads = [(task, resolved, calibration) for task in tasks]
                    if calibration is None:
                        nested = [
                            _calibration_failure_task_rows(
                                task, resolved, calibration_failures[calibration_key]
                            )
                            for task in tasks
                        ]
                    elif resolved["run"]["parallel_level"] == "replications" and jobs != 1:
                        with ProcessPoolExecutor(max_workers=None if jobs < 1 else jobs) as executor:
                            nested = list(executor.map(_worker, payloads))
                    else:
                        nested = [_worker(payload) for payload in payloads]
                    flat = [row for rows in nested for row in rows]
                    target_rows = [
                        row
                        for row in flat
                        if row["record_type"] in {"target", "failure"}
                    ]
                    rank_rows = [row for row in flat if row["record_type"] == "rank"]
                    fit_rows = [row for row in flat if row["record_type"] == "fit_diagnostic"]
                    inference_rows = [row for row in flat if row["record_type"] == "inference_diagnostic"]
                    replication_rows = [row for row in flat if row["record_type"] == "replication"]
                    _atomic_parquet(target_rows, target_destination)
                    _atomic_parquet(rank_rows, rank_destination)
                    _atomic_parquet(fit_rows, fit_destination)
                    _atomic_parquet(inference_rows, inference_destination)
                    _atomic_parquet(replication_rows, replication_destination)
    manifest = {
        "config_hash": digest,
        "git_commit": _git_commit(),
        "cli_argv": cli_argv or [],
        "rank_mode": resolved["run"]["rank_mode"],
        "method": method_name(resolved["run"]["rank_mode"]),
        "failure_codes": FAILURE_CODES,
        "frozen_group_means_a": [
            resolved["dgp"]["mu_lambda_a_1"],
            resolved["dgp"]["mu_lambda_a_2"],
        ],
        "requested_replications_per_cell": replications,
        "true_rank_designs": _task_designs(resolved),
        "coefficient_bound_B": coefficient_bound,
        "required_simulation_interior_margin": required_margin,
        "maximum_deterministic_coefficient_envelope": maximum_envelope,
        "deterministic_interior_margin": (
            coefficient_bound - maximum_envelope
            if maximum_envelope is not None
            else None
        ),
    }
    (root / "run_manifest.json").write_text(
        json.dumps(manifest, indent=2) + "\n", encoding="utf-8"
    )
    (root / "manifest.json").write_text(
        json.dumps(manifest, indent=2) + "\n", encoding="utf-8"
    )
    _consolidate_chunks(root, "replications", "attempted_replications", REPLICATION_COLUMNS)
    _consolidate_chunks(root, "fit", "fit_diagnostics", FIT_DIAGNOSTIC_COLUMNS)
    _consolidate_chunks(root, "inference", "inference_diagnostics", INFERENCE_DIAGNOSTIC_COLUMNS)
    attempts = pd.read_parquet(root / "attempted_replications.parquet")
    raw_files = sorted((root / "raw").glob("*.parquet"))
    raw = pd.concat([pd.read_parquet(path) for path in raw_files], ignore_index=True)
    targets = raw.loc[raw["record_type"].eq("target")].copy()
    replication_records = build_replication_records(
        attempts,
        targets,
        targets=[str(value) for value in resolved["inference"]["targets"]],
    )
    replication_records.to_parquet(root / "replication_records.parquet", index=False)
    replication_records.to_csv(root / "replication_records.csv", index=False)
    return resolved, root
