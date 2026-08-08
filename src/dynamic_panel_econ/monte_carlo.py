"""Deterministic, parallel, resume-safe Monte Carlo orchestration."""

from __future__ import annotations

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
from .inference import (
    infer_corrected_target,
    infer_target,
    prepare_riesz_system,
    prepare_split_fits,
)
from .rank_selection import RankPilotFailure, RankSelectionFailure, select_ranks
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
) -> dict[tuple[int, int, int, tuple[int, ...] | None], CalibrationResult]:
    params = _params(config)
    results = {}
    for dgp in config["run"]["dgps"]:
        for n, t in config["run"]["cells"]:
            for true_rank in _task_designs(config):
                common = {
                    "params": params,
                    "pi_h": float(config["dgp"]["pi_h"]),
                    "target_r2": float(config["dgp"]["target_r2"]),
                    "draws": int(config["dgp"]["calibration_draws"]),
                }
                if true_rank is None:
                    result = calibrate_cell(
                        int(dgp),
                        int(n),
                        int(t),
                        int(config["run"]["master_seed"]),
                        **common,
                    )
                else:
                    result = calibrate_rank_stress_cell(
                        int(dgp),
                        int(n),
                        int(t),
                        true_rank,
                        int(config["run"]["master_seed"]),
                        component_strengths=tuple(config["rank_stress"]["component_strengths"]),
                        **common,
                    )
                results[(int(dgp), int(n), int(t), true_rank)] = result
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


def _failure_record(task: Task, config: dict[str, Any], status: str, detail: str) -> dict[str, Any]:
    dgp, n, t, replication, true_rank = _task_parts(task)
    return {
        "record_type": "failure",
        "dgp": dgp,
        "N": n,
        "T": t,
        "replication": replication,
        "target": "__replication__",
        "true_rank_vector": json.dumps(true_rank),
        "status": status,
        "failure_detail": detail,
        "config_hash": config_hash(config),
        "git_commit": _git_commit(),
    }


def _rank_record(
    task: Task,
    config: dict[str, Any],
    selection: Any,
    rank_runtime: float,
    status: str,
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
        "dgp": dgp,
        "N": n,
        "T": t,
        "replication": replication,
        "status": status,
        "true_rank_vector": json.dumps(true_rank),
        "selected_rank_vector": json.dumps(selected),
        "rank_cap_thresholded_vector": json.dumps(diagnostics["rank_cap_thresholded_vector"]),
        "true_rank_in_candidates": true_rank in candidates,
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
        "cap_pilot_start_attempts": json.dumps(
            cap_pilot["outer_start_attempts"], sort_keys=True, default=_json_default
        ),
        "rank_runtime_seconds": rank_runtime,
        "rank_diagnostics_json": json.dumps(diagnostics, sort_keys=True, default=_json_default),
        "config_hash": config_hash(config),
        "git_commit": _git_commit(),
    }
    for multiplier, result in sensitivity.get("ic_multipliers", {}).items():
        record[f"ic_multiplier_{multiplier}_selected_rank"] = json.dumps(
            result["selected_rank"]
        )
    for multiplier, result in sensitivity.get("threshold_multipliers", {}).items():
        record[f"threshold_multiplier_{multiplier}_selected_rank"] = json.dumps(
            result["selected_rank"]
        )
    for index, block in enumerate(("A", "B", "H")):
        record[f"{block}_underselected"] = selected[index] < true_rank[index]
        record[f"{block}_overselected"] = selected[index] > true_rank[index]
        record[f"{block}_true_rank"] = true_rank[index]
        record[f"{block}_selected_rank"] = selected[index]
    return record


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

        selection_started = time.perf_counter()
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
            return [_failure_record(task, config, "rank_pilot_failure", str(exc))]
        except RankSelectionFailure as exc:
            return [_failure_record(task, config, "rank_selection_failure", str(exc))]
        rank_runtime = time.perf_counter() - selection_started
        rank_status = "rank_at_cap" if selection.diagnostics["selected_rank_at_cap"] else "success"
        rank_row = _rank_record(task, config, selection, rank_runtime, rank_status)
        if rank_status != "success":
            return [rank_row, _failure_record(task, config, rank_status, "selected rank equals imposed cap")]

        fit = selection.selected.fit
        if not fit.converged:
            return [rank_row, _failure_record(task, config, "full_fit_not_converged", f"rank={fit.ranks}")]
        if fit.stationarity_residual > float(config["estimation"]["stationarity_tol"]):
            return [rank_row, _failure_record(task, config, "first_order_residual_high", str(fit.stationarity_residual))]
        if fit.max_envelope_ratio >= 1.0:
            return [rank_row, _failure_record(task, config, "coefficient_bound_active", str(fit.max_envelope_ratio))]

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
        rows: list[dict[str, Any]] = [rank_row]
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
                    spatial=dgp >= 2,
                    c_sp=float(config["inference"]["spatial_c"]),
                    riesz_tolerance=float(config["inference"]["riesz_tol"]),
                    riesz_max_iter=int(config["inference"]["riesz_max_iter"]),
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
                    spatial=dgp >= 2,
                    c_sp=float(config["inference"]["spatial_c"]),
                    riesz_tolerance=float(config["inference"]["riesz_tol"]),
                    riesz_max_iter=int(config["inference"]["riesz_max_iter"]),
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
            valid_interval = bool(
                status == "success" and np.isfinite(primary_se) and primary_se > 0.0
            )
            z_true = abs(primary_estimate - truth) / primary_se if valid_interval else np.nan
            z_zero = abs(primary_estimate) / primary_se if valid_interval else np.nan
            best_neighbor_rank = selection.diagnostics.get("best_neighbor_rank")
            neighbor_change = None
            if best_neighbor_rank is not None:
                neighbor = selection.candidates.get(tuple(best_neighbor_rank))
                if neighbor is not None and neighbor.valid:
                    neighbor_change = target_value(spec.direction, neighbor.fit.theta) - target_value(
                        spec.direction, fit.theta
                    )
            row: dict[str, Any] = {
                "record_type": "target",
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
                "selected_rank": json.dumps(fit.ranks),
                "selected_ic": selection.selected.ic,
                "rank_at_cap": selection.diagnostics["selected_rank_at_cap"],
                "candidate_count": selection.diagnostics["candidate_count_final"],
                "ic_gap": selection.diagnostics["ic_gap"],
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
            row.update({key: value for key, value in panel.diagnostics.items() if np.isscalar(value)})
            row.update(regularity)
            row.update(panel.truths)
            rows.append(row)
        return rows
    except Exception as exc:  # Every requested replication leaves an auditable record.
        return [_failure_record(task, config, "unexpected_exception", f"{type(exc).__name__}: {exc}")]


def _worker(payload: tuple[Task, dict[str, Any], dict[str, Any]]) -> list[dict[str, Any]]:
    task, config, calibration = payload
    with threadpool_limits(limits=int(config["run"]["blas_threads"])):
        return run_replication(task, config, calibration)


def _atomic_parquet(rows: list[dict[str, Any]], destination: Path) -> None:
    destination.parent.mkdir(parents=True, exist_ok=True)
    temporary = destination.with_suffix(destination.suffix + ".partial")
    pd.DataFrame(rows).to_parquet(temporary, index=False)
    os.replace(temporary, destination)


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
) -> tuple[dict[str, Any], Path]:
    resolved, gap_reports = resolve_group_gap(config)
    digest = config_hash(resolved)
    root = Path(resolved["run"]["output_root"]) / resolved["run"]["name"] / digest
    if root.exists() and not resume and not overwrite:
        raise FileExistsError(f"output exists; use --resume or --overwrite: {root}")
    root.mkdir(parents=True, exist_ok=True)
    write_resolved_config(resolved, root / "resolved_config.json")
    (root / "group_gap_pilot.json").write_text(
        json.dumps(gap_reports, indent=2) + "\n", encoding="utf-8"
    )
    calibrations = calibrate_design(resolved)
    calibration_rows = []
    for key, result in calibrations.items():
        row = asdict(result)
        row.update(row.pop("diagnostics"))
        row["true_rank_vector"] = json.dumps(key[3])
        calibration_rows.append(row)
    coefficient_bound = float(resolved["estimation"]["coefficient_bound"])
    required_margin = float(resolved["estimation"]["simulation_interior_margin"])
    maximum_envelope = max(
        float(row["theoretical_coefficient_envelope"]) for row in calibration_rows
    )
    if maximum_envelope > coefficient_bound - required_margin:
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
                calibration = asdict(calibrations[(int(dgp), int(n), int(t), true_rank)])
                for begin in range(0, replications, chunk_size):
                    end = min(replications, begin + chunk_size)
                    stem = f"dgp{dgp}_N{n}_T{t}{rank_suffix}_r{begin:05d}-{end - 1:05d}.parquet"
                    target_destination = root / "raw" / stem
                    rank_destination = root / "rank" / stem
                    if target_destination.exists() and rank_destination.exists() and resume:
                        continue
                    tasks: list[Task] = [
                        (int(dgp), int(n), int(t), replication, true_rank)
                        for replication in range(begin, end)
                    ]
                    payloads = [(task, resolved, calibration) for task in tasks]
                    if resolved["run"]["parallel_level"] == "replications" and jobs != 1:
                        with ProcessPoolExecutor(max_workers=None if jobs < 1 else jobs) as executor:
                            nested = list(executor.map(_worker, payloads))
                    else:
                        nested = [_worker(payload) for payload in payloads]
                    flat = [row for rows in nested for row in rows]
                    target_rows = [row for row in flat if row["record_type"] != "rank"]
                    rank_rows = [row for row in flat if row["record_type"] == "rank"]
                    _atomic_parquet(target_rows, target_destination)
                    _atomic_parquet(rank_rows, rank_destination)
    manifest = {
        "config_hash": digest,
        "git_commit": _git_commit(),
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
        "deterministic_interior_margin": coefficient_bound - maximum_envelope,
    }
    (root / "run_manifest.json").write_text(
        json.dumps(manifest, indent=2) + "\n", encoding="utf-8"
    )
    return resolved, root
