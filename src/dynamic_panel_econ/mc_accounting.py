"""Lossless Monte Carlo schemas, retention rules, summaries, and reconciliation."""

from __future__ import annotations

import json
from collections.abc import Sequence
from copy import deepcopy
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd

PRIMARY_STATUSES = (
    "success",
    "calibration_failure",
    "full_fit_failure",
    "rank_pilot_failure",
    "rank_selection_numerically_unresolved",
    "selected_rank_post_refit_numerically_unresolved",
    "rank_selection_failure",
    "rank_at_cap",
    "candidate_numerically_unresolved",
    "boundary_interiority_failure",
    "coefficient_bound_hit",
    "constrained_solver_failure",
    "constrained_feasibility_failure",
    "constrained_optimality_failure",
    "nonfinite_constrained_solution",
    "target_unsupported_selected_rank",
    "split_fit_failure",
    "split_rank_loss",
    "split_target_support_loss",
    "tangent_gram_eigensolver_failure",
    "tangent_gram_nearly_singular",
    "riesz_solver_failure",
    "riesz_target_instability",
    "nonfinite_estimate",
    "nonfinite_standard_error",
    "invalid_variance",
    "software_exception",
)

PRIMARY_STATUS_DEFINITIONS = {
    "success": "requested target completed with valid point and inference outputs",
    "calibration_failure": "DGP calibration or generated-input validity failed",
    "full_fit_failure": "supplied or selected full-panel coefficient fit failed numerical checks",
    "rank_pilot_failure": "selected-mode rank-at-most-cap pilot was unresolved",
    "rank_selection_numerically_unresolved": (
        "Revision-10 cap+1 pilot failed maintained numerical acceptance diagnostics"
    ),
    "selected_rank_post_refit_numerically_unresolved": (
        "Revision-10 final literal selected-rank post-refit was unresolved"
    ),
    "rank_selection_failure": "selected-mode candidate comparison could not return a rank",
    "rank_at_cap": "selected rank reached an imposed cap",
    "candidate_numerically_unresolved": "no stable candidate post-refit was eligible",
    "boundary_interiority_failure": (
        "a required coefficient fit was boundary-active; the point estimate remains valid "
        "but theorem-based inference is suppressed"
    ),
    "coefficient_bound_hit": "a legacy fitted coefficient reached the old diagnostic boundary",
    "constrained_solver_failure": "the literal box-constrained numerical solver failed",
    "constrained_feasibility_failure": "a constrained fit violated the entrywise box tolerance",
    "constrained_optimality_failure": "a constrained fit failed its prespecified KKT tolerance",
    "nonfinite_constrained_solution": "the constrained objective or coefficient solution was nonfinite",
    "target_unsupported_selected_rank": "full fitted tangent space did not support the target",
    "split_fit_failure": "one of exactly four split coefficient fits failed",
    "split_rank_loss": "a split fit lost a supplied numerical rank",
    "split_target_support_loss": "a split tangent space did not support the target",
    "tangent_gram_eigensolver_failure": "tangent-Gram spectrum computation failed",
    "tangent_gram_nearly_singular": "tangent-Gram conditioning violated its floor",
    "riesz_solver_failure": "the target-specific Riesz equation did not converge",
    "riesz_target_instability": "the target-specific Rayleigh diagnostic violated its floor",
    "nonfinite_estimate": "the primary point estimate was nonfinite",
    "nonfinite_standard_error": "the primary standard error was nonfinite",
    "invalid_variance": "the variance was nonfinite or nonpositive",
    "software_exception": "an otherwise unclassified software exception was recorded",
}

STATUS_ALIASES = {
    "nonfinite_input": "calibration_failure",
    "full_fit_not_converged": "full_fit_failure",
    "first_order_residual_high": "full_fit_failure",
    "coefficient_bound_active": "boundary_interiority_failure",
    "split_target_unsupported_selected_rank": "split_target_support_loss",
    "split_tangent_gram_eigensolver_failure": "tangent_gram_eigensolver_failure",
    "split_tangent_gram_nearly_singular": "tangent_gram_nearly_singular",
    "riesz_not_converged": "riesz_solver_failure",
    "split_fit_not_converged": "split_fit_failure",
    "split_riesz_target_instability": "riesz_target_instability",
    "spatial_variance_failure": "invalid_variance",
    "nonpositive_variance": "invalid_variance",
    "unexpected_exception": "software_exception",
}

REPLICATION_COLUMNS = (
    "run_id", "semantic_replication_id", "dgp_realization_hash", "dgp", "N", "T", "replication", "method", "target",
    "primary_status", "completed_dgp_replication", "point_estimate_valid",
    "inference_valid", "retained_for_bias_rmse", "retained_for_coverage",
    "retained_for_rejection", "supplied_rank_vector", "selected_rank_vector",
    "cap_pilot_rank", "candidate_coverage", "rank_selection_diagnostics",
    "estimate", "truth", "standard_error", "variance", "covered_95pct",
    "reject_zero_5pct", "extreme_estimate_flag", "replication_runtime_seconds",
    "exception_type", "exception_message", "warning_flags",
    "unconstrained_outside_box", "constrained_fallback_used", "boundary_active",
    "constrained_solver_status",
    "expected_fit_count", "nominal_delta", "realized_true_contrast",
    "null_or_alternative", "power_block",
)

FIT_DIAGNOSTIC_COLUMNS = (
    "run_id", "semantic_replication_id", "dgp_realization_hash", "dgp", "N", "T", "replication", "method", "fit_type", "target",
    "requested_rank", "numerical_rank", "initialization_route", "start_number",
    "objective_initial", "objective_final", "relative_objective_change", "iterations",
    "convergence_flag", "iteration_cap_hit", "stationarity_residual",
    "stationarity_pass", "coefficient_envelope", "coefficient_envelope_ratio",
    "max_abs_coefficient", "coefficient_bound_hit", "sigma_1", "sigma_r",
    "sigma_r_over_sigma_1",
    "unconstrained_max_abs", "unconstrained_inside_box", "unconstrained_outside_box",
    "constrained_fallback_used", "boundary_active", "max_constraint_violation",
    "constrained_KKT_residual", "constrained_iterations", "constrained_runtime",
    "constrained_runtime_seconds",
    "constrained_solver_status", "constrained_objective",
    "best_start_objective", "second_start_objective", "objective_stability_gap",
    "objective_stability_pass", "runtime_seconds", "exception_type",
    "exception_message", "nuclear_path_index", "lambda", "thresholded_rank",
    "candidate_source", "IC", "IC_valid", "diagnostic_context",
    "initial_coefficient_envelope", "final_coefficient_envelope",
    "coefficient_envelope_history",
)

INFERENCE_DIAGNOSTIC_COLUMNS = (
    "run_id", "semantic_replication_id", "dgp_realization_hash", "dgp", "N", "T", "replication", "method", "target",
    "primary_status", "target_tangent_norm", "target_supported", "riesz_iterations",
    "riesz_residual", "riesz_converged", "riesz_target_rayleigh_quotient",
    "tangent_gram_min_eigenvalue", "tangent_gram_max_eigenvalue",
    "tangent_gram_condition_number", "variance_estimate", "standard_error",
    "interval_length", "estimate_finite", "standard_error_finite", "phi_full",
    "phi_time_sum", "phi_unit_sum", "phi_corrected", "time_split_0_status",
    "time_split_1_status", "unit_split_0_status", "unit_split_1_status",
    "true_target_projection_ratio", "true_entry_unit_leverage_scaled",
    "true_entry_time_leverage_scaled", "target_applicability", "headline_theorem_target",
)


def canonical_status(status: str) -> str:
    value = STATUS_ALIASES.get(status, status)
    if value not in PRIMARY_STATUSES:
        return "software_exception"
    return value


def method_name(rank_mode: str) -> str:
    return "fixed_rank" if rank_mode == "fixed" else "selected_rank"


def semantic_replication_id(dgp: int, n: int, t: int, replication: int, true_rank: Sequence[int]) -> str:
    """ID shared by fixed and selected methods for matched DGP draws."""

    ranks = "-".join(str(int(value)) for value in true_rank)
    return f"dgp{dgp}_N{n}_T{t}_r{replication:05d}_truth{ranks}"


def empty_schema(columns: Sequence[str]) -> pd.DataFrame:
    return pd.DataFrame({column: pd.Series(dtype="object") for column in columns})


def write_schema_files(root: str | Path) -> list[Path]:
    root = Path(root)
    root.mkdir(parents=True, exist_ok=True)
    paths = []
    for name, columns in (
        ("replication_records", REPLICATION_COLUMNS),
        ("fit_diagnostics", FIT_DIAGNOSTIC_COLUMNS),
        ("inference_diagnostics", INFERENCE_DIAGNOSTIC_COLUMNS),
    ):
        destination = root / "schemas" / f"{name}.json"
        destination.parent.mkdir(parents=True, exist_ok=True)
        destination.write_text(json.dumps({"columns": list(columns)}, indent=2) + "\n", encoding="utf-8")
        paths.append(destination)
    status_path = root / "schemas" / "primary_statuses.json"
    status_path.write_text(
        json.dumps(
            {
                "mutually_exclusive": True,
                "statuses": list(PRIMARY_STATUSES),
                "definitions": PRIMARY_STATUS_DEFINITIONS,
                "aliases_from_legacy_records": STATUS_ALIASES,
            },
            indent=2,
        )
        + "\n",
        encoding="utf-8",
    )
    paths.append(status_path)
    return paths


def apply_retention_flags(records: pd.DataFrame) -> pd.DataFrame:
    """Flag validity without deleting, trimming, or winsorizing any attempted row."""

    result = records.reset_index(drop=True).copy()
    estimate = pd.to_numeric(
        result["estimate"]
        if "estimate" in result
        else pd.Series(np.nan, index=result.index),
        errors="coerce",
    )
    standard_error = pd.to_numeric(
        result["standard_error"]
        if "standard_error" in result
        else pd.Series(np.nan, index=result.index),
        errors="coerce",
    )
    variance = pd.to_numeric(
        result["variance"]
        if "variance" in result
        else pd.Series(np.nan, index=result.index),
        errors="coerce",
    )
    point = np.isfinite(estimate)
    finite_inference = (
        point
        & np.isfinite(standard_error)
        & (standard_error > 0.0)
        & np.isfinite(variance)
        & (variance > 0.0)
    )
    status = result.get("primary_status", pd.Series("success", index=result.index)).astype(str)
    status = status.map(canonical_status)
    status = status.mask(~point & status.eq("success"), "nonfinite_estimate")
    status = status.mask(point & ~np.isfinite(standard_error) & status.eq("success"), "nonfinite_standard_error")
    status = status.mask(
        point & np.isfinite(standard_error) & (standard_error <= 0.0) & status.eq("success"),
        "invalid_variance",
    )
    status = status.mask(point & np.isfinite(standard_error) & (~np.isfinite(variance) | (variance <= 0.0)) & status.eq("success"), "invalid_variance")
    inference = finite_inference & status.eq("success")
    result["point_estimate_valid"] = point
    result["inference_valid"] = inference
    result["retained_for_bias_rmse"] = point
    result["retained_for_coverage"] = inference
    result["retained_for_rejection"] = inference
    result["primary_status"] = status
    if "warning_flags" not in result:
        result["warning_flags"] = "[]"
    if not result["inference_valid"].eq(result["primary_status"].eq("success")).all():
        raise ValueError("primary-status and inference-retention invariant failed")
    result["extreme_estimate_flag"] = False
    group_columns = [column for column in ("dgp", "N", "T", "method", "target") if column in result]
    for _, indices in result.groupby(group_columns, dropna=False, sort=False).groups.items():
        values = estimate.loc[indices]
        finite_values = values.loc[np.isfinite(values)]
        if finite_values.empty:
            continue
        median = float(finite_values.median())
        mad = float((finite_values - median).abs().median())
        if mad > 0.0:
            result.loc[finite_values.index, "extreme_estimate_flag"] = (
                (finite_values - median).abs() > 10.0 * 1.4826 * mad
            )
    return result


def _cross_attempts_targets(attempts: pd.DataFrame, targets: Sequence[str]) -> pd.DataFrame:
    identity = ["run_id", "dgp", "N", "T", "replication", "method"]
    base = attempts.drop_duplicates(identity).drop(columns="target", errors="ignore").copy()
    base["_join"] = 1
    target_frame = pd.DataFrame({"target": list(targets), "_join": 1})
    return base.merge(target_frame, on="_join", how="left", validate="many_to_many").drop(columns="_join")


def build_replication_records(
    attempts: pd.DataFrame,
    target_records: pd.DataFrame,
    *,
    targets: Sequence[str],
) -> pd.DataFrame:
    """Create one lossless accounting row per attempted replication and requested target."""

    if not targets:
        records = attempts.copy()
        records["target"] = "__replication__"
        for column in ("estimate", "truth", "standard_error", "variance"):
            records[column] = np.nan
        records["point_estimate_valid"] = False
        records["inference_valid"] = False
        records["retained_for_bias_rmse"] = False
        records["retained_for_coverage"] = False
        records["retained_for_rejection"] = False
        records["extreme_estimate_flag"] = False
        return records
    identity = ["run_id", "dgp", "N", "T", "replication", "method", "target"]
    target_records = target_records.copy()
    for column in identity:
        if column not in target_records:
            target_records[column] = pd.Series(dtype="object")
    expected = _cross_attempts_targets(attempts, targets)
    duplicate = target_records.duplicated(identity, keep=False)
    if duplicate.any():
        raise ValueError("multiple target records exist for one attempted replication")
    records = expected.merge(
        target_records,
        on=identity,
        how="left",
        validate="one_to_one",
        suffixes=("_attempt", ""),
    )
    missing = records["primary_status"].isna()
    attempt_status = records.get(
        "primary_status_attempt", pd.Series("software_exception", index=records.index)
    )
    records.loc[missing, "primary_status"] = attempt_status.loc[missing].map(canonical_status)
    return apply_retention_flags(records)


def summarize_accounting(
    attempts: pd.DataFrame,
    target_records: pd.DataFrame,
    *,
    targets: Sequence[str],
) -> pd.DataFrame:
    """Summarize against attempted denominators while retaining every failed attempt."""

    joined = build_replication_records(attempts, target_records, targets=targets)
    rows: list[dict[str, Any]] = []
    keys = ["dgp", "N", "T", "method", "target"]
    for cell, group in joined.groupby(keys, dropna=False, sort=True):
        point = group.loc[group["retained_for_bias_rmse"]]
        inference = group.loc[group["retained_for_coverage"]]
        attempted = len(group)
        errors = pd.to_numeric(point["estimate"], errors="coerce") - pd.to_numeric(point["truth"], errors="coerce")
        runtimes = pd.to_numeric(group.get("replication_runtime_seconds"), errors="coerce")
        statuses = group["primary_status"].value_counts()
        record: dict[str, Any] = {
            "dgp": cell[0], "N": cell[1], "T": cell[2], "method": cell[3], "target": cell[4],
            "attempted_replications": attempted,
            "completed_dgp_replications": int(group.get("completed_dgp_replication", False).fillna(False).sum()),
            "point_estimate_valid": int(group["point_estimate_valid"].sum()),
            "inference_valid": int(group["inference_valid"].sum()),
            "retained_for_bias_rmse": len(point),
            "retained_for_coverage": len(inference),
            "retained_for_rejection": int(group["retained_for_rejection"].sum()),
            "R_attempted": attempted, "R_point": len(point), "R_inference": len(inference),
            "success_count": int(statuses.get("success", 0)),
            "point_retained_share": len(point) / attempted,
            "inference_retained_share": len(inference) / attempted,
            "mean_truth": float(pd.to_numeric(point["truth"], errors="coerce").mean()) if len(point) else np.nan,
            "mean_estimate": float(pd.to_numeric(point["estimate"], errors="coerce").mean()) if len(point) else np.nan,
            "bias": float(errors.mean()) if len(point) else np.nan,
            "rmse": float(np.sqrt(np.mean(errors**2))) if len(point) else np.nan,
            "mc_sd": float(pd.to_numeric(point["estimate"], errors="coerce").std(ddof=1)) if len(point) > 1 else np.nan,
            "mean_se": float(pd.to_numeric(inference["standard_error"], errors="coerce").mean()) if len(inference) else np.nan,
            "coverage": float(inference["covered_95pct"].astype(bool).mean()) if len(inference) else np.nan,
            "rejection_probability": float(inference["reject_zero_5pct"].astype(bool).mean()) if len(inference) else np.nan,
            "mean_interval_length": float((2 * 1.959963984540054 * pd.to_numeric(inference["standard_error"], errors="coerce")).mean()) if len(inference) else np.nan,
            "runtime_median": float(runtimes.median()), "runtime_mean": float(runtimes.mean()),
            "runtime_p10": float(runtimes.quantile(0.10)), "runtime_p90": float(runtimes.quantile(0.90)),
            "runtime_p95": float(runtimes.quantile(0.95)),
            "replication_ids_json": json.dumps(
                [
                    f"{item.run_id}:{int(item.replication)}"
                    for item in group[["run_id", "replication"]]
                    .drop_duplicates()
                    .itertuples(index=False)
                ]
            ),
        }
        for status in PRIMARY_STATUSES:
            record[f"failure_{status}"] = int(statuses.get(status, 0)) if status != "success" else 0
        support_failures = int(statuses.get("target_unsupported_selected_rank", 0) + statuses.get("split_target_support_loss", 0))
        total_failures = attempted - len(inference)
        record["target_support_failure_rate"] = support_failures / attempted
        record["total_inference_failure_rate"] = total_failures / attempted
        record["numerical_failure_rate"] = (total_failures - support_failures) / attempted
        rows.append(record)
    summary = pd.DataFrame(rows)
    reconcile_summary(summary)
    return summary


def reconcile_summary(summary: pd.DataFrame) -> None:
    for row in summary.to_dict("records"):
        attempted = int(row["R_attempted"])
        failures = sum(int(row.get(f"failure_{status}", 0)) for status in PRIMARY_STATUSES if status != "success")
        successes = int(row.get("success_count", 0))
        if attempted != successes + failures:
            raise ValueError("attempted replication accounting does not reconcile")
        if not 0 <= int(row["R_inference"]) <= int(row["R_point"]) <= attempted:
            raise ValueError("retention denominators do not reconcile")
        if int(row["R_inference"]) != successes:
            raise ValueError("successful inference count does not equal R_inference")


def reconcile_fit_rows(fit_rows: pd.DataFrame, expected_fit_counts: pd.DataFrame) -> None:
    keys = ["run_id", "dgp", "N", "T", "replication", "method"]
    actual = fit_rows.groupby(keys, dropna=False).size().rename("actual_fit_count").reset_index()
    check = expected_fit_counts.merge(actual, on=keys, how="left", validate="one_to_one")
    if not (check["expected_fit_count"] == check["actual_fit_count"].fillna(0)).all():
        raise ValueError("fit diagnostic rows do not reconcile with executed fits")


def reconcile_matched_draws(attempts: pd.DataFrame) -> None:
    """Require fixed and selected methods to use identical semantic replication IDs."""

    if "semantic_replication_id" not in attempts:
        raise ValueError("semantic_replication_id is required for method matching")
    methods = set(attempts["method"].unique())
    if not {"fixed_rank", "selected_rank"} <= methods:
        return
    keys = ["dgp", "N", "T"]
    for cell, group in attempts.groupby(keys, dropna=False, sort=True):
        fixed = set(group.loc[group["method"].eq("fixed_rank"), "semantic_replication_id"])
        selected = set(group.loc[group["method"].eq("selected_rank"), "semantic_replication_id"])
        if fixed != selected:
            raise ValueError(f"fixed and selected DGP draws do not match for cell {cell}")


def power_loading_means(delta: float) -> tuple[float, float]:
    return 1.0 - float(delta) / 2.0, 1.0 + float(delta) / 2.0


def power_design_configs(base_config: dict[str, Any]) -> list[dict[str, Any]]:
    """Expand nominal alternatives into matched fixed/selected DGP-4 designs."""

    designs = []
    block = str(base_config["run"].get("power_block", "A"))
    for delta in base_config["run"].get("alternative_grid", [0.0]):
        first, second = power_loading_means(float(delta))
        for rank_mode in ("fixed", "selected"):
            config = deepcopy(base_config)
            config["run"]["experiment"] = "power"
            config["run"]["rank_mode"] = rank_mode
            config["run"]["dgps"] = [4]
            config["run"]["nominal_delta"] = float(delta)
            if block == "A":
                config["dgp"]["mu_lambda_a_1"] = first
                config["dgp"]["mu_lambda_a_2"] = second
            else:
                config["dgp"]["mu_lambda_b_1"] = first
                config["dgp"]["mu_lambda_b_2"] = second
            designs.append(config)
    return designs


def summarize_power(records: pd.DataFrame) -> pd.DataFrame:
    valid = apply_retention_flags(records)
    rows = []
    keys = ["method", "target", "nominal_delta"]
    for cell, group in valid.groupby(keys, dropna=False, sort=True):
        inference = group.loc[group["retained_for_rejection"]]
        denominator = len(inference)
        rejections = int(inference["reject_zero_5pct"].astype(bool).sum())
        probability = rejections / denominator if denominator else np.nan
        rows.append(
            {
                "method": cell[0], "target": cell[1], "nominal_delta": cell[2],
                "null_or_alternative": "size/null" if float(cell[2]) == 0.0 else "alternative",
                "realized_true_contrast": float(pd.to_numeric(group["realized_true_contrast"], errors="coerce").mean()),
                "attempted_replications": len(group), "valid_inference_replications": denominator,
                "rejection_count": rejections, "rejection_probability": probability,
                "mc_binomial_se": float(np.sqrt(probability * (1.0 - probability) / denominator)) if denominator else np.nan,
            }
        )
    return pd.DataFrame(rows)
