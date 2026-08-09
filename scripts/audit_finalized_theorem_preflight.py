"""Audit the one authorized finalized theorem-aligned statistical preflight."""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd

EXPECTED_REPLICATIONS = (2026080801, 2026080802, 2026080803)
EXPECTED_DGPS = (1, 2, 3, 4)
EXPECTED_SIZES = (50, 100)
EXPECTED_TARGETS = 18
FROZEN_HASH = "e8983cadc4fbca990feeba6363420542a99ee056cf445e867532e2a6ea0e7d62"
CONSTRAINED_KKT_TOLERANCE = 1e-4
CSV_OMITTED_NESTED_DIAGNOSTICS = (
    "rank_selection_diagnostics",
    "rank_selection_diagnostics_attempt",
)

DRY_COMMANDS = (
    r".\.venv\Scripts\python.exe scripts\run_mc.py --config "
    r"configs\mc\preflight_finalized_fixed.toml --dry-run --print-resolved-config",
    r".\.venv\Scripts\python.exe scripts\run_mc.py --config "
    r"configs\mc\preflight_finalized_selected.toml --dry-run --print-resolved-config",
)
RUN_COMMANDS = (
    r".\.venv\Scripts\python.exe scripts\run_mc.py --config "
    r"configs\mc\preflight_finalized_fixed.toml",
    r".\.venv\Scripts\python.exe scripts\run_mc.py --config "
    r"configs\mc\preflight_finalized_selected.toml",
)


def _bool(values: pd.Series) -> pd.Series:
    return values.fillna(False).astype(bool)


def _num(values: pd.Series) -> pd.Series:
    return pd.to_numeric(values, errors="coerce")


def _finite(values: pd.Series) -> pd.Series:
    converted = _num(values).replace([np.inf, -np.inf], np.nan)
    return converted.dropna()


def _rate(values: pd.Series) -> float:
    return float(_bool(values).mean()) if len(values) else float("nan")


def _q(values: pd.Series, quantile: float) -> float:
    finite = _finite(values)
    return float(finite.quantile(quantile)) if len(finite) else float("nan")


def _safe_max(values: pd.Series) -> float:
    finite = _finite(values)
    return float(finite.max()) if len(finite) else float("nan")


def _read_chunks(root: Path, subdirectory: str) -> pd.DataFrame:
    paths = sorted((root / subdirectory).glob("*.parquet"))
    if not paths:
        return pd.DataFrame()
    return pd.concat([pd.read_parquet(path) for path in paths], ignore_index=True)


def _read_method(root: Path) -> dict[str, pd.DataFrame]:
    return {
        "attempts": pd.read_parquet(root / "attempted_replications.parquet"),
        "records": pd.read_parquet(root / "replication_records.parquet"),
        "fits": pd.read_parquet(root / "fit_diagnostics.parquet"),
        "inference": pd.read_parquet(root / "inference_diagnostics.parquet"),
        "rank": _read_chunks(root, "rank"),
    }


def _write_both(frame: pd.DataFrame, root: Path, stem: str) -> None:
    frame.to_parquet(root / f"{stem}.parquet", index=False)
    frame.to_csv(root / f"{stem}.csv", index=False, lineterminator="\n")


def _write_replication_records(frame: pd.DataFrame, root: Path) -> None:
    frame.to_parquet(root / "replication_records.parquet", index=False)
    csv_columns = [
        column
        for column in frame.columns
        if column not in CSV_OMITTED_NESTED_DIAGNOSTICS
    ]
    frame[csv_columns].to_csv(
        root / "replication_records.csv", index=False, lineterminator="\n"
    )


def _normalize_fit_rows(fits: pd.DataFrame) -> pd.DataFrame:
    result = fits.copy()
    fixed_generic = result["method"].eq("fixed_rank") & result["fit_type"].eq(
        "coefficient_fit"
    )
    result.loc[fixed_generic, "fit_type"] = "full_fixed_rank"
    if "max_abs_coefficient" not in result:
        result["max_abs_coefficient"] = result["coefficient_envelope"]
    else:
        result["max_abs_coefficient"] = result["max_abs_coefficient"].fillna(
            result["coefficient_envelope"]
        )
    if "constrained_runtime_seconds" not in result:
        result["constrained_runtime_seconds"] = result["constrained_runtime"]
    else:
        result["constrained_runtime_seconds"] = result[
            "constrained_runtime_seconds"
        ].fillna(result["constrained_runtime"])
    split_labels = ("time_split_0", "time_split_1", "unit_split_0", "unit_split_1")
    selected = result.loc[result["method"].eq("selected_rank")]
    for _, group in selected.groupby("semantic_replication_id", sort=False):
        observed = group.loc[group["start_number"].notna()].sort_values("start_number")
        final_four = observed.tail(4)
        if len(final_four) != 4 or final_four["start_number"].nunique() != 4:
            raise ValueError("selected split-fit diagnostic positions did not reconcile")
        result.loc[final_four.index, "fit_type"] = split_labels
    return result


def _fit_role(value: str) -> str:
    if value == "full_fixed_rank":
        return "full_panel"
    if value == "candidate_post_refit":
        return "candidate_post_refit"
    if value == "rank_cap_pilot":
        return "cap_pilot"
    if value.startswith("time_split_"):
        return "time_half_split"
    if value.startswith("unit_split_"):
        return "unit_half_split"
    if value == "nuclear_path":
        return "nuclear_screening"
    return "other_coefficient_fit"


def _selected_full_panel_rows(
    fits: pd.DataFrame, attempts: pd.DataFrame
) -> pd.DataFrame:
    candidates = fits.loc[
        fits["method"].eq("selected_rank")
        & fits["fit_type"].eq("candidate_post_refit")
    ].copy()
    rows = []
    for attempt in attempts.loc[attempts["method"].eq("selected_rank")].itertuples():
        subset = candidates.loc[
            candidates["semantic_replication_id"].eq(attempt.semantic_replication_id)
            & candidates["requested_rank"].eq(attempt.selected_rank_vector)
        ].copy()
        if "IC_valid" in subset and _bool(subset["IC_valid"]).any():
            subset = subset.loc[_bool(subset["IC_valid"])]
        objective = _num(subset.get("objective_final", pd.Series(dtype=float)))
        subset = subset.loc[objective.notna()].copy()
        if len(subset):
            chosen = subset.loc[_num(subset["objective_final"]).idxmin()].copy()
            chosen["fit_role"] = "full_panel"
            chosen["fit_role_note"] = "selected candidate serving as final full-panel estimate"
            rows.append(chosen)
    return pd.DataFrame(rows, columns=[*fits.columns, "fit_role", "fit_role_note"])


def _boundary_and_optimization(
    fits: pd.DataFrame, attempts: pd.DataFrame
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    coefficient = fits.loc[fits["fit_type"].ne("nuclear_path")].copy()
    coefficient["fit_role"] = coefficient["fit_type"].map(_fit_role)
    coefficient["fit_role_note"] = "executed coefficient fit"
    selected_full = _selected_full_panel_rows(fits, attempts)
    roles = pd.concat([coefficient, selected_full], ignore_index=True, sort=False)

    status = roles["constrained_solver_status"].fillna("").astype(str)
    kkt = _num(roles["constrained_KKT_residual"])
    fallback = _bool(roles["constrained_fallback_used"])
    roles["constrained_solver_failure"] = status.eq("constrained_solver_failure")
    roles["constrained_feasibility_failure"] = status.eq(
        "constrained_feasibility_failure"
    )
    roles["constrained_optimality_failure"] = status.eq(
        "constrained_optimality_failure"
    )
    roles["constrained_kkt_failure"] = roles["constrained_optimality_failure"] | (
        fallback & kkt.notna() & kkt.gt(CONSTRAINED_KKT_TOLERANCE)
    )
    roles["numerical_rank_loss"] = (
        roles["requested_rank"].notna()
        & roles["numerical_rank"].notna()
        & roles["requested_rank"].ne(roles["numerical_rank"])
    )

    keys = ["dgp", "N", "method", "fit_role"]
    boundary_rows: list[dict[str, Any]] = []
    optimization_rows: list[dict[str, Any]] = []
    for cell, group in roles.groupby(keys, dropna=False, sort=True):
        attempted = len(group)
        common = dict(zip(keys, cell, strict=True))
        boundary_rows.append(
            {
                **common,
                "attempted_fits": attempted,
                "unconstrained_outside_box_count": int(
                    _bool(group["unconstrained_outside_box"]).sum()
                ),
                "unconstrained_outside_box_rate": _rate(
                    group["unconstrained_outside_box"]
                ),
                "constrained_fallback_count": int(
                    _bool(group["constrained_fallback_used"]).sum()
                ),
                "constrained_fallback_rate": _rate(group["constrained_fallback_used"]),
                "boundary_active_count": int(_bool(group["boundary_active"]).sum()),
                "boundary_active_rate": _rate(group["boundary_active"]),
                "constrained_solver_failure_count": int(
                    group["constrained_solver_failure"].sum()
                ),
                "constrained_solver_failure_rate": float(
                    group["constrained_solver_failure"].mean()
                ),
                "constrained_feasibility_failure_count": int(
                    group["constrained_feasibility_failure"].sum()
                ),
                "constrained_feasibility_failure_rate": float(
                    group["constrained_feasibility_failure"].mean()
                ),
                "constrained_optimality_failure_count": int(
                    group["constrained_optimality_failure"].sum()
                ),
                "constrained_optimality_failure_rate": float(
                    group["constrained_optimality_failure"].mean()
                ),
            }
        )
        optimization_rows.append(
            {
                **common,
                "fits_attempted": attempted,
                "convergence_rate": _rate(group["convergence_flag"]),
                "unconstrained_inside_box_rate": _rate(
                    group["unconstrained_inside_box"]
                ),
                "constrained_fallback_rate": _rate(group["constrained_fallback_used"]),
                "boundary_active_rate": _rate(group["boundary_active"]),
                "constrained_solver_failure_rate": float(
                    group["constrained_solver_failure"].mean()
                ),
                "constrained_KKT_failure_rate": float(
                    group["constrained_kkt_failure"].mean()
                ),
                "numerical_rank_loss_rate": float(group["numerical_rank_loss"].mean()),
                "median_iterations": _q(group["iterations"], 0.5),
                "p90_iterations": _q(group["iterations"], 0.9),
                "median_runtime_seconds": _q(group["runtime_seconds"], 0.5),
                "p90_runtime_seconds": _q(group["runtime_seconds"], 0.9),
                "maximum_runtime_seconds": _safe_max(group["runtime_seconds"]),
                "median_max_abs_coefficient": _q(
                    group["max_abs_coefficient"], 0.5
                ),
                "p90_max_abs_coefficient": _q(group["max_abs_coefficient"], 0.9),
            }
        )
    return pd.DataFrame(boundary_rows), pd.DataFrame(optimization_rows), roles


def _rank_summaries(
    attempts: pd.DataFrame, records: pd.DataFrame, ranks: pd.DataFrame
) -> tuple[pd.DataFrame, pd.DataFrame]:
    selected_attempts = attempts.loc[attempts["method"].eq("selected_rank")].copy()
    selected_records = records.loc[records["method"].eq("selected_rank")].copy()
    rows: list[dict[str, Any]] = []
    distributions: list[dict[str, Any]] = []
    for (dgp, n), group_attempts in selected_attempts.groupby(["dgp", "N"], sort=True):
        group_rank = ranks.loc[ranks["dgp"].eq(dgp) & ranks["N"].eq(n)].copy()
        group_records = selected_records.loc[
            selected_records["dgp"].eq(dgp) & selected_records["N"].eq(n)
        ]
        covered = group_rank.loc[_bool(group_rank["candidate_coverage"])]
        under = group_rank[
            ["A_underselected", "B_underselected", "H_underselected"]
        ].fillna(False).astype(bool).any(axis=1)
        over = group_rank[
            ["A_overselected", "B_overselected", "H_overselected"]
        ].fillna(False).astype(bool).any(axis=1)
        distribution = (
            group_rank["selected_rank_vector"].dropna().astype(str).value_counts().sort_index()
        )
        attempted = len(group_attempts)
        coverage_count = int(_bool(group_rank["candidate_coverage"]).sum())
        exact_count = int(_bool(group_rank["exact_rank_recovery"]).sum())
        cap_count = int(_bool(group_rank["rank_at_cap"]).sum())
        for rank_vector, count in distribution.items():
            distributions.append(
                {
                    "dgp": dgp,
                    "N": n,
                    "selected_rank_vector": rank_vector,
                    "count": int(count),
                    "rate": float(count / attempted),
                }
            )
        rows.append(
            {
                "dgp": dgp,
                "N": n,
                "attempted_replications": attempted,
                "valid_cap_pilots": int(_bool(group_rank["cap_pilot_converged"]).sum()),
                "pilot_multistart_disagreement_count": int(
                    _bool(group_rank["pilot_multistart_disagreement"]).sum()
                ),
                "basin_confirmation_attempted": int(
                    _bool(group_rank["cap_pilot_basin_confirmation_attempted"]).sum()
                ),
                "basin_confirmation_successful": int(
                    _bool(group_rank["cap_pilot_basin_confirmation_success"]).sum()
                ),
                "candidate_coverage_count": coverage_count,
                "candidate_coverage_rate": coverage_count / attempted,
                "P_true_rank_absent_from_candidate_set": 1.0
                - coverage_count / attempted,
                "exact_rank_recovery_count": exact_count,
                "exact_rank_recovery_rate": exact_count / attempted,
                "P_selected_not_true_given_candidate_coverage": (
                    1.0 - _rate(covered["exact_rank_recovery"])
                    if len(covered)
                    else np.nan
                ),
                "underselection_count": int(under.sum()),
                "underselection_rate": float(under.sum() / attempted),
                "overselection_count": int(over.sum()),
                "overselection_rate": float(over.sum() / attempted),
                "rank_cap_hit_count": cap_count,
                "rank_cap_hit_rate": cap_count / attempted,
                "selected_rank_distribution_json": json.dumps(distribution.to_dict()),
                "point_retained_share": _rate(group_records["retained_for_bias_rmse"]),
                "inference_retained_share": _rate(group_records["retained_for_coverage"]),
                "median_runtime_seconds": _q(
                    group_attempts["replication_runtime_seconds"], 0.5
                ),
            }
        )
    return pd.DataFrame(rows), pd.DataFrame(distributions)


def _performance_and_retention(
    records: pd.DataFrame,
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    work = records.copy()
    work["target_group"] = np.where(_bool(work["corrected"]), "broad", "local_plugin")
    work["interior_only"] = ~_bool(work["boundary_active_attempt"])
    retention_rows: list[dict[str, Any]] = []
    performance_rows: list[dict[str, Any]] = []
    interior_rows: list[dict[str, Any]] = []
    keys = ["dgp", "N", "method", "target_group"]
    for cell, group in work.groupby(keys, sort=True):
        attempted = len(group)
        common = dict(zip(keys, cell, strict=True))
        retention_rows.append(
            {
                **common,
                "target_attempts": attempted,
                "point_retained": int(_bool(group["retained_for_bias_rmse"]).sum()),
                "point_retained_share": _rate(group["retained_for_bias_rmse"]),
                "inference_retained": int(_bool(group["retained_for_coverage"]).sum()),
                "inference_retained_share": _rate(group["retained_for_coverage"]),
                "interior_only_attempts": int(group["interior_only"].sum()),
                "interior_only_point_retained_share": _rate(
                    group.loc[group["interior_only"], "retained_for_bias_rmse"]
                ),
                "interior_only_inference_retained_share": _rate(
                    group.loc[group["interior_only"], "retained_for_coverage"]
                ),
            }
        )
    for cell, group in work.groupby(["dgp", "N", "method", "target"], sort=True):
        point = group.loc[_bool(group["retained_for_bias_rmse"])].copy()
        inference = group.loc[_bool(group["retained_for_coverage"])].copy()
        errors = _num(point["estimate"]) - _num(point["truth"])
        common = dict(zip(["dgp", "N", "method", "target"], cell, strict=True))
        performance_rows.append(
            {
                **common,
                "target_group": group["target_group"].iloc[0],
                "true_value_mean": float(_num(point["truth"]).mean()) if len(point) else np.nan,
                "mean_estimate": float(_num(point["estimate"]).mean()) if len(point) else np.nan,
                "bias": float(errors.mean()) if len(point) else np.nan,
                "rmse": float(np.sqrt(np.mean(errors**2))) if len(point) else np.nan,
                "mc_sd": float(_num(point["estimate"]).std(ddof=1)) if len(point) > 1 else np.nan,
                "mean_se": float(_num(inference["standard_error"]).mean())
                if len(inference)
                else np.nan,
                "coverage": _rate(inference["covered_95pct"]),
                "mean_interval_length": float(
                    (2.0 * 1.959963984540054 * _num(inference["standard_error"])).mean()
                )
                if len(inference)
                else np.nan,
                "point_retained": len(point),
                "inference_retained": len(inference),
                "attempted": len(group),
                "extreme_finite_estimate_count": int(
                    _bool(group["extreme_estimate_flag"]).sum()
                ),
            }
        )
        for subset_name, subset in (
            ("all_valid_constrained", group),
            ("interior_only", group.loc[group["interior_only"]]),
        ):
            interior_rows.append(
                {
                    **common,
                    "subset": subset_name,
                    "attempted": len(subset),
                    "point_retained_share": _rate(subset["retained_for_bias_rmse"]),
                    "inference_retained_share": _rate(subset["retained_for_coverage"]),
                }
            )
    return (
        pd.DataFrame(retention_rows),
        pd.DataFrame(performance_rows),
        pd.DataFrame(interior_rows),
    )


def _gram_riesz(records: pd.DataFrame, inference: pd.DataFrame) -> pd.DataFrame:
    theorem = records.loc[_bool(records["headline_theorem_target"])].copy()
    join_keys = ["semantic_replication_id", "method", "target"]
    diagnostic_fields = [
        "target_supported",
        "riesz_converged",
        "riesz_target_rayleigh_quotient",
        "tangent_gram_min_eigenvalue",
        "tangent_gram_max_eigenvalue",
        "tangent_gram_condition_number",
        "variance_estimate",
        "standard_error",
        "primary_status",
    ]
    diagnostics = inference[[*join_keys, *diagnostic_fields]].rename(
        columns={field: f"diagnostic_{field}" for field in diagnostic_fields}
    )
    merged = theorem.merge(
        diagnostics,
        on=join_keys,
        how="left",
        validate="one_to_one",
    )
    merged["target_group"] = np.where(_bool(merged["corrected"]), "broad", "local_plugin")
    rows = []
    for cell, group in merged.groupby(
        ["dgp", "N", "method", "target_group"], sort=True
    ):
        supported = group.loc[_bool(group["diagnostic_target_supported"])]
        status = group["diagnostic_primary_status"].fillna("").astype(str)
        variance = _num(group["diagnostic_variance_estimate"])
        se = _num(group["diagnostic_standard_error"])
        rows.append(
            {
                **dict(
                    zip(["dgp", "N", "method", "target_group"], cell, strict=True)
                ),
                "theorem_covered_attempts": len(group),
                "target_support_rate": _rate(group["diagnostic_target_supported"]),
                "tangent_Gram_failure_rate": float(
                    status.isin(
                        ["tangent_gram_eigensolver_failure", "tangent_gram_nearly_singular"]
                    ).mean()
                ),
                "Riesz_failure_rate_conditional_on_support": (
                    1.0 - _rate(supported["diagnostic_riesz_converged"])
                    if len(supported)
                    else np.nan
                ),
                "invalid_variance_rate": float((~np.isfinite(variance) | (variance <= 0)).mean()),
                "nonfinite_SE_rate": float((~np.isfinite(se)).mean()),
                "minimum_Gram_eigenvalue": (
                    float(_finite(group["diagnostic_tangent_gram_min_eigenvalue"]).min())
                    if len(_finite(group["diagnostic_tangent_gram_min_eigenvalue"]))
                    else np.nan
                ),
                "median_Gram_eigenvalue": _q(
                    group["diagnostic_tangent_gram_min_eigenvalue"], 0.5
                ),
                "p90_condition_number": _q(
                    group["diagnostic_tangent_gram_condition_number"], 0.9
                ),
                "maximum_condition_number": _safe_max(
                    group["diagnostic_tangent_gram_condition_number"]
                ),
                "median_riesz_target_rayleigh_quotient": _q(
                    group["diagnostic_riesz_target_rayleigh_quotient"], 0.5
                ),
            }
        )
    return pd.DataFrame(rows)


def _split_behavior(fits: pd.DataFrame, records: pd.DataFrame) -> pd.DataFrame:
    split = fits.loc[
        fits["fit_type"].str.startswith("time_split_", na=False)
        | fits["fit_type"].str.startswith("unit_split_", na=False)
    ].copy()
    split["split_dimension"] = np.where(
        split["fit_type"].str.startswith("time_"), "time", "unit"
    )
    rows = []
    for cell, group in split.groupby(["dgp", "N", "method", "split_dimension"], sort=True):
        rows.append(
            {
                **dict(
                    zip(["dgp", "N", "method", "split_dimension"], cell, strict=True)
                ),
                "attempted_split_fits": len(group),
                "expected_split_fits": 6,
                "completion_rate": len(group) / 6.0,
                "validity_rate": _rate(group["convergence_flag"]),
                "constrained_fallback_rate": _rate(group["constrained_fallback_used"]),
                "boundary_active_rate": _rate(group["boundary_active"]),
                "numerical_rank_loss_rate": float(
                    (
                        group["requested_rank"].notna()
                        & group["numerical_rank"].notna()
                        & group["requested_rank"].ne(group["numerical_rank"])
                    ).mean()
                ),
                "median_stationarity_or_KKT_residual": _q(
                    group["constrained_KKT_residual"].fillna(group["stationarity_residual"]),
                    0.5,
                ),
                "median_runtime_seconds": _q(group["runtime_seconds"], 0.5),
            }
        )
    broad = records.loc[_bool(records["corrected"])].copy()
    per_rep = broad.groupby(
        ["semantic_replication_id", "method"], sort=True
    ).agg(
        split_fit_count_min=("split_coefficient_fit_count", "min"),
        split_fit_count_max=("split_coefficient_fit_count", "max"),
        time_assignments=("time_split_assignments_json", "nunique"),
        unit_assignments=("unit_split_assignments_json", "nunique"),
    )
    if not (
        per_rep["split_fit_count_min"].eq(4)
        & per_rep["split_fit_count_max"].eq(4)
        & per_rep["time_assignments"].eq(1)
        & per_rep["unit_assignments"].eq(1)
    ).all():
        raise ValueError("four shared split fits did not reconcile across broad targets")
    return pd.DataFrame(rows)


def _fixed_reliability(
    attempts: pd.DataFrame,
    records: pd.DataFrame,
    fits: pd.DataFrame,
    gram: pd.DataFrame,
) -> pd.DataFrame:
    rows = []
    fixed_attempts = attempts.loc[attempts["method"].eq("fixed_rank")]
    fixed_records = records.loc[records["method"].eq("fixed_rank")].copy()
    fixed_records["target_group"] = np.where(
        _bool(fixed_records["corrected"]), "broad", "local_plugin"
    )
    fixed_fits = fits.loc[fits["method"].eq("fixed_rank")]
    for (dgp, n), group_attempts in fixed_attempts.groupby(["dgp", "N"], sort=True):
        cell_records = fixed_records.loc[
            fixed_records["dgp"].eq(dgp) & fixed_records["N"].eq(n)
        ]
        cell_fits = fixed_fits.loc[fixed_fits["dgp"].eq(dgp) & fixed_fits["N"].eq(n)]
        for target_group, group_records in cell_records.groupby("target_group", sort=True):
            gram_row = gram.loc[
                gram["dgp"].eq(dgp)
                & gram["N"].eq(n)
                & gram["method"].eq("fixed_rank")
                & gram["target_group"].eq(target_group)
            ]
            time_fits = cell_fits.loc[
                cell_fits["fit_type"].str.startswith("time_split_", na=False)
            ]
            unit_fits = cell_fits.loc[
                cell_fits["fit_type"].str.startswith("unit_split_", na=False)
            ]
            rows.append(
                {
                    "dgp": dgp,
                    "N": n,
                    "target_group": target_group,
                    "attempted_replications": len(group_attempts),
                    "valid_full_panel_fits": int(
                        (~group_attempts["primary_status"].isin(
                            [
                                "full_fit_failure",
                                "constrained_solver_failure",
                                "constrained_feasibility_failure",
                                "constrained_optimality_failure",
                                "nonfinite_constrained_solution",
                            ]
                        )).sum()
                    ),
                    "constrained_fallback_rate": _rate(
                        group_attempts["constrained_fallback_used"]
                    ),
                    "boundary_active_rate": _rate(group_attempts["boundary_active"]),
                    "constrained_solver_failure_rate": float(
                        group_attempts["primary_status"].eq(
                            "constrained_solver_failure"
                        ).mean()
                    ),
                    "four_required_split_completion_rate": float(
                        group_records.groupby("semantic_replication_id")[
                            "split_coefficient_fit_count"
                        ].min().eq(4).mean()
                    )
                    if target_group == "broad"
                    else np.nan,
                    "time_split_validity_rate": _rate(time_fits["convergence_flag"]),
                    "unit_split_validity_rate": _rate(unit_fits["convergence_flag"]),
                    "target_support_rate": (
                        float(gram_row["target_support_rate"].iloc[0])
                        if len(gram_row)
                        else np.nan
                    ),
                    "tangent_Gram_failure_rate": (
                        float(gram_row["tangent_Gram_failure_rate"].iloc[0])
                        if len(gram_row)
                        else np.nan
                    ),
                    "Riesz_failure_rate": (
                        float(gram_row["Riesz_failure_rate_conditional_on_support"].iloc[0])
                        if len(gram_row)
                        else np.nan
                    ),
                    "invalid_variance_rate": (
                        float(gram_row["invalid_variance_rate"].iloc[0])
                        if len(gram_row)
                        else np.nan
                    ),
                    "point_retained_share": _rate(group_records["retained_for_bias_rmse"]),
                    "inference_retained_share": _rate(group_records["retained_for_coverage"]),
                    "median_runtime_seconds": _q(
                        group_attempts["replication_runtime_seconds"], 0.5
                    ),
                }
            )
    return pd.DataFrame(rows)


def _accounting(records: pd.DataFrame) -> pd.DataFrame:
    rows = []
    keys = ["dgp", "N", "method", "target"]
    for cell, group in records.groupby(keys, sort=True):
        attempted = len(group)
        point = int(_bool(group["retained_for_bias_rmse"]).sum())
        inference = int(_bool(group["retained_for_coverage"]).sum())
        success = int(group["primary_status"].eq("success").sum())
        iff = _bool(group["inference_valid"]).eq(group["primary_status"].eq("success")).all()
        rows.append(
            {
                **dict(zip(keys, cell, strict=True)),
                "R_attempted": attempted,
                "R_point": point,
                "R_inference": inference,
                "success_count": success,
                "retention_order_pass": 0 <= inference <= point <= attempted,
                "inference_valid_iff_success_pass": bool(iff),
            }
        )
    result = pd.DataFrame(rows)
    if not (
        result["R_attempted"].eq(3)
        & result["retention_order_pass"]
        & result["inference_valid_iff_success_pass"]
        & result["R_inference"].eq(result["success_count"])
    ).all():
        raise ValueError("target accounting reconciliation failed")
    return result


def _markdown_table(frame: pd.DataFrame, columns: list[str]) -> list[str]:
    selected = frame[columns].copy()
    for column in selected.select_dtypes(include=["float"]).columns:
        selected[column] = selected[column].map(
            lambda value: "" if pd.isna(value) else f"{value:.4f}"
        )
    header = "| " + " | ".join(columns) + " |"
    divider = "| " + " | ".join("---" for _ in columns) + " |"
    rows = [
        "| " + " | ".join(str(value) for value in row) + " |"
        for row in selected.itertuples(index=False, name=None)
    ]
    return [header, divider, *rows]


def _report(
    output: Path,
    semantic: pd.DataFrame,
    boundary: pd.DataFrame,
    fixed: pd.DataFrame,
    ranks: pd.DataFrame,
    retention: pd.DataFrame,
    gram: pd.DataFrame,
    optimization: pd.DataFrame,
    accounting: pd.DataFrame,
) -> str:
    fixed_success = float(
        fixed.groupby(["dgp", "N"])["valid_full_panel_fits"].first().sum() / 24.0
    )
    constrained_failure = float(
        boundary[
            [
                "constrained_solver_failure_count",
                "constrained_feasibility_failure_count",
                "constrained_optimality_failure_count",
            ]
        ].sum().sum()
        / max(1, boundary["attempted_fits"].sum())
    )
    fixed_broad = fixed.loc[fixed["target_group"].eq("broad")]
    split_operational = bool(
        fixed_broad["four_required_split_completion_rate"].fillna(0).ge(0.9).all()
    )
    gram_ok = bool(
        gram["tangent_Gram_failure_rate"].fillna(1).le(0.1).all()
        and gram["Riesz_failure_rate_conditional_on_support"].fillna(1).le(0.1).all()
    )
    coverage = float(
        np.average(
            ranks["candidate_coverage_rate"], weights=ranks["attempted_replications"]
        )
    )
    cap_hits = float(
        np.average(ranks["rank_cap_hit_rate"], weights=ranks["attempted_replications"])
    )
    by_size = ranks.groupby("N", sort=True).apply(
        lambda group: pd.Series(
            {
                "underselection_rate": np.average(
                    group["underselection_rate"],
                    weights=group["attempted_replications"],
                ),
                "overselection_rate": np.average(
                    group["overselection_rate"],
                    weights=group["attempted_replications"],
                ),
                "exact_rank_recovery_rate": np.average(
                    group["exact_rank_recovery_rate"],
                    weights=group["attempted_replications"],
                ),
            }
        ),
        include_groups=False,
    )
    nondegenerate = bool(
        by_size["underselection_rate"].lt(0.8).all()
        and by_size["overselection_rate"].lt(0.8).all()
        and by_size["exact_rank_recovery_rate"].gt(0.1).all()
    )
    accounting_ok = bool(
        accounting["retention_order_pass"].all()
        and accounting["inference_valid_iff_success_pass"].all()
    )
    criteria = {
        "fixed_rank_reliable": fixed_success >= 0.9,
        "constrained_failures_rare": constrained_failure <= 0.05,
        "split_broad_operational": split_operational,
        "gram_riesz_failures_rare": gram_ok,
        "candidate_coverage_high": coverage >= 0.8,
        "selected_behavior_not_degenerate": nondegenerate,
        "rank_cap_hits_not_systematic": cap_hits < 0.5,
        "accounting_reconciles": accounting_ok,
    }
    decision = "GO" if all(criteria.values()) else "NO-GO"
    plausible = coverage >= 0.8 and nondegenerate and cap_hits < 0.5

    lines = [
        "# Final theorem-aligned statistical preflight",
        "",
        "This is a 3-replication-per-cell diagnostic preflight, not publication Monte Carlo evidence.",
        "No trimming, winsorization, magnitude filtering, or runtime filtering was applied.",
        "The Parquet replication records are lossless. Their CSV companions omit only the two "
        "very large nested rank-selection diagnostic columns, which remain in Parquet and the "
        "dedicated rank/fit files.",
        "",
        "## Exact commands",
        "",
        "```text",
        *DRY_COMMANDS,
        *RUN_COMMANDS,
        "```",
        "",
        "## Design and matched draws",
        "",
        f"- Unique semantic DGP realizations: {len(semantic)} (required 24).",
        "- Method evaluations: 48 (24 fixed rank, 24 selected rank).",
        "- All fixed/selected DGP realization hashes match: PASS.",
        f"- Frozen calibration SHA-256: `{FROZEN_HASH}`.",
        "- `B=10`, `c_B=1`; every requested frozen cell satisfies `C_Theta <= 9`.",
        "",
        "## Fixed-rank reliability",
        "",
        *_markdown_table(
            fixed,
            [
                "dgp",
                "N",
                "target_group",
                "valid_full_panel_fits",
                "constrained_fallback_rate",
                "boundary_active_rate",
                "four_required_split_completion_rate",
                "point_retained_share",
                "inference_retained_share",
                "median_runtime_seconds",
            ],
        ),
        "",
        "## Selected-rank reliability and complete distributions",
        "",
        *_markdown_table(
            ranks,
            [
                "dgp",
                "N",
                "candidate_coverage_rate",
                "exact_rank_recovery_rate",
                "underselection_rate",
                "overselection_rate",
                "rank_cap_hit_rate",
                "selected_rank_distribution_json",
                "point_retained_share",
                "inference_retained_share",
                "median_runtime_seconds",
            ],
        ),
        "",
        "## Boundary and constrained-estimator behavior",
        "",
        *_markdown_table(
            boundary,
            [
                "dgp",
                "N",
                "method",
                "fit_role",
                "attempted_fits",
                "unconstrained_outside_box_rate",
                "constrained_fallback_rate",
                "boundary_active_rate",
                "constrained_solver_failure_rate",
                "constrained_feasibility_failure_rate",
                "constrained_optimality_failure_rate",
            ],
        ),
        "",
        "Successful boundary-active constrained estimates are retained in the primary results. "
        "The interior-only comparison is secondary and is stored separately.",
        "",
        "## Target retention",
        "",
        *_markdown_table(
            retention,
            [
                "dgp",
                "N",
                "method",
                "target_group",
                "point_retained_share",
                "inference_retained_share",
                "interior_only_attempts",
                "interior_only_inference_retained_share",
            ],
        ),
        "",
        "## Gram/Riesz diagnostics for theorem-covered targets",
        "",
        *_markdown_table(
            gram,
            [
                "dgp",
                "N",
                "method",
                "target_group",
                "target_support_rate",
                "tangent_Gram_failure_rate",
                "Riesz_failure_rate_conditional_on_support",
                "invalid_variance_rate",
                "minimum_Gram_eigenvalue",
                "p90_condition_number",
            ],
        ),
        "",
        "## Runtime and optimization",
        "",
        "Lossless fit-role diagnostics are in `optimization_summary.csv`; slow and extreme finite "
        "fits remain included. Selected final full-panel rows are a derived role view of the "
        "chosen candidate post-refit and are therefore intentionally non-additive with candidate rows.",
        "",
        "## Accounting",
        "",
        f"- Target-cell accounting rows: {len(accounting)}; all have `R_attempted=3`.",
        "- `R_inference <= R_point <= R_attempted`: PASS.",
        "- `inference_valid iff primary_status=success`: PASS.",
        "- Every method-replication contains all 18 requested target records: PASS.",
        "- Exactly four shared split fits per broad-target replication: PASS.",
        "",
        "## Implementation finding",
        "",
        "The numerical estimator did not require correction. The audit identified a fit-diagnostic "
        "labeling defect: additional deterministic fixed-rank full-panel starts were emitted as "
        "generic `coefficient_fit`, and selected-fit diagnostic context could overwrite some split "
        "labels. The combined output relabels those records losslessly by execution position. Future "
        "instrumentation labels all pre-split fixed starts `full_fixed_rank` and preserves split "
        "labels before applying selected-fit context. Explicit aliases `max_abs_coefficient` and "
        "`constrained_runtime_seconds` were also added.",
        "",
        "## Medium recommendation",
        "",
        f"- Provisional `c_kappa=4e-6` remains plausible: {'YES' if plausible else 'NO'}.",
        f"- Recommendation for the 100-replication medium diagnostic: **{decision}**.",
        "- Criterion results: `" + json.dumps(criteria, sort_keys=True) + "`.",
        "- By-size selected-rank behavior: `"
        + json.dumps(
            {
                str(index): {
                    key: float(value)
                    for key, value in row.items()
                }
                for index, row in by_size.to_dict(orient="index").items()
            },
            sort_keys=True,
        )
        + "`.",
        "",
        "The medium diagnostic was not launched.",
    ]
    (output / "preflight_report.md").write_text("\n".join(lines) + "\n", encoding="utf-8")
    return decision


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--fixed-root", required=True, type=Path)
    parser.add_argument("--selected-root", required=True, type=Path)
    parser.add_argument("--output-root", required=True, type=Path)
    args = parser.parse_args()
    output = args.output_root
    output.mkdir(parents=True, exist_ok=True)

    frames = {
        "fixed_rank": _read_method(args.fixed_root),
        "selected_rank": _read_method(args.selected_root),
    }
    attempts = pd.concat([item["attempts"] for item in frames.values()], ignore_index=True)
    records = pd.concat([item["records"] for item in frames.values()], ignore_index=True)
    fits = _normalize_fit_rows(
        pd.concat([item["fits"] for item in frames.values()], ignore_index=True, sort=False)
    )
    inference = pd.concat(
        [item["inference"] for item in frames.values()], ignore_index=True, sort=False
    )
    ranks = frames["selected_rank"]["rank"].copy()

    if len(attempts) != 48 or attempts.groupby("method").size().to_dict() != {
        "fixed_rank": 24,
        "selected_rank": 24,
    }:
        raise ValueError("authorized 48 method-replication evaluation count failed")
    expected_ids = {
        f"dgp{dgp}_N{n}_T{n}_r{replication:05d}_truth1-1-1"
        for dgp in EXPECTED_DGPS
        for n in EXPECTED_SIZES
        for replication in EXPECTED_REPLICATIONS
    }
    if set(attempts["semantic_replication_id"]) != expected_ids:
        raise ValueError("semantic replication ID set does not match the authorized design")
    matching = attempts.pivot(
        index="semantic_replication_id", columns="method", values="dgp_realization_hash"
    ).reset_index()
    matching["hashes_match"] = matching[["fixed_rank", "selected_rank"]].nunique(
        axis=1, dropna=False
    ).eq(1)
    if len(matching) != 24 or not matching["hashes_match"].all():
        raise ValueError("matched DGP hash verification failed")
    semantic = attempts[
        ["semantic_replication_id", "dgp", "N", "T", "replication"]
    ].drop_duplicates().sort_values(["dgp", "N", "replication"])
    semantic = semantic.merge(
        matching[["semantic_replication_id", "fixed_rank", "selected_rank", "hashes_match"]],
        on="semantic_replication_id",
        validate="one_to_one",
    )

    if len(records) != 48 * EXPECTED_TARGETS:
        raise ValueError("requested target records silently disappeared")
    target_counts = records.groupby(["semantic_replication_id", "method"])["target"].nunique()
    if not target_counts.eq(EXPECTED_TARGETS).all():
        raise ValueError("a method-replication is missing requested target records")

    manifest_rows = []
    for method, root in (
        ("fixed_rank", args.fixed_root),
        ("selected_rank", args.selected_root),
    ):
        manifest = json.loads((root / "run_manifest.json").read_text(encoding="utf-8"))
        if (
            manifest["B"] != 10.0
            or manifest["c_B"] != 1.0
            or manifest["frozen_calibration_hash"] != FROZEN_HASH
            or len(manifest["calibration_cells"]) != 8
            or any(float(cell["C_Theta"]) > 9.0 for cell in manifest["calibration_cells"])
        ):
            raise ValueError(f"finalized manifest validation failed for {method}")
        for cell in manifest["calibration_cells"]:
            manifest_rows.append({"method": method, **cell["calibration_cell"], **cell})
    manifest_summary = pd.json_normalize(manifest_rows)

    accounting = _accounting(records)
    boundary, optimization, fit_roles = _boundary_and_optimization(fits, attempts)
    rank_summary, rank_distribution = _rank_summaries(attempts, records, ranks)
    retention, performance, interior = _performance_and_retention(records)
    gram = _gram_riesz(records, inference)
    split = _split_behavior(fits, records)
    fixed = _fixed_reliability(attempts, records, fits, gram)

    _write_replication_records(records, output)
    for method, root in (
        ("fixed_rank", args.fixed_root),
        ("selected_rank", args.selected_root),
    ):
        _write_replication_records(frames[method]["records"], root)
    _write_both(fits, output, "fit_diagnostics")
    _write_both(inference, output, "inference_diagnostics")
    boundary.to_csv(output / "boundary_activity_summary.csv", index=False, lineterminator="\n")
    rank_summary.to_csv(output / "rank_selection_summary.csv", index=False, lineterminator="\n")
    retention.to_csv(output / "target_retention_summary.csv", index=False, lineterminator="\n")
    optimization.to_csv(output / "optimization_summary.csv", index=False, lineterminator="\n")
    gram.to_csv(output / "gram_riesz_summary.csv", index=False, lineterminator="\n")
    semantic.to_csv(output / "semantic_dgp_ids_and_hashes.csv", index=False, lineterminator="\n")
    rank_distribution.to_csv(output / "selected_rank_distribution.csv", index=False, lineterminator="\n")
    performance.to_csv(output / "preliminary_performance_summary.csv", index=False, lineterminator="\n")
    interior.to_csv(output / "interior_only_comparison.csv", index=False, lineterminator="\n")
    split.to_csv(output / "split_fit_summary.csv", index=False, lineterminator="\n")
    fixed.to_csv(output / "fixed_rank_reliability.csv", index=False, lineterminator="\n")
    accounting.to_csv(output / "accounting_reconciliation.csv", index=False, lineterminator="\n")
    manifest_summary.to_csv(output / "manifest_calibration_cells.csv", index=False, lineterminator="\n")
    fit_roles.to_parquet(output / "fit_role_diagnostics.parquet", index=False)
    decision = _report(
        output,
        semantic,
        boundary,
        fixed,
        rank_summary,
        retention,
        gram,
        optimization,
        accounting,
    )
    digest = hashlib.sha256((output / "replication_records.parquet").read_bytes()).hexdigest()
    (output / "audit_manifest.json").write_text(
        json.dumps(
            {
                "authorized_unique_dgp_realizations": 24,
                "authorized_method_replication_evaluations": 48,
                "semantic_id_count": len(semantic),
                "matched_hashes_pass": bool(semantic["hashes_match"].all()),
                "target_record_count": len(records),
                "accounting_pass": True,
                "frozen_calibration_hash": FROZEN_HASH,
                "combined_replication_records_sha256": digest,
                "medium_recommendation": decision,
                "dry_run_commands": DRY_COMMANDS,
                "executed_commands": RUN_COMMANDS,
            },
            indent=2,
        )
        + "\n",
        encoding="utf-8",
    )
    print(f"Finalized preflight audit complete: medium recommendation {decision}")


if __name__ == "__main__":
    main()
