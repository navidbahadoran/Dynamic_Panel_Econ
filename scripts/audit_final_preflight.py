"""Audit the final matched fixed/selected Monte Carlo preflight."""

from __future__ import annotations

import argparse
import json
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd


def _read_chunks(root: Path, subdirectory: str) -> pd.DataFrame:
    paths = sorted((root / subdirectory).glob("*.parquet"))
    return pd.concat([pd.read_parquet(path) for path in paths], ignore_index=True)


def _as_bool(values: pd.Series) -> pd.Series:
    return values.fillna(False).astype(bool)


def _rate(values: pd.Series) -> float:
    return float(_as_bool(values).mean()) if len(values) else float("nan")


def _finite(values: pd.Series) -> pd.Series:
    return pd.to_numeric(values, errors="coerce").replace([np.inf, -np.inf], np.nan).dropna()


def _quantile(values: pd.Series, q: float) -> float:
    finite = _finite(values)
    return float(finite.quantile(q)) if len(finite) else float("nan")


def _method_frames(label: str, root: Path) -> dict[str, pd.DataFrame]:
    result = {
        "attempts": pd.read_parquet(root / "attempted_replications.parquet"),
        "records": pd.read_parquet(root / "replication_records.parquet"),
        "fits": pd.read_parquet(root / "fit_diagnostics.parquet"),
        "inference": pd.read_parquet(root / "inference_diagnostics.parquet"),
        "rank": _read_chunks(root, "rank"),
        "raw": _read_chunks(root, "raw"),
    }
    for frame in result.values():
        frame["preflight_method"] = label
    return result


def _optimization_summary(label: str, fits: pd.DataFrame) -> tuple[dict[str, Any], pd.DataFrame]:
    coefficient = fits.loc[fits["fit_type"].ne("nuclear_path")].copy()
    stable = coefficient["objective_stability_pass"].dropna()
    accepted_stationarity = _finite(
        coefficient.loc[_as_bool(coefficient["stationarity_pass"]), "stationarity_residual"]
    )
    accepted_gap = _finite(
        coefficient.loc[_as_bool(coefficient["objective_stability_pass"]), "objective_stability_gap"]
    )
    numerical_loss = (
        coefficient["requested_rank"].notna()
        & coefficient["numerical_rank"].notna()
        & coefficient["requested_rank"].ne(coefficient["numerical_rank"])
    )
    row = {
        "method": label,
        "coefficient_fit_count": len(coefficient),
        "convergence_rate": _rate(coefficient["convergence_flag"]),
        "iteration_cap_hit_rate": _rate(coefficient["iteration_cap_hit"]),
        "stationarity_failure_rate": 1.0 - _rate(coefficient["stationarity_pass"]),
        "coefficient_bound_hit_rate": _rate(coefficient["coefficient_bound_hit"]),
        "objective_stability_evaluated_count": len(stable),
        "objective_stability_failure_rate": 1.0 - _rate(stable),
        "numerical_rank_loss_rate": float(numerical_loss.mean()) if len(coefficient) else np.nan,
        "median_runtime_seconds": _quantile(coefficient["runtime_seconds"], 0.5),
        "p90_runtime_seconds": _quantile(coefficient["runtime_seconds"], 0.9),
        "maximum_runtime_seconds": float(_finite(coefficient["runtime_seconds"]).max()),
        "median_stationarity_residual": _quantile(coefficient["stationarity_residual"], 0.5),
        "maximum_accepted_stationarity_residual": (
            float(accepted_stationarity.max()) if len(accepted_stationarity) else np.nan
        ),
        "median_objective_stability_gap": _quantile(coefficient["objective_stability_gap"], 0.5),
        "maximum_accepted_objective_stability_gap": (
            float(accepted_gap.max()) if len(accepted_gap) else np.nan
        ),
    }
    valid = coefficient.loc[
        _as_bool(coefficient["convergence_flag"])
        & _as_bool(coefficient["stationarity_pass"])
        & ~_as_bool(coefficient["coefficient_bound_hit"])
        & coefficient["objective_stability_pass"].fillna(True).astype(bool)
    ].copy()
    valid["method"] = label
    return row, valid


def _inference_summary(label: str, records: pd.DataFrame, inference: pd.DataFrame) -> dict[str, Any]:
    support = _as_bool(inference["target_supported"])
    gram_failed = (
        inference["tangent_gram_min_eigenvalue"].isna()
        | (pd.to_numeric(inference["tangent_gram_min_eigenvalue"], errors="coerce") <= 1e-10)
    )
    supported = inference.loc[support]
    variance = pd.to_numeric(inference["variance_estimate"], errors="coerce")
    standard_error = pd.to_numeric(inference["standard_error"], errors="coerce")
    return {
        "method": label,
        "target_attempts": len(records),
        "target_support_rate": _rate(inference["target_supported"]),
        "tangent_gram_failure_rate": float(gram_failed.mean()) if len(inference) else np.nan,
        "riesz_failure_rate_conditional_on_support": (
            1.0 - _rate(supported["riesz_converged"]) if len(supported) else np.nan
        ),
        "invalid_variance_rate": float((~np.isfinite(variance) | (variance <= 0)).mean())
        if len(inference)
        else np.nan,
        "nonfinite_standard_error_rate": float((~np.isfinite(standard_error)).mean())
        if len(inference)
        else np.nan,
        "point_retained_count": int(_as_bool(records["retained_for_bias_rmse"]).sum()),
        "point_retained_share": _rate(records["retained_for_bias_rmse"]),
        "inference_retained_count": int(_as_bool(records["retained_for_coverage"]).sum()),
        "inference_retained_share": _rate(records["retained_for_coverage"]),
        "extreme_estimate_flag_count": int(_as_bool(records["extreme_estimate_flag"]).sum()),
    }


def _selected_summary(label: str, attempts: pd.DataFrame, rank: pd.DataFrame) -> dict[str, Any]:
    selected = rank["selected_rank_vector"].dropna().astype(str)
    true_present = rank.loc[_as_bool(rank["candidate_coverage"])]
    exact = _as_bool(rank["exact_rank_recovery"])
    under = _as_bool(rank["A_underselected"]) | _as_bool(rank["B_underselected"]) | _as_bool(
        rank["H_underselected"]
    )
    over = _as_bool(rank["A_overselected"]) | _as_bool(rank["B_overselected"]) | _as_bool(
        rank["H_overselected"]
    )
    return {
        "method": label,
        "attempted": len(attempts),
        "pilot_success": len(rank),
        "candidate_coverage_count": int(_as_bool(rank["candidate_coverage"]).sum()),
        "candidate_coverage_share_of_pilot_success": _rate(rank["candidate_coverage"]),
        "true_rank_absent_count": int((~_as_bool(rank["candidate_coverage"])).sum()),
        "ic_choice_failure_given_coverage_count": int(
            (~_as_bool(true_present["exact_rank_recovery"])).sum()
        ),
        "ic_choice_failure_given_coverage_share": (
            1.0 - _rate(true_present["exact_rank_recovery"]) if len(true_present) else np.nan
        ),
        "exact_recovery_count": int(exact.sum()),
        "underselection_count": int(under.sum()),
        "overselection_count": int(over.sum()),
        "selected_zero_vector_count": int(selected.eq("[0, 0, 0]").sum()),
        "rank_cap_hit_count": int(_as_bool(rank["rank_at_cap"]).sum()),
        "selected_rank_distribution": json.dumps(selected.value_counts().sort_index().to_dict()),
        "median_replication_runtime_seconds": _quantile(
            attempts["replication_runtime_seconds"], 0.5
        ),
        "total_replication_runtime_seconds": float(
            _finite(attempts["replication_runtime_seconds"]).sum()
        ),
    }


def _replication_diagnostics(
    label: str, item: dict[str, pd.DataFrame]
) -> pd.DataFrame:
    keys = ["dgp", "N", "T", "replication"]
    attempts = item["attempts"].copy()
    rank = item["rank"].copy()
    records = item["records"].copy()
    inference = item["inference"].copy()
    full_fit_types = ["full_fixed_rank"] if label == "fixed_rank" else ["coefficient_fit"]
    full_fits = item["fits"].loc[item["fits"]["fit_type"].isin(full_fit_types)].copy()
    if label != "fixed_rank" and len(full_fits):
        selected_rank = rank[keys + ["selected_rank_vector"]]
        full_fits = full_fits.merge(selected_rank, on=keys, how="inner")
        full_fits = full_fits.loc[
            full_fits["requested_rank"].eq(full_fits["selected_rank_vector"])
        ]
        full_fits = full_fits.sort_values("objective_final").drop_duplicates(keys)
    elif len(full_fits):
        full_fits = full_fits.sort_values("objective_final").drop_duplicates(keys)
    fit_columns = keys + [
        "requested_rank", "numerical_rank", "convergence_flag", "iteration_cap_hit",
        "stationarity_residual", "stationarity_pass", "coefficient_envelope_ratio",
        "coefficient_bound_hit", "objective_stability_gap", "objective_stability_pass",
    ]
    full_fits = full_fits[[column for column in fit_columns if column in full_fits]]

    retained = records.groupby(keys, sort=True).agg(
        valid_point_estimate_count=("retained_for_bias_rmse", lambda x: int(_as_bool(x).sum())),
        valid_inference_count=("retained_for_coverage", lambda x: int(_as_bool(x).sum())),
        target_count=("target", "size"),
        extreme_estimate_flag_count=("extreme_estimate_flag", lambda x: int(_as_bool(x).sum())),
    ).reset_index()
    if len(inference):
        inference_agg = inference.groupby(keys, sort=True).agg(
            target_support_all=("target_supported", lambda x: bool(_as_bool(x).all())),
            tangent_gram_min_eigenvalue=("tangent_gram_min_eigenvalue", "min"),
            tangent_gram_condition_number_max=("tangent_gram_condition_number", "max"),
            riesz_converged_all=("riesz_converged", lambda x: bool(_as_bool(x).all())),
        ).reset_index()
    else:
        inference_agg = pd.DataFrame(columns=keys)

    diagnostic = attempts.merge(full_fits, on=keys, how="left")
    diagnostic = diagnostic.merge(retained, on=keys, how="left")
    diagnostic = diagnostic.merge(inference_agg, on=keys, how="left")
    rank_columns = keys + [
        "cap_pilot_attempted_route_count", "cap_pilot_stable_route_count",
        "cap_pilot_rank", "candidate_coverage", "selected_rank_vector",
        "exact_rank_recovery", "rank_at_cap", "A_underselected", "B_underselected",
        "H_underselected", "A_overselected", "B_overselected", "H_overselected",
    ]
    available_rank_columns = [column for column in rank_columns if column in rank]
    diagnostic = diagnostic.merge(rank[available_rank_columns], on=keys, how="left", suffixes=("", "_rank"))
    failures = item["raw"].loc[item["raw"]["record_type"].eq("failure")]
    failure_columns = keys + [
        "cap_pilot_attempted_route_count", "cap_pilot_stable_route_count"
    ]
    available_failure_columns = [column for column in failure_columns if column in failures]
    if len(failures) and len(available_failure_columns) > len(keys):
        failure_routes = failures[available_failure_columns].drop_duplicates(keys)
        diagnostic = diagnostic.merge(failure_routes, on=keys, how="left", suffixes=("", "_failure"))
        for column in ("cap_pilot_attempted_route_count", "cap_pilot_stable_route_count"):
            failure_column = f"{column}_failure"
            if failure_column in diagnostic:
                diagnostic[column] = diagnostic[column].fillna(diagnostic[failure_column])
                diagnostic = diagnostic.drop(columns=failure_column)
    diagnostic["method"] = label
    diagnostic["dgp_generation_calibration_success"] = diagnostic["completed_dgp_replication"]
    if label != "fixed_rank":
        diagnostic["cap_pilot_success"] = diagnostic["cap_pilot_rank"].notna()
        diagnostic["underselection"] = diagnostic[
            ["A_underselected", "B_underselected", "H_underselected"]
        ].fillna(False).any(axis=1)
        diagnostic["overselection"] = diagnostic[
            ["A_overselected", "B_overselected", "H_overselected"]
        ].fillna(False).any(axis=1)
    return diagnostic


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--fixed", required=True, type=Path)
    parser.add_argument("--selected-3e-6", required=True, type=Path)
    parser.add_argument("--selected-4e-6", required=True, type=Path)
    parser.add_argument("--output-dir", required=True, type=Path)
    args = parser.parse_args()
    roots = {
        "fixed_rank": args.fixed,
        "selected_3e-6": args.selected_3e_6,
        "selected_4e-6": args.selected_4e_6,
    }
    frames = {label: _method_frames(label, root) for label, root in roots.items()}
    output = args.output_dir
    output.mkdir(parents=True, exist_ok=True)

    attempts = pd.concat([item["attempts"] for item in frames.values()], ignore_index=True)
    matching = attempts.pivot(
        index="semantic_replication_id", columns="preflight_method", values="dgp_realization_hash"
    ).reset_index()
    matching["all_hashes_equal"] = matching[list(roots)].nunique(axis=1, dropna=False).eq(1)
    matching.to_csv(output / "matched_dgp_hashes.csv", index=False)
    if len(matching) != 12 or not matching["all_hashes_equal"].all():
        raise ValueError("matched DGP realization hash audit failed")

    accounting_rows = []
    for (method, dgp), group in attempts.groupby(["preflight_method", "dgp"], sort=True):
        counts = group["primary_status"].value_counts().to_dict()
        row = {"method": method, "dgp": dgp, "attempted": len(group)}
        row.update({f"status_{key}": value for key, value in counts.items()})
        row["reconciled"] = len(group) == sum(counts.values())
        accounting_rows.append(row)
    accounting = pd.DataFrame(accounting_rows).fillna(0)
    accounting.to_csv(output / "failure_accounting.csv", index=False)
    if not accounting["reconciled"].all() or len(attempts) != 36:
        raise ValueError("attempted replication accounting failed")

    optimization_rows = []
    valid_fits = []
    inference_rows = []
    selected_rows = []
    for label, item in frames.items():
        optimization, valid = _optimization_summary(label, item["fits"])
        optimization_rows.append(optimization)
        valid_fits.append(valid)
        inference_rows.append(_inference_summary(label, item["records"], item["inference"]))
        if label != "fixed_rank":
            selected_rows.append(_selected_summary(label, item["attempts"], item["rank"]))
    pd.DataFrame(optimization_rows).to_csv(output / "optimization_summary.csv", index=False)
    pd.DataFrame(inference_rows).to_csv(output / "inference_summary.csv", index=False)
    selected_summary = pd.DataFrame(selected_rows)
    selected_summary.to_csv(output / "selected_rank_summary.csv", index=False)

    valid_fit_frame = pd.concat(valid_fits, ignore_index=True)
    identification = [
        "method", "semantic_replication_id", "dgp", "replication", "fit_type",
        "requested_rank", "numerical_rank", "runtime_seconds", "stationarity_residual",
        "objective_stability_gap",
    ]
    timed_valid_fits = valid_fit_frame.loc[
        pd.to_numeric(valid_fit_frame["runtime_seconds"], errors="coerce").notna()
    ]
    timed_valid_fits.sort_values("runtime_seconds", ascending=False).head(5)[
        identification
    ].to_csv(output / "five_slowest_valid_fits.csv", index=False)
    valid_fit_frame.sort_values("stationarity_residual", ascending=False).head(5)[
        identification
    ].to_csv(output / "five_largest_stationarity_residuals.csv", index=False)

    per_replication = attempts[
        [
            "preflight_method", "semantic_replication_id", "dgp", "replication",
            "primary_status", "dgp_realization_hash", "supplied_rank_vector",
            "selected_rank_vector", "cap_pilot_rank", "candidate_coverage",
            "expected_fit_count", "replication_runtime_seconds",
        ]
    ].copy()
    per_replication.to_csv(output / "method_replications.csv", index=False)
    replication_diagnostics = pd.concat(
        [_replication_diagnostics(label, item) for label, item in frames.items()],
        ignore_index=True,
    )
    replication_diagnostics.to_csv(output / "replication_diagnostics.csv", index=False)
    replication_diagnostics.loc[
        replication_diagnostics["method"].eq("fixed_rank")
    ].to_csv(output / "fixed_rank_replications.csv", index=False)
    replication_diagnostics.loc[
        replication_diagnostics["method"].ne("fixed_rank")
    ].to_csv(output / "selected_rank_replications.csv", index=False)

    split_failures = []
    for label, item in frames.items():
        broad = item["records"].loc[item["records"]["corrected"].fillna(False).astype(bool)]
        successful = broad.loc[broad["primary_status"].eq("success")]
        bad = successful.loc[successful["split_coefficient_fit_count"].ne(4)]
        if len(bad):
            split_failures.append(label)
    if split_failures:
        raise ValueError(f"broad-target split coefficient-fit audit failed: {split_failures}")

    fixed_attempts = frames["fixed_rank"]["attempts"]
    fixed_records = frames["fixed_rank"]["records"]
    fixed_inference = frames["fixed_rank"]["inference"]
    audit = {
        "attempted_method_replications": len(attempts),
        "unique_dgp_draws": len(matching),
        "matched_hashes_pass": bool(matching["all_hashes_equal"].all()),
        "fixed_success": int(fixed_attempts["primary_status"].eq("success").sum()),
        "fixed_failures": fixed_attempts.loc[
            fixed_attempts["primary_status"].ne("success"), "primary_status"
        ].value_counts().to_dict(),
        "fixed_point_retained": int(_as_bool(fixed_records["retained_for_bias_rmse"]).sum()),
        "fixed_inference_retained": int(_as_bool(fixed_records["retained_for_coverage"]).sum()),
        "fixed_target_attempts": len(fixed_records),
        "fixed_gram_min": float(_finite(fixed_inference["tangent_gram_min_eigenvalue"]).min()),
        "fixed_gram_condition_max": float(
            _finite(fixed_inference["tangent_gram_condition_number"]).max()
        ),
        "fixed_riesz_failure_count": int((~_as_bool(fixed_inference["riesz_converged"])).sum()),
        "selected": selected_rows,
        "all_accounting_reconciled": bool(accounting["reconciled"].all()),
        "exactly_four_shared_split_fits_pass": True,
        "fit_runtime_available_in_executed_outputs": False,
        "go_no_go": "NO-GO: return evidence to author; do not launch medium or production",
    }
    (output / "audit.json").write_text(json.dumps(audit, indent=2) + "\n", encoding="utf-8")
    lines = [
        "# Final small independent preflight audit",
        "",
        "Preliminary diagnostics only; three replications per DGP are not substantive Monte Carlo evidence.",
        "",
        f"- Matched DGP hashes: PASS ({len(matching)}/12 semantic draws, all three methods equal).",
        f"- Attempt accounting: PASS ({len(attempts)} = 36 method-replication evaluations).",
        f"- Fixed rank: {audit['fixed_success']}/12 successful; failures {audit['fixed_failures']}.",
        f"- Fixed target retention: point {audit['fixed_point_retained']}/{len(fixed_records)}; inference {audit['fixed_inference_retained']}/{len(fixed_records)}.",
        f"- Selected 3e-6 pilot success: {selected_rows[0]['pilot_success']}/12.",
        f"- Selected 3e-6 ranks: {selected_rows[0]['selected_rank_distribution']}.",
        f"- Selected 4e-6 pilot success: {selected_rows[1]['pilot_success']}/12.",
        f"- Selected 4e-6 ranks: {selected_rows[1]['selected_rank_distribution']}.",
        "- Broad-target split reuse: PASS (exactly four coefficient fits per successful replication).",
        "- Fit-level runtime: UNAVAILABLE in these executed outputs; instrumentation was corrected only for future runs, and no extra evaluations were launched.",
        "- Decision: NO-GO because both pilots completed only 9/12 and fit-runtime reporting was incomplete; no medium or production run was launched.",
        "",
        "See the adjacent CSV and JSON files for lossless detail.",
    ]
    (output / "REPORT.md").write_text("\n".join(lines) + "\n", encoding="utf-8")
    print(json.dumps(audit, indent=2))


if __name__ == "__main__":
    main()
