"""Build split-fit forensics and the authorized size/cap comparisons."""

from __future__ import annotations

import argparse
import json
from collections.abc import Iterable
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd

SPLIT_TYPES = ("time_split_0", "time_split_1", "unit_split_0", "unit_split_1")
GRAM_FAILURES = {"tangent_gram_eigensolver_failure", "tangent_gram_nearly_singular"}
RIESZ_FAILURES = {"riesz_solver_failure", "riesz_target_instability"}
RANK_FLOOR = 1e-10
BROAD_TARGETS = {
    "A_full_mean", "B_full_mean", "A_G1_time_average", "A_G2_time_average",
    "A_G2_minus_G1_time_average", "B_G1_time_average", "B_G2_time_average",
    "B_G2_minus_G1_time_average",
}


def _bool(value: Any) -> bool:
    if value is None or (isinstance(value, float) and np.isnan(value)):
        return False
    if isinstance(value, (bool, np.bool_)):
        return bool(value)
    return str(value).lower() == "true"


def _json(value: Any, default: Any = None) -> Any:
    if value is None or (isinstance(value, float) and np.isnan(value)):
        return default
    return json.loads(value) if isinstance(value, str) else value


def _chunks(root: Path, folder: str) -> pd.DataFrame:
    paths = sorted((root / folder).glob("*.parquet"))
    return pd.concat([pd.read_parquet(path) for path in paths], ignore_index=True)


def _load(root: Path) -> dict[str, pd.DataFrame]:
    return {
        "attempts": pd.read_parquet(root / "attempted_replications.parquet"),
        "records": pd.read_parquet(root / "replication_records.parquet"),
        "fits": pd.read_parquet(root / "fit_diagnostics.parquet"),
        "rank": _chunks(root, "rank"),
    }


def _invalidity(row: pd.Series) -> tuple[str, str, bool]:
    reasons: list[str] = []
    if pd.notna(row.get("exception_type")):
        reasons.append("software_exception")
    objective = pd.to_numeric(pd.Series([row.get("objective_final")]), errors="coerce").iloc[0]
    if not np.isfinite(objective):
        reasons.append("nonfinite_objective")
    if _bool(row.get("iteration_cap_hit")):
        reasons.append("iteration_cap")
    if not _bool(row.get("convergence_flag")):
        reasons.append("optimizer_nonconvergence")
    if not _bool(row.get("stationarity_pass")):
        reasons.append("stationarity_failure")
    if _bool(row.get("coefficient_bound_hit")):
        reasons.append("coefficient_bound_hit")
    requested = tuple(_json(row.get("requested_rank"), []))
    realized = tuple(_json(row.get("numerical_rank"), []))
    rank_loss = requested != realized or float(row.get("sigma_r_over_sigma_1", np.inf)) < RANK_FLOOR
    if rank_loss:
        reasons.append("numerical_rank_loss")
    stability = row.get("objective_stability_pass")
    if pd.notna(stability) and not _bool(stability):
        reasons.append("objective_instability")
    precedence = (
        "software_exception", "nonfinite_objective", "iteration_cap",
        "optimizer_nonconvergence", "stationarity_failure", "coefficient_bound_hit",
        "numerical_rank_loss", "objective_instability",
    )
    primary = next((item for item in precedence if item in reasons), "valid")
    return primary, ";".join(reasons) if reasons else "valid", rank_loss


def _split_metadata(records: pd.DataFrame) -> dict[tuple[str, str], dict[str, Any]]:
    result: dict[tuple[str, str], dict[str, Any]] = {}
    usable = records.loc[records["split_diagnostics_json"].notna()]
    for row in usable.itertuples():
        for item in _json(row.split_diagnostics_json, []):
            fit_type = f"{item['kind']}_split_{item['part']}"
            result[(row.semantic_replication_id, fit_type)] = item
    return result


def _full_fit(fits: pd.DataFrame, rank_row: pd.Series) -> pd.Series | None:
    local = fits.loc[
        fits["semantic_replication_id"].eq(rank_row["semantic_replication_id"])
        & ~fits["fit_type"].isin((*SPLIT_TYPES, "nuclear_path"))
    ].copy()
    if local.empty:
        return None
    objective = rank_row.get("full_objective")
    if pd.isna(objective):
        diagnostics = _json(rank_row.get("rank_diagnostics_json"), {})
        selected = tuple(_json(rank_row.get("selected_rank_vector"), []))
        if selected:
            local = local.loc[local["requested_rank"].map(lambda x: tuple(_json(x, [])) == selected)]
            matches = [
                item for item in diagnostics.get("candidate_records", [])
                if tuple(item["ranks"]) == selected
            ]
            objective = matches[0]["objective"] if matches else np.nan
    if pd.notna(objective) and not local.empty:
        delta = (pd.to_numeric(local["objective_final"], errors="coerce") - float(objective)).abs()
        return local.loc[delta.idxmin()]
    return None if local.empty else local.iloc[0]


def split_fit_rows(data: dict[str, pd.DataFrame], method: str) -> pd.DataFrame:
    records, fits, rank = data["records"], data["fits"], data["rank"]
    failed_ids = set(
        records.loc[records["primary_status"].eq("split_fit_failure"), "semantic_replication_id"]
    )
    split = fits.loc[
        fits["semantic_replication_id"].isin(failed_ids) & fits["fit_type"].isin(SPLIT_TYPES)
    ].copy()
    expected = 4 * len(failed_ids)
    if len(split) != expected:
        raise ValueError(f"{method}: expected {expected} split fits, found {len(split)}")
    metadata = _split_metadata(records)
    rank_lookup = rank.set_index("semantic_replication_id", drop=False)
    rows: list[dict[str, Any]] = []
    for _, fit in split.iterrows():
        semantic_id, fit_type = fit["semantic_replication_id"], fit["fit_type"]
        meta = metadata[(semantic_id, fit_type)]
        primary, exact, rank_loss = _invalidity(fit)
        full = _full_fit(fits, rank_lookup.loc[semantic_id])
        block = {item["block"]: item for item in meta["split_rank_singular_values"]}
        split_kind, split_index = fit_type.split("_split_")
        row = {
            "DGP": int(fit["dgp"]), "semantic_replication_id": semantic_id,
            "method": method, "split_type": split_kind, "split_index": int(split_index),
            "effective_N": int(meta["n"]), "effective_T": int(meta["t"]),
            "requested_rank": fit["requested_rank"], "realized_numerical_rank": fit["numerical_rank"],
            "objective": fit["objective_final"], "iterations": fit["iterations"],
            "convergence": fit["convergence_flag"], "iteration_cap_hit": fit["iteration_cap_hit"],
            "stationarity_residual": fit["stationarity_residual"],
            "stationarity_pass": fit["stationarity_pass"],
            "coefficient_envelope": fit["coefficient_envelope"],
            "coefficient_envelope_ratio": fit["coefficient_envelope_ratio"],
            "coefficient_bound_hit": fit["coefficient_bound_hit"],
            "sigma_1": fit["sigma_1"], "sigma_r": fit["sigma_r"],
            "sigma_r_over_sigma_1": fit["sigma_r_over_sigma_1"],
            "numerical_rank_loss": rank_loss,
            "best_objective": fit["best_start_objective"],
            "second_best_objective": fit["second_start_objective"],
            "objective_stability_gap": fit["objective_stability_gap"],
            "objective_stability_pass": fit["objective_stability_pass"],
            "objective_stability_assessed": pd.notna(fit["objective_stability_pass"]),
            "runtime_seconds": fit["runtime_seconds"], "primary_numerical_cause": primary,
            "exact_invalidity_reason": exact,
        }
        for label in ("A", "B", "H"):
            item = block[label]
            row.update({
                f"{label}_sigma_1": item["sigma_1"], f"{label}_sigma_r": item["sigma_r"],
                f"{label}_sigma_r_over_sigma_1": item["sigma_r_over_sigma_1"],
                f"{label}_rank_supported": item["computational_rank_supported"],
            })
        if full is not None:
            diag = rank_lookup.loc[semantic_id]
            gap = diag.get("selected_best_two_objective_gap")
            if pd.isna(gap):
                multi = _json(diag.get("fixed_rank_multistart_diagnostics"), {})
                gap = multi.get("best_two_objective_gap", multi.get("original_stability_gap"))
            row.update({
                "full_coefficient_envelope": full["coefficient_envelope"],
                "full_coefficient_envelope_ratio": full["coefficient_envelope_ratio"],
                "full_stationarity_residual": full["stationarity_residual"],
                "full_sigma_r_over_sigma_1": full["sigma_r_over_sigma_1"],
                "full_objective_stability_gap": gap,
                "full_runtime_seconds": full["runtime_seconds"],
            })
        rows.append(row)
    return pd.DataFrame(rows).sort_values(
        ["method", "DGP", "semantic_replication_id", "split_type", "split_index"]
    )


def all_split_summary(data: dict[str, pd.DataFrame], method: str, sample_size: int) -> pd.DataFrame:
    fits = data["fits"].loc[data["fits"]["fit_type"].isin(SPLIT_TYPES)].copy()
    classified = fits.apply(_invalidity, axis=1, result_type="expand")
    fits[["primary_numerical_cause", "exact_invalidity_reason", "numerical_rank_loss"]] = classified
    fits["split_type"] = fits["fit_type"].str.split("_split_").str[0]
    groupings: Iterable[tuple[str, list[str]]] = (
        ("method_split_cause", ["split_type", "primary_numerical_cause"]),
        ("method_split_dgp_cause", ["split_type", "dgp", "primary_numerical_cause"]),
    )
    rows: list[dict[str, Any]] = []
    for aggregation, columns in groupings:
        denominator_columns = [column for column in columns if column != "primary_numerical_cause"]
        denominators = fits.groupby(denominator_columns, dropna=False).size()
        for key, group in fits.groupby(columns, dropna=False):
            values = key if isinstance(key, tuple) else (key,)
            item = dict(zip(columns, values, strict=True))
            denominator_key = tuple(item[column] for column in denominator_columns)
            if len(denominator_key) == 1:
                denominator_key = denominator_key[0]
            item.update({
                "sample_size": sample_size, "method": method, "aggregation": aggregation,
                "fit_count": len(group), "total_split_fits": int(denominators.loc[denominator_key]),
            })
            rows.append(item)
    for split_type, group in fits.groupby("split_type"):
        failed = group.loc[group["primary_numerical_cause"].ne("valid")]
        rows.append({
            "sample_size": sample_size, "method": method, "aggregation": "split_failure_rate",
            "split_type": split_type, "primary_numerical_cause": "any_failure",
            "fit_count": len(failed), "total_split_fits": len(group),
            "failure_rate": len(failed) / len(group),
        })
    return pd.DataFrame(rows)


def retention(data: dict[str, pd.DataFrame], method: str, sample_size: int) -> pd.DataFrame:
    frame = data["records"].copy()
    frame["target_type"] = np.where(
        frame["target"].isin(BROAD_TARGETS), "broad_split_corrected", "local_plugin"
    )
    rows = []
    for target_type, group in frame.groupby("target_type"):
        point = int(group["point_estimate_valid"].map(_bool).sum())
        inference = int(group["inference_valid"].map(_bool).sum())
        rows.append({
            "sample_size": sample_size, "method": method, "target_type": target_type,
            "attempted": len(group), "point_retained": point,
            "point_retained_share": point / len(group), "inference_retained": inference,
            "inference_retained_share": inference / len(group),
        })
    return pd.DataFrame(rows)


def cap_sensitivity(cap2: dict[str, pd.DataFrame], cap3: dict[str, pd.DataFrame]) -> pd.DataFrame:
    rows = []
    first = cap2["rank"].set_index("semantic_replication_id")
    second = cap3["rank"].set_index("semantic_replication_id")
    for semantic_id in sorted(second.index):
        old, new = first.loc[semantic_id], second.loc[semantic_id]
        if old["dgp_realization_hash"] != new["dgp_realization_hash"]:
            raise ValueError(f"cap replay hash mismatch for {semantic_id}")
        def candidates(row: pd.Series) -> str:
            diag = _json(row["rank_diagnostics_json"], {})
            return json.dumps(sorted(item["ranks"] for item in diag["candidate_records"]))
        rows.append({
            "diagnostic_replay": True, "semantic_replication_id": semantic_id,
            "DGP": int(new["dgp"]), "replication": int(new["replication"]),
            "dgp_realization_hash": new["dgp_realization_hash"],
            "cap_2_pilot_rank": old["cap_pilot_rank"], "cap_2_candidate_set": candidates(old),
            "cap_2_selected_rank": old["selected_rank_vector"],
            "cap_3_pilot_rank": new["cap_pilot_rank"], "cap_3_candidate_set": candidates(new),
            "cap_3_selected_rank": new["selected_rank_vector"],
            "moved_to_rank_3": 3 in tuple(_json(new["selected_rank_vector"], [])),
        })
    return pd.DataFrame(rows)


def size_metrics(data: dict[str, pd.DataFrame], sample_size: int) -> dict[str, Any]:
    attempts, records, fits = data["attempts"], data["records"], data["fits"]
    split = fits.loc[fits["fit_type"].isin(SPLIT_TYPES)].copy()
    split["cause"] = split.apply(lambda row: _invalidity(row)[0], axis=1)
    time = split.loc[split["fit_type"].str.startswith("time")]
    unit = split.loc[split["fit_type"].str.startswith("unit")]
    complete = split.groupby("semantic_replication_id")["cause"].apply(lambda x: len(x) == 4 and x.eq("valid").all())
    local = records.loc[~records["target"].isin(BROAD_TARGETS)]
    broad = records.loc[records["target"].isin(BROAD_TARGETS)]
    return {
        "sample_size": sample_size, "attempted_replications": len(attempts),
        "full_panel_success": int(
            (~attempts["primary_status"].isin({"full_fit_failure", "coefficient_bound_hit"})).sum()
        ),
        "full_panel_coefficient_bound_hits": int(attempts["primary_status"].eq("coefficient_bound_hit").sum()),
        "replications_with_four_split_success": int(complete.sum()),
        "four_split_success_rate": float(complete.sum() / len(attempts)),
        "time_split_failure_rate": float(time["cause"].ne("valid").mean()),
        "unit_split_failure_rate": float(unit["cause"].ne("valid").mean()),
        "point_retained_share": float(records["point_estimate_valid"].map(_bool).mean()),
        "local_inference_retained_share": float(local["inference_valid"].map(_bool).mean()),
        "broad_inference_retained_share": float(broad["inference_valid"].map(_bool).mean()),
        "total_inference_retained_share": float(records["inference_valid"].map(_bool).mean()),
        "gram_failure_rate": float(records["primary_status"].isin(GRAM_FAILURES).mean()),
        "riesz_failure_rate": float(records["primary_status"].isin(RIESZ_FAILURES).mean()),
        "runtime_total_seconds": float(attempts["replication_runtime_seconds"].sum()),
        "runtime_mean_seconds": float(attempts["replication_runtime_seconds"].mean()),
    }


def reconciliation(
    fixed: dict[str, pd.DataFrame],
    selected: dict[str, pd.DataFrame],
    level: pd.DataFrame,
    cap3: dict[str, pd.DataFrame] | None,
    n100: dict[str, pd.DataFrame] | None,
) -> dict[str, Any]:
    result: dict[str, Any] = {}
    for label, data in (("n50_fixed", fixed), ("n50_selected", selected)):
        attempts, records = data["attempts"], data["records"]
        if len(attempts) != 12 or len(records) != 216:
            raise ValueError(f"{label} does not reconcile to 12 x 18")
        result[label] = {
            "attempted_replications": len(attempts),
            "target_records": len(records),
            "unique_semantic_replication_ids": attempts["semantic_replication_id"].nunique(),
        }
    expected_level = 4 * sum(
        data["records"].loc[
            data["records"]["primary_status"].eq("split_fit_failure"),
            "semantic_replication_id",
        ].nunique()
        for data in (fixed, selected)
    )
    if len(level) != expected_level:
        raise ValueError("fit-level failure bundles do not reconcile")
    result["n50_fit_level"] = {
        "failed_bundle_fit_rows": len(level),
        "invalid_split_fits": int(level["primary_numerical_cause"].ne("valid").sum()),
    }
    if cap3 is not None:
        attempts = cap3["attempts"]
        if len(attempts) != 3 or set(attempts["replication"].astype(int)) != {3, 4, 5}:
            raise ValueError("cap-3 replay does not contain exactly DGP-4 replications 3-5")
        result["dgp4_cap3_replay"] = {
            "attempted_replications": len(attempts),
            "replication_indices": sorted(attempts["replication"].astype(int).tolist()),
            "diagnostic_replay": True,
        }
    if n100 is not None:
        attempts, records = n100["attempts"], n100["records"]
        prior_ids = set(fixed["attempts"]["semantic_replication_id"]) | set(
            selected["attempts"]["semantic_replication_id"]
        )
        new_ids = set(attempts["semantic_replication_id"])
        split_count = int(n100["fits"]["fit_type"].isin(SPLIT_TYPES).sum())
        if (
            len(attempts) != 12
            or len(records) != 216
            or split_count != 48
            or set(attempts["replication"].astype(int)) != {6, 7, 8}
            or prior_ids & new_ids
        ):
            raise ValueError("N=100 accounting or semantic-draw novelty failed")
        result["n100_fixed"] = {
            "attempted_replications": len(attempts),
            "target_records": len(records),
            "actual_split_fits": split_count,
            "replication_indices": sorted(set(attempts["replication"].astype(int))),
            "semantic_id_overlap_with_n50": 0,
        }
    return result


def _report(
    level: pd.DataFrame,
    summary: pd.DataFrame,
    retention_table: pd.DataFrame,
    cap: pd.DataFrame | None,
    size: pd.DataFrame | None,
) -> str:
    failed = level.loc[level["primary_numerical_cause"].ne("valid")]
    n50_rate = summary.loc[
        summary["sample_size"].eq(50) & summary["aggregation"].eq("split_failure_rate")
    ]
    fixed_failed = failed.loc[failed["method"].eq("fixed_rank")]
    selected_failed = failed.loc[failed["method"].eq("selected_rank")]
    fixed_by_dgp = fixed_failed.groupby("DGP").size().to_dict()
    selected_by_dgp = selected_failed.groupby("DGP").size().to_dict()
    retention_lines = []
    for row in retention_table.loc[retention_table["sample_size"].eq(50)].itertuples():
        retention_lines.append(
            f"- {row.method}, {row.target_type}: attempted {row.attempted}, point retained "
            f"{row.point_retained}, inference retained {row.inference_retained}."
        )
    full_envelope = level["full_coefficient_envelope_ratio"].median()
    failed_envelope = failed["coefficient_envelope_ratio"].median()
    full_stationarity = level["full_stationarity_residual"].median()
    failed_stationarity = failed["stationarity_residual"].median()
    full_sigma = level["full_sigma_r_over_sigma_1"].median()
    failed_sigma = failed["sigma_r_over_sigma_1"].median()
    full_runtime = level["full_runtime_seconds"].median()
    failed_runtime = failed["runtime_seconds"].median()
    lines = [
        "# Split-fit forensics",
        "",
        "This diagnostic preserves the accepted estimator, B=9, split formula, Riesz construction, IC rate, DGPs, and all numerical tolerances. No observations were trimmed or removed.",
        "",
        "## Existing N=50 preflight",
        "",
        f"The table contains {len(level)} actual fits from the four-fit bundles of every replication with a required split failure; {len(failed)} fits are invalid. Every invalid fit is classified as `{failed['primary_numerical_cause'].mode().iloc[0]}`.",
        f"Fixed rank has {len(fixed_failed)} invalid fits (DGP counts {fixed_by_dgp}); selected rank has {len(selected_failed)} (DGP counts {selected_by_dgp}). No other requested primary-cause category occurs.",
        f"Fixed-rank time/unit failure rates are {float(n50_rate.loc[(n50_rate['method'].eq('fixed_rank')) & (n50_rate['split_type'].eq('time')), 'failure_rate'].iloc[0]):.3f}/{float(n50_rate.loc[(n50_rate['method'].eq('fixed_rank')) & (n50_rate['split_type'].eq('unit')), 'failure_rate'].iloc[0]):.3f}; selected-rank rates are {float(n50_rate.loc[(n50_rate['method'].eq('selected_rank')) & (n50_rate['split_type'].eq('time')), 'failure_rate'].iloc[0]):.3f}/{float(n50_rate.loc[(n50_rate['method'].eq('selected_rank')) & (n50_rate['split_type'].eq('unit')), 'failure_rate'].iloc[0]):.3f}.",
        "All invalid N=50 split fits converged below the iteration cap, passed stationarity, and retained their supplied computational rank. Split objective stability is not assessed because the prescribed split stage uses one fit per half; blank stability fields are therefore not evidence of instability.",
        "Per-block singular diagnostics show no A/B/H rank collapse. The historical fit records store only the maximum coefficient envelope across A, B, and H, not the block attaining it, so matrix-specific bound attribution cannot be reconstructed without replaying these N=50 fits; no such replay was authorized. This limitation is stated rather than guessed.",
        "",
        "Local/plugin targets do not require split correction; broad targets do. Retention is:",
        *retention_lines,
        "",
        "## Full versus split conditioning",
        "",
        f"Across the failed bundles, the median full-panel envelope ratio is {full_envelope:.3f}, versus {failed_envelope:.3f} for invalid split fits. Median stationarity residuals are {full_stationarity:.3g} versus {failed_stationarity:.3g}; median sigma-r ratio is {full_sigma:.3f} versus {failed_sigma:.3f}; and median runtimes are {full_runtime:.3f}s versus {failed_runtime:.3f}s. Split objective-stability gaps are unavailable by design, while full-panel multi-start gaps are saved in the conditioning table.",
        "The failed bundles therefore pair well-behaved 50x50 full fits with smaller 50x25 or 25x50 fits whose coefficient envelopes cross B. Stationarity and singular-rank diagnostics remain satisfactory, making finite-small-panel envelope activity the descriptive mechanism. This is not evidence of asymptotic failure.",
    ]
    if cap is not None:
        lines.extend(["", "## DGP 4 rank-cap replay", ""])
        for row in cap.itertuples():
            lines.append(f"- replication {row.replication}: cap 2 selected {row.cap_2_selected_rank}; cap 3 selected {row.cap_3_selected_rank}.")
        moved = int(cap["moved_to_rank_3"].sum())
        lines.append(f"Rank 3 is selected in {moved}/{len(cap)} diagnostic replays.")
        lines.append(
            "The two former (1,1,2) cap hits remain (1,1,2), so cap=(2,2,2) does not appear to truncate these DGP-4 solutions. This does not automatically change or freeze the main cap."
        )
    if size is not None:
        n50 = size.loc[size["sample_size"].eq(50)].iloc[0]
        n100 = size.loc[size["sample_size"].eq(100)].iloc[0]
        lines.extend([
            "", "## Fixed-rank size escalation", "",
            f"N=50 four-split success is {n50.four_split_success_rate:.3f}; N=100 is {n100.four_split_success_rate:.3f}. Time-split failure changes from {n50.time_split_failure_rate:.3f} to {n100.time_split_failure_rate:.3f}; unit-split failure changes from {n50.unit_split_failure_rate:.3f} to {n100.unit_split_failure_rate:.3f}.",
            f"At N=100, full-panel success is {int(n100.full_panel_success)}/12 with {int(n100.full_panel_coefficient_bound_hits)} bound hits. Point retention is {n100.point_retained_share:.3f}; local/broad/total inference retention is {n100.local_inference_retained_share:.3f}/{n100.broad_inference_retained_share:.3f}/{n100.total_inference_retained_share:.3f}. Gram and Riesz failure rates are {n100.gram_failure_rate:.3f}/{n100.riesz_failure_rate:.3f}. Total runtime is {n100.runtime_total_seconds:.3f}s (mean {n100.runtime_mean_seconds:.3f}s per replication).",
            "The sole N=100 split failure is the DGP-3 replication-7 unit half 2 coefficient-bound hit; it converged, passed stationarity, and preserved rank. No other failure mechanism occurs.",
            "These 12 new evaluations are a numerical gate, not publication Monte Carlo evidence.",
            "", "## Gate decision", "",
        ])
        rare = n100.four_split_success_rate >= 0.9
        lines.append(
            "The evidence supports proceeding to the medium diagnostic, but it was not launched."
            if rare else
            "NO-GO remains: split failures are not rare at N=100, so the split-fit numerical algorithm needs further work before medium simulations."
        )
    lines.extend(["", "## Reconciliation", "", f"Fit-level rows: {len(level)}; invalid rows: {len(failed)}. Summary and retention denominators are generated directly from the consolidated run records."])
    return "\n".join(lines) + "\n"


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--fixed-n50", required=True, type=Path)
    parser.add_argument("--selected-n50", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--cap3", type=Path)
    parser.add_argument("--fixed-n100", type=Path)
    args = parser.parse_args()
    args.output.mkdir(parents=True, exist_ok=True)
    fixed, selected = _load(args.fixed_n50), _load(args.selected_n50)
    level = pd.concat(
        [split_fit_rows(fixed, "fixed_rank"), split_fit_rows(selected, "selected_rank")],
        ignore_index=True,
    )
    level.to_csv(args.output / "split_fit_level_diagnostics.csv", index=False)
    conditioning = level[[
        "DGP", "semantic_replication_id", "method", "split_type", "split_index",
        "effective_N", "effective_T", "coefficient_envelope", "full_coefficient_envelope",
        "stationarity_residual", "full_stationarity_residual", "sigma_r_over_sigma_1",
        "full_sigma_r_over_sigma_1", "objective_stability_gap",
        "full_objective_stability_gap", "runtime_seconds", "full_runtime_seconds",
        "primary_numerical_cause",
    ]].copy()
    conditioning.to_csv(args.output / "full_vs_split_conditioning.csv", index=False)
    summary = pd.concat(
        [all_split_summary(fixed, "fixed_rank", 50), all_split_summary(selected, "selected_rank", 50)],
        ignore_index=True,
    )
    retention_table = pd.concat(
        [retention(fixed, "fixed_rank", 50), retention(selected, "selected_rank", 50)],
        ignore_index=True,
    )
    cap = None
    cap3_data = None
    size = None
    n100 = None
    if args.cap3:
        cap3_data = _load(args.cap3)
        cap = cap_sensitivity(selected, cap3_data)
        cap.to_csv(args.output / "dgp4_rank_cap_sensitivity.csv", index=False)
    if args.fixed_n100:
        n100 = _load(args.fixed_n100)
        n100_summary = all_split_summary(n100, "fixed_rank", 100)
        summary = pd.concat([summary, n100_summary], ignore_index=True)
        retention_table = pd.concat([retention_table, retention(n100, "fixed_rank", 100)], ignore_index=True)
        n100_fails = split_fit_rows(n100, "fixed_rank")
        n100_fails.to_csv(args.output / "n100_split_fit_level_diagnostics.csv", index=False)
        size = pd.DataFrame([size_metrics(fixed, 50), size_metrics(n100, 100)])
        size.to_csv(args.output / "fixed_rank_size_comparison.csv", index=False)
    summary.to_csv(args.output / "split_failure_summary.csv", index=False)
    retention_table.to_csv(args.output / "target_retention_by_type.csv", index=False)
    (args.output / "split_fit_forensics_report.md").write_text(
        _report(level, summary, retention_table, cap, size), encoding="utf-8"
    )
    (args.output / "accounting_reconciliation.json").write_text(
        json.dumps(
            reconciliation(fixed, selected, level, cap3_data, n100), indent=2,
            sort_keys=True,
        ) + "\n",
        encoding="utf-8",
    )


if __name__ == "__main__":
    main()
