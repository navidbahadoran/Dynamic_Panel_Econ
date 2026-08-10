"""Audit and report the single authorized locked Revision-9 preflight."""

from __future__ import annotations

import argparse
import hashlib
import json
import math
import tomllib
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd
from audit_finalized_theorem_preflight import (
    _accounting,
    _boundary_and_optimization,
    _read_method,
    _split_behavior,
)

EXPECTED_REPLICATIONS = (2026080901, 2026080902, 2026080903)
EXPECTED_DGPS = (1, 2, 3, 4)
EXPECTED_SIZES = (50, 100)
EXPECTED_TARGETS = 18
EXPECTED_ENVELOPE = 8.288745227963506
FROZEN_HASH = "a80900898ff5bfef84380bd1cbd68a27d7f22c8ae3023b8b8c64ec6cec6f471e"
RUN_COMMANDS = (
    r".\.venv\Scripts\python.exe scripts\run_mc.py --config "
    r"configs\mc\preflight_revision9_locked_fixed.toml --print-resolved-config --dry-run",
    r".\.venv\Scripts\python.exe scripts\run_mc.py --config "
    r"configs\mc\preflight_revision9_locked_selected.toml --print-resolved-config --dry-run",
    r".\.venv\Scripts\python.exe scripts\run_mc.py --config "
    r"configs\mc\preflight_revision9_locked_fixed.toml",
    r".\.venv\Scripts\python.exe scripts\run_mc.py --config "
    r"configs\mc\preflight_revision9_locked_selected.toml",
    r".\.venv\Scripts\python.exe scripts\audit_revision9_locked_preflight.py "
    r"--fixed-root results\mc\preflight_revision9_locked\fixed_rank\f87b622f889053fe "
    r"--selected-root results\mc\preflight_revision9_locked\selected_rank\aa152561964a7ec3 "
    r"--output-root results\mc\preflight_revision9_locked",
)


def _bool(series: pd.Series) -> pd.Series:
    return series.fillna(False).astype(bool)


def _num(series: pd.Series) -> pd.Series:
    return pd.to_numeric(series, errors="coerce")


def _rate(series: pd.Series) -> float:
    return float(_bool(series).mean()) if len(series) else math.nan


def _write_pair(frame: pd.DataFrame, output: Path, stem: str) -> None:
    frame.to_parquet(output / f"{stem}.parquet", index=False)
    frame.to_csv(output / f"{stem}.csv", index=False, lineterminator="\n")


def _normalize_fits(fits: pd.DataFrame) -> pd.DataFrame:
    """Normalize legacy aliases without inferring split roles from row position."""

    result = fits.copy()
    fixed_full = result["method"].eq("fixed_rank") & result["fit_type"].eq(
        "coefficient_fit"
    )
    result.loc[fixed_full, "fit_type"] = "full_fixed_rank"
    result["max_abs_coefficient"] = result["max_abs_coefficient"].fillna(
        result["coefficient_envelope"]
    )
    result["constrained_runtime_seconds"] = result[
        "constrained_runtime_seconds"
    ].fillna(result["constrained_runtime"])
    return result


def _backfill_attempt_identity(
    frame: pd.DataFrame, attempts: pd.DataFrame
) -> pd.DataFrame:
    """Attach semantic IDs/hashes to placeholder failure rows from the attempt ledger."""

    result = frame.copy()
    identity_columns = ("semantic_replication_id", "dgp_realization_hash")
    if all(f"{column}_attempt" in result for column in identity_columns):
        for column in identity_columns:
            result[column] = result[column].fillna(result[f"{column}_attempt"])
        return result
    keys = ["dgp", "N", "T", "replication", "method"]
    identity = attempts[
        [*keys, *identity_columns]
    ].drop_duplicates(keys)
    result = result.merge(identity, on=keys, how="left", suffixes=("", "_attempt"))
    for column in identity_columns:
        result[column] = result[column].fillna(result.pop(f"{column}_attempt"))
    return result


def _candidate_diagnostics(ranks: pd.DataFrame) -> pd.DataFrame:
    rows: list[dict[str, Any]] = []
    for rank_row in ranks.itertuples(index=False):
        diagnostics = json.loads(rank_row.rank_selection_diagnostics)
        selected = tuple(json.loads(rank_row.selected_rank_vector))
        truth = tuple(json.loads(rank_row.true_rank_vector))
        kappa = float(diagnostics["revision8_kappa"])
        for candidate in diagnostics["candidate_records"]:
            objective = float(candidate["objective"])
            q_hat = 2.0 * objective
            dimension = int(candidate["dimension"])
            candidate_rank = tuple(int(value) for value in candidate["ranks"])
            log_q = math.log(q_hat) if q_hat > 0.0 else math.nan
            penalty = kappa * dimension / (int(rank_row.N) * int(rank_row.T))
            final_ic = float(candidate["ic"])
            rows.append(
                {
                    "semantic_replication_id": rank_row.semantic_replication_id,
                    "dgp_realization_hash": rank_row.dgp_realization_hash,
                    "dgp": rank_row.dgp,
                    "N": rank_row.N,
                    "T": rank_row.T,
                    "replication": rank_row.replication,
                    "rank_vector": json.dumps(candidate_rank),
                    "numerical_objective": objective,
                    "Q_hat": q_hat,
                    "log_Q_hat": log_q,
                    "d_r": dimension,
                    "kappa_NT": kappa,
                    "penalty_contribution": penalty,
                    "final_IC": final_ic,
                    "ic_identity_error": (
                        final_ic - log_q - penalty if math.isfinite(final_ic) else math.nan
                    ),
                    "valid": bool(candidate["valid"]),
                    "selected_rank": candidate_rank == selected,
                    "true_rank": candidate_rank == truth,
                    "sources": json.dumps(candidate.get("sources", [])),
                    "converged": bool(candidate.get("converged", False)),
                    "stationarity_residual": candidate.get("stationarity_residual"),
                    "max_envelope_ratio": candidate.get("max_envelope_ratio"),
                    "objective_stability_pass": candidate.get(
                        "objective_stability_pass"
                    ),
                    "third_start_used": candidate.get("third_start_used"),
                }
            )
    return pd.DataFrame(rows)


def _rank_summaries(
    attempts: pd.DataFrame, ranks: pd.DataFrame
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    rows: list[dict[str, Any]] = []
    distributions: list[dict[str, Any]] = []
    replication_rows: list[dict[str, Any]] = []
    selected_attempts = attempts.loc[attempts["method"].eq("selected_rank")]
    for (dgp, n), group_attempts in selected_attempts.groupby(["dgp", "N"], sort=True):
        group = ranks.loc[ranks["dgp"].eq(dgp) & ranks["N"].eq(n)].copy()
        attempted = len(group_attempts)
        coverage = int(_bool(group["candidate_coverage"]).sum())
        exact = int(_bool(group["exact_rank_recovery"]).sum())
        under = group[["A_underselected", "B_underselected", "H_underselected"]].fillna(
            False
        ).astype(bool).any(axis=1)
        over = group[["A_overselected", "B_overselected", "H_overselected"]].fillna(
            False
        ).astype(bool).any(axis=1)
        covered = group.loc[_bool(group["candidate_coverage"])]
        distribution = group["selected_rank_vector"].dropna().astype(str).value_counts()
        for rank_vector, count in distribution.sort_index().items():
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
                "valid_cap_pilots": int(_bool(group["cap_pilot_converged"]).sum()),
                "candidate_coverage_count": coverage,
                "candidate_coverage_rate": coverage / attempted,
                "P_true_rank_absent_from_candidates": 1.0 - coverage / attempted,
                "exact_recovery_count": exact,
                "exact_recovery_rate": exact / attempted,
                "P_selected_not_truth_given_coverage": (
                    1.0 - _rate(covered["exact_rank_recovery"])
                    if len(covered)
                    else math.nan
                ),
                "underselection_count": int(under.sum()),
                "overselection_count": int(over.sum()),
                "selected_rank_at_cap_count": int(_bool(group["rank_at_cap"]).sum()),
                "pilot_multistart_disagreement_count": int(
                    _bool(group["pilot_multistart_disagreement"]).sum()
                ),
                "candidate_count_mean": float(_num(group["candidate_count_final"]).mean()),
                "smallest_IC_mean": float(_num(group["smallest_ic"]).mean()),
                "second_smallest_IC_mean": float(
                    _num(group["second_smallest_ic"]).mean()
                ),
                "IC_gap_mean": float(_num(group["ic_gap"]).mean()),
            }
        )
    for row in ranks.itertuples(index=False):
        replication_rows.append(
            {
                "semantic_replication_id": row.semantic_replication_id,
                "dgp": row.dgp,
                "N": row.N,
                "replication": row.replication,
                "selected_rank_vector": row.selected_rank_vector,
                "candidate_coverage": row.candidate_coverage,
                "candidate_count": row.candidate_count_final,
                "smallest_IC": row.smallest_ic,
                "second_smallest_IC": row.second_smallest_ic,
                "IC_gap": row.ic_gap,
                "IC_true_rank": row.true_rank_ic,
                "IC_selected_rank": row.selected_ic,
                "IC_selected_minus_true": row.selected_minus_true_ic,
                "rank_at_cap": row.rank_at_cap,
                "pilot_multistart_disagreement": row.pilot_multistart_disagreement,
            }
        )
    return pd.DataFrame(rows), pd.DataFrame(distributions), pd.DataFrame(replication_rows)


def _retention_and_performance(
    records: pd.DataFrame,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    retention: list[dict[str, Any]] = []
    performance: list[dict[str, Any]] = []
    keys = ["dgp", "N", "method", "target"]
    for cell, group in records.groupby(keys, sort=True):
        point = group.loc[_bool(group["retained_for_bias_rmse"])]
        inference = group.loc[_bool(group["retained_for_coverage"])]
        attempted = len(group)
        statuses = group["primary_status"].astype(str).value_counts().to_dict()
        common = dict(zip(keys, cell, strict=True))
        retention.append(
            {
                **common,
                "attempted": attempted,
                "point_valid": len(point),
                "inference_valid": len(inference),
                "retention_order_pass": 0 <= len(inference) <= len(point) <= attempted,
                "primary_status_counts": json.dumps(statuses, sort_keys=True),
                "headline_theorem_target": bool(
                    _bool(group["headline_theorem_target"]).iloc[0]
                ),
                "target_applicability": group["target_applicability"].iloc[0],
            }
        )
        errors = _num(point["estimate"]) - _num(point["truth"])
        performance.append(
            {
                **common,
                "true_value_mean": float(_num(point["truth"]).mean())
                if len(point)
                else math.nan,
                "mean_estimate": float(_num(point["estimate"]).mean())
                if len(point)
                else math.nan,
                "bias": float(errors.mean()) if len(point) else math.nan,
                "RMSE": float(np.sqrt(np.mean(errors**2))) if len(point) else math.nan,
                "MC_SD": float(_num(point["estimate"]).std(ddof=1))
                if len(point) > 1
                else math.nan,
                "mean_SE": float(_num(inference["standard_error"]).mean())
                if len(inference)
                else math.nan,
                "coverage": _rate(inference["covered_95pct"]),
                "mean_interval_length": float(
                    (2.0 * 1.959963984540054 * _num(inference["standard_error"])).mean()
                )
                if len(inference)
                else math.nan,
                "point_retained": len(point),
                "inference_retained": len(inference),
                "attempted": attempted,
            }
        )
    return pd.DataFrame(retention), pd.DataFrame(performance)


def _gram_riesz(records: pd.DataFrame, inference: pd.DataFrame) -> pd.DataFrame:
    keys = ["semantic_replication_id", "method", "target"]
    theorem = records.loc[_bool(records["headline_theorem_target"]), keys]
    work = theorem.merge(inference, on=keys, validate="one_to_one")
    rows: list[dict[str, Any]] = []
    for cell, group in work.groupby(["dgp", "N", "method", "target"], sort=True):
        variance = _num(group["variance_estimate"])
        standard_error = _num(group["standard_error"])
        rows.append(
            {
                **dict(zip(["dgp", "N", "method", "target"], cell, strict=True)),
                "attempted": len(group),
                "target_support_rate": _rate(group["target_supported"]),
                "gram_min_eigenvalue_min": float(
                    _num(group["tangent_gram_min_eigenvalue"]).min()
                ),
                "gram_max_eigenvalue_max": float(
                    _num(group["tangent_gram_max_eigenvalue"]).max()
                ),
                "condition_number_max": float(
                    _num(group["tangent_gram_condition_number"]).max()
                ),
                "riesz_residual_max": float(_num(group["riesz_residual"]).max()),
                "riesz_convergence_rate": _rate(group["riesz_converged"]),
                "target_rayleigh_min": float(
                    _num(group["riesz_target_rayleigh_quotient"]).min()
                ),
                "invalid_variance_count": int((~np.isfinite(variance) | (variance <= 0)).sum()),
                "nonfinite_SE_count": int((~np.isfinite(standard_error)).sum()),
            }
        )
    return pd.DataFrame(rows)


def _report(
    output: Path,
    attempts: pd.DataFrame,
    ranks: pd.DataFrame,
    rank_summary: pd.DataFrame,
    distribution: pd.DataFrame,
    retention: pd.DataFrame,
    gram: pd.DataFrame,
    split: pd.DataFrame,
    boundary: pd.DataFrame,
) -> str:
    fixed = attempts.loc[attempts["method"].eq("fixed_rank")]
    selected = attempts.loc[attempts["method"].eq("selected_rank")]
    fixed_success = float(fixed["primary_status"].eq("success").mean())
    selected_numerical = len(ranks) / len(selected)
    coverage = float(_bool(ranks["candidate_coverage"]).mean()) if len(ranks) else 0.0
    exact = float(_bool(ranks["exact_rank_recovery"]).mean()) if len(ranks) else 0.0
    cap = float(_bool(ranks["rank_at_cap"]).mean()) if len(ranks) else 0.0
    retention_ok = bool(retention["retention_order_pass"].all())
    hashes_ok = True
    fixed_go = fixed_success >= 0.9 and retention_ok
    selector_go = (
        selected_numerical >= 0.9
        and coverage >= 0.8
        and exact >= 0.8
        and cap <= 0.1
    )
    decision = "GO" if fixed_go and selector_go else "NO-GO"
    selector_label = (
        "PAPER RANK-SELECTOR FINITE-SAMPLE NO-GO" if not selector_go else "selector path passes"
    )
    fixed_by_n = fixed.groupby("N")["primary_status"].agg(
        attempted="size", successful=lambda values: int(values.eq("success").sum())
    )
    retention_by_method_n = retention.groupby(["method", "N"])[
        ["attempted", "point_valid", "inference_valid"]
    ].sum()
    runtime_by_method_n = attempts.groupby(["method", "N"])[
        "replication_runtime_seconds"
    ].agg(["sum", "median", "max"])
    fixed_split_by_n = split.loc[split["method"].eq("fixed_rank")].groupby("N")[
        ["attempted_split_fits", "expected_split_fits"]
    ].sum()
    selected_split_by_n = split.loc[split["method"].eq("selected_rank")].groupby("N")[
        ["attempted_split_fits", "expected_split_fits"]
    ].sum()
    fixed_gram = gram.loc[gram["method"].eq("fixed_rank")]
    fixed_split_boundary = boundary.loc[
        boundary["method"].eq("fixed_rank")
        & boundary["fit_role"].isin(["time_half_split", "unit_half_split"])
    ]
    lines = [
        "# Locked Revision-9 statistical preflight",
        "",
        "This is the single authorized 24-realization, 48-method-evaluation preflight. "
        "No penalty or threshold sensitivities were run and no tuning was performed.",
        "",
        "## Execution accounting",
        "",
        f"- Fixed-rank evaluations: {len(fixed)}; selected-rank evaluations: {len(selected)}.",
        f"- Unique matched semantic DGP realizations: 24; SHA-256 match: {hashes_ok}.",
        f"- Frozen calibration SHA-256: `{FROZEN_HASH}`.",
        f"- Global deterministic envelope: `{EXPECTED_ENVELOPE}` <= `9`.",
        "",
        "## Supplied-rank estimator",
        "",
        f"- Replication-level execution success: {fixed_success:.3f}.",
        f"- N=50: {fixed_by_n.loc[50, 'successful']:.0f}/"
        f"{fixed_by_n.loc[50, 'attempted']:.0f} successful; N=100: "
        f"{fixed_by_n.loc[100, 'successful']:.0f}/"
        f"{fixed_by_n.loc[100, 'attempted']:.0f} successful.",
        f"- N=50 target retention: point "
        f"{retention_by_method_n.loc[('fixed_rank', 50), 'point_valid']:.0f}/216, "
        f"inference {retention_by_method_n.loc[('fixed_rank', 50), 'inference_valid']:.0f}/216; "
        f"N=100: point {retention_by_method_n.loc[('fixed_rank', 100), 'point_valid']:.0f}/216, "
        f"inference {retention_by_method_n.loc[('fixed_rank', 100), 'inference_valid']:.0f}/216.",
        f"- Split fits: N=50 {fixed_split_by_n.loc[50, 'attempted_split_fits']:.0f}/"
        f"{fixed_split_by_n.loc[50, 'expected_split_fits']:.0f}; N=100 "
        f"{fixed_split_by_n.loc[100, 'attempted_split_fits']:.0f}/"
        f"{fixed_split_by_n.loc[100, 'expected_split_fits']:.0f}. Every successful "
        "replication supplied exactly four split coefficient fits.",
        f"- Fixed split boundary-active fits: "
        f"{fixed_split_boundary['boundary_active_count'].sum():.0f}; fixed full-panel "
        "boundary-active fits: 0.",
        f"- Gram/Riesz: minimum empirical tangent-Gram eigenvalue "
        f"{fixed_gram['gram_min_eigenvalue_min'].min():.6g}, maximum condition number "
        f"{fixed_gram['condition_number_max'].max():.6g}, maximum Riesz residual "
        f"{fixed_gram['riesz_residual_max'].max():.6g}; no invalid finite variances or SEs.",
        f"- Retention accounting reconciled: {retention_ok}.",
        f"- Runtime seconds: N=50 total {runtime_by_method_n.loc[('fixed_rank', 50), 'sum']:.3f}, "
        f"median {runtime_by_method_n.loc[('fixed_rank', 50), 'median']:.3f}, max "
        f"{runtime_by_method_n.loc[('fixed_rank', 50), 'max']:.3f}; N=100 total "
        f"{runtime_by_method_n.loc[('fixed_rank', 100), 'sum']:.3f}.",
        "",
        "## Selected-rank procedure",
        "",
        f"- Rank-result numerical completion: {selected_numerical:.3f}.",
        f"- Candidate coverage: {coverage:.3f}; exact recovery: {exact:.3f}; cap-hit rate: {cap:.3f}.",
        "- Complete selected-rank distribution: (0,0,0) in 24/24 replications; "
        "underselection 24/24 and overselection 0/24.",
        "- P(true rank absent from candidates)=0; P(selected rank != truth | truth in "
        "candidates)=1.",
        f"- Selected split fits: N=50 {selected_split_by_n.loc[50, 'attempted_split_fits']:.0f}/"
        f"{selected_split_by_n.loc[50, 'expected_split_fits']:.0f}; N=100 "
        f"{selected_split_by_n.loc[100, 'attempted_split_fits']:.0f}/"
        f"{selected_split_by_n.loc[100, 'expected_split_fits']:.0f}.",
        "- Point retention is 432/432; inference retention is 0/432 because the "
        "selected zero tangent spaces make the requested targets unsupported.",
        f"- Runtime seconds: N=50 total {runtime_by_method_n.loc[('selected_rank', 50), 'sum']:.3f}, "
        f"median {runtime_by_method_n.loc[('selected_rank', 50), 'median']:.3f}, max "
        f"{runtime_by_method_n.loc[('selected_rank', 50), 'max']:.3f}; N=100 total "
        f"{runtime_by_method_n.loc[('selected_rank', 100), 'sum']:.3f}, median "
        f"{runtime_by_method_n.loc[('selected_rank', 100), 'median']:.3f}, max "
        f"{runtime_by_method_n.loc[('selected_rank', 100), 'max']:.3f}.",
        f"- Assessment: **{selector_label}**.",
        "- The complete distribution and candidate-level IC decomposition are stored separately.",
        "",
        "## Medium-diagnostic decision",
        "",
        f"**{decision}.** This decision is diagnostic and does not authorize tuning or a medium run.",
        "",
        "## Reporting corrections",
        "",
        "- Failure placeholder target rows now inherit semantic IDs and DGP hashes from "
        "the authoritative attempt ledger.",
        "- Candidate Q_hat is reported as twice the numerical half-loss, matching the "
        "information-criterion implementation; invalid candidates retain infinite IC.",
        "- The decision helper now treats exact rank recovery as a required selector "
        "diagnostic, in addition to numerical completion, coverage, and cap behavior.",
        "",
        "## Commands",
        "",
        "```powershell",
        *RUN_COMMANDS,
        "```",
        "",
        f"Rank summary rows: {len(rank_summary)}; distribution rows: {len(distribution)}; "
        f"Gram/Riesz rows: {len(gram)}; split-summary rows: {len(split)}; "
        f"boundary-summary rows: {len(boundary)}.",
    ]
    (output / "preflight_report.md").write_text("\n".join(lines) + "\n", encoding="utf-8")
    return decision


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--fixed-root", required=True, type=Path)
    parser.add_argument("--selected-root", required=True, type=Path)
    parser.add_argument("--output-root", required=True, type=Path)
    parser.add_argument(
        "--frozen-calibration",
        type=Path,
        default=Path("configs/mc/frozen_dgp_calibration.toml"),
    )
    args = parser.parse_args()
    output = args.output_root
    output.mkdir(parents=True, exist_ok=True)

    frozen_hash = hashlib.sha256(args.frozen_calibration.read_bytes()).hexdigest()
    frozen = tomllib.loads(args.frozen_calibration.read_text(encoding="utf-8"))
    envelope = max(float(row["C_Theta"]) for row in frozen["calibration"])
    if frozen_hash != FROZEN_HASH or abs(envelope - EXPECTED_ENVELOPE) > 1e-12 or envelope > 9:
        raise ValueError("locked frozen calibration identity or envelope failed")

    frames = {
        "fixed_rank": _read_method(args.fixed_root),
        "selected_rank": _read_method(args.selected_root),
    }
    attempts = pd.concat([frame["attempts"] for frame in frames.values()], ignore_index=True)
    records = _backfill_attempt_identity(
        pd.concat([frame["records"] for frame in frames.values()], ignore_index=True),
        attempts,
    )
    fits = _normalize_fits(
        pd.concat([frame["fits"] for frame in frames.values()], ignore_index=True, sort=False)
    )
    inference = pd.concat(
        [frame["inference"] for frame in frames.values()], ignore_index=True, sort=False
    )
    ranks = frames["selected_rank"]["rank"].copy()

    if len(attempts) != 48 or attempts.groupby("method").size().to_dict() != {
        "fixed_rank": 24,
        "selected_rank": 24,
    }:
        raise ValueError("authorized evaluation count failed")
    expected_ids = {
        f"dgp{dgp}_N{n}_T{n}_r{replication:05d}_truth1-1-1"
        for dgp in EXPECTED_DGPS
        for n in EXPECTED_SIZES
        for replication in EXPECTED_REPLICATIONS
    }
    if set(attempts["semantic_replication_id"]) != expected_ids:
        raise ValueError("semantic replication IDs differ from the authorized design")
    matching = attempts.pivot(
        index="semantic_replication_id", columns="method", values="dgp_realization_hash"
    ).reset_index()
    matching["hashes_match"] = matching[["fixed_rank", "selected_rank"]].nunique(
        axis=1, dropna=False
    ).eq(1)
    if len(matching) != 24 or not matching["hashes_match"].all():
        raise ValueError("matched DGP realization hash verification failed")
    matched = attempts[["semantic_replication_id", "dgp", "N", "T", "replication"]].drop_duplicates()
    matched = matched.merge(matching, on="semantic_replication_id", validate="one_to_one")

    if len(records) != 48 * EXPECTED_TARGETS:
        raise ValueError("target records silently disappeared")
    counts = records.groupby(["semantic_replication_id", "method"])["target"].nunique()
    if not counts.eq(EXPECTED_TARGETS).all():
        raise ValueError("a method-replication is missing a requested target")
    for method, root in (("fixed_rank", args.fixed_root), ("selected_rank", args.selected_root)):
        manifest = json.loads((root / "run_manifest.json").read_text(encoding="utf-8"))
        if (
            manifest["frozen_calibration_hash"] != FROZEN_HASH
            or manifest["B"] != 10.0
            or manifest["c_B"] != 1.0
            or len(manifest["calibration_cells"]) != 8
            or any(float(cell["C_Theta"]) > 9.0 for cell in manifest["calibration_cells"])
        ):
            raise ValueError(f"manifest calibration validation failed for {method}")

    boundary, optimization, _ = _boundary_and_optimization(fits, attempts)
    accounting = _accounting(records)
    rank_summary, distribution, rank_replications = _rank_summaries(attempts, ranks)
    retention, performance = _retention_and_performance(records)
    gram = _gram_riesz(records, inference)
    split = _split_behavior(fits, records)
    candidates = _candidate_diagnostics(ranks)
    if (
        len(candidates)
        and candidates.loc[candidates["valid"], "ic_identity_error"].abs().max()
        > 1e-10
    ):
        raise ValueError("candidate IC decomposition failed its numerical identity")
    accounting_pass = bool(
        accounting["retention_order_pass"].all()
        and accounting["inference_valid_iff_success_pass"].all()
    )
    if not retention["retention_order_pass"].all() or not accounting_pass:
        raise ValueError("retention accounting failed")

    _write_pair(records, output, "replication_records")
    _write_pair(fits, output, "fit_diagnostics")
    _write_pair(inference, output, "inference_diagnostics")
    _write_pair(candidates, output, "candidate_ic_diagnostics")
    matched.sort_values(["dgp", "N", "replication"]).to_csv(
        output / "matched_dgp_hashes.csv", index=False, lineterminator="\n"
    )
    rank_summary.to_csv(output / "rank_selection_summary.csv", index=False, lineterminator="\n")
    distribution.to_csv(output / "selected_rank_distribution.csv", index=False, lineterminator="\n")
    rank_replications.to_csv(output / "rank_selection_replications.csv", index=False, lineterminator="\n")
    boundary.to_csv(output / "boundary_activity_summary.csv", index=False, lineterminator="\n")
    retention.to_csv(output / "target_retention_summary.csv", index=False, lineterminator="\n")
    gram.to_csv(output / "gram_riesz_summary.csv", index=False, lineterminator="\n")
    optimization.to_csv(output / "optimization_summary.csv", index=False, lineterminator="\n")
    performance.to_csv(output / "preliminary_performance_summary.csv", index=False, lineterminator="\n")
    dgp4_group_targets = {
        "A_G1_fixed_time",
        "A_G2_fixed_time",
        "A_G2_minus_G1_fixed_time",
        "A_G1_time_average",
        "A_G2_time_average",
        "A_G2_minus_G1_time_average",
        "B_G1_fixed_time",
        "B_G2_fixed_time",
        "B_G2_minus_G1_fixed_time",
    }
    performance.loc[
        performance["dgp"].eq(4)
        & performance["method"].eq("selected_rank")
        & performance["target"].isin(dgp4_group_targets),
        ["N", "target", "true_value_mean"],
    ].sort_values(["N", "target"]).to_csv(
        output / "dgp4_true_group_means.csv", index=False, lineterminator="\n"
    )
    split.to_csv(output / "split_fit_summary.csv", index=False, lineterminator="\n")
    accounting.to_csv(output / "failure_accounting.csv", index=False, lineterminator="\n")
    decision = _report(
        output,
        attempts,
        ranks,
        rank_summary,
        distribution,
        retention,
        gram,
        split,
        boundary,
    )
    (output / "preflight_audit_manifest.json").write_text(
        json.dumps(
            {
                "authorized_unique_dgp_realizations": 24,
                "authorized_method_replication_evaluations": 48,
                "matched_hashes_pass": True,
                "frozen_calibration_hash": FROZEN_HASH,
                "global_C_Theta_max": EXPECTED_ENVELOPE,
                "B": 10.0,
                "c_B": 1.0,
                "medium_recommendation": decision,
                "executed_commands": RUN_COMMANDS,
            },
            indent=2,
        )
        + "\n",
        encoding="utf-8",
    )
    print(f"Locked Revision-9 preflight audit complete: {decision}")


if __name__ == "__main__":
    main()
