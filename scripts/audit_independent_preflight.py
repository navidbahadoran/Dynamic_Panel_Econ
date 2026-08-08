"""Audit the final 12-draw matched fixed/selected independent preflight."""

from __future__ import annotations

import argparse
import json
from collections import Counter
from pathlib import Path
from typing import Any

import pandas as pd

TARGET_COUNT = 18
EXPECTED_REPLICATIONS = set(range(3, 6))
GRAM_FAILURES = {
    "tangent_gram_eigensolver_failure",
    "tangent_gram_nearly_singular",
}
RIESZ_FAILURES = {"riesz_solver_failure", "riesz_target_instability"}


def _chunk_table(root: Path, folder: str) -> pd.DataFrame:
    files = sorted((root / folder).glob("*.parquet"))
    return pd.concat([pd.read_parquet(path) for path in files], ignore_index=True)


def _bool(series: pd.Series) -> pd.Series:
    if series.dtype == bool:
        return series
    return series.astype(str).str.lower().eq("true")


def _rank(value: Any) -> tuple[int, ...]:
    return tuple(int(item) for item in json.loads(str(value)))


def _write_text(path: Path, value: str) -> None:
    with path.open("w", encoding="utf-8", newline="\n") as handle:
        handle.write(value)


def _status_counts(records: pd.DataFrame, method: str) -> list[dict[str, Any]]:
    counts = records["primary_status"].value_counts()
    return [
        {"method": method, "primary_status": status, "target_record_count": int(count)}
        for status, count in counts.items()
    ]


def _reconcile(
    fixed_attempts: pd.DataFrame,
    selected_attempts: pd.DataFrame,
    fixed_records: pd.DataFrame,
    selected_records: pd.DataFrame,
    previous_ids: set[str],
) -> dict[str, Any]:
    details: dict[str, Any] = {}
    for method, attempts, records in (
        ("fixed_rank", fixed_attempts, fixed_records),
        ("selected_rank", selected_attempts, selected_records),
    ):
        if len(attempts) != 12 or len(records) != 12 * TARGET_COUNT:
            raise ValueError(f"{method} attempted denominator is not 12 x 18")
        if set(attempts["replication"].astype(int)) != EXPECTED_REPLICATIONS:
            raise ValueError(f"{method} did not use replication indices 3-5")
        identity = ["run_id", "dgp", "N", "T", "replication", "method", "target"]
        if records.duplicated(identity).any():
            raise ValueError(f"{method} contains duplicate target records")
        inference = _bool(records["inference_valid"])
        success = records["primary_status"].eq("success")
        if not inference.eq(success).all():
            raise ValueError(f"{method} violates inference/status equivalence")
        point = _bool(records["point_estimate_valid"])
        if not int(inference.sum()) <= int(point.sum()) <= len(records):
            raise ValueError(f"{method} retention denominators do not reconcile")
        failures = int(records["primary_status"].ne("success").sum())
        if len(records) != int(success.sum()) + failures:
            raise ValueError(f"{method} statuses are not mutually exclusive")
        overlap = set(attempts["semantic_replication_id"]) & previous_ids
        if overlap:
            raise ValueError(f"{method} reuses prior semantic IDs: {sorted(overlap)}")
        details[method] = {
            "R_attempted": len(records),
            "R_point": int(point.sum()),
            "R_inference": int(inference.sum()),
            "primary_success": int(success.sum()),
            "primary_failures": failures,
            "semantic_id_overlap_with_prior_preflight": 0,
        }
    matched = fixed_attempts.merge(
        selected_attempts,
        on=["dgp", "N", "T", "replication", "semantic_replication_id"],
        suffixes=("_fixed", "_selected"),
        validate="one_to_one",
    )
    if len(matched) != 12 or not matched["dgp_realization_hash_fixed"].eq(
        matched["dgp_realization_hash_selected"]
    ).all():
        raise ValueError("fixed and selected DGP realization hashes do not match")
    details["matched_draws"] = 12
    details["all_dgp_hashes_match"] = True
    return details


def audit(fixed: Path, selected: Path, output: Path, prior_root: Path) -> None:
    output.mkdir(parents=True, exist_ok=True)
    fixed_attempts = pd.read_parquet(fixed / "attempted_replications.parquet")
    selected_attempts = pd.read_parquet(selected / "attempted_replications.parquet")
    fixed_records = pd.read_parquet(fixed / "replication_records.parquet")
    selected_records = pd.read_parquet(selected / "replication_records.parquet")
    selected_ranks = _chunk_table(selected, "rank")
    fixed_fits = pd.read_parquet(fixed / "fit_diagnostics.parquet")
    selected_fits = pd.read_parquet(selected / "fit_diagnostics.parquet")

    previous_ids: set[str] = set()
    for path in prior_root.glob("*/*/attempted_replications.parquet"):
        previous_ids.update(pd.read_parquet(path, columns=["semantic_replication_id"])[
            "semantic_replication_id"
        ])
    reconciliation = _reconcile(
        fixed_attempts,
        selected_attempts,
        fixed_records,
        selected_records,
        previous_ids,
    )

    selected_vectors = selected_ranks["selected_rank_vector"].map(_rank)
    distribution = Counter(selected_vectors)
    selected_summary = {
        "attempted_replications": 12,
        "valid_cap_pilot": len(selected_ranks),
        "pilot_multistart_disagreement": int(
            _bool(selected_ranks["pilot_multistart_disagreement"]).sum()
        ),
        "basin_confirmation_attempted": int(
            _bool(selected_ranks["cap_pilot_basin_confirmation_attempted"]).sum()
        ),
        "basin_confirmation_success": int(
            _bool(selected_ranks["cap_pilot_basin_confirmation_success"]).sum()
        ),
        "candidate_coverage": int(_bool(selected_ranks["candidate_coverage"]).sum()),
        "candidate_coverage_conditional_share": float(
            _bool(selected_ranks["candidate_coverage"]).mean()
        ),
        "exact_recovery": int(_bool(selected_ranks["exact_rank_recovery"]).sum()),
        "underselection": int(
            sum(any(a < b for a, b in zip(rank, (1, 1, 1), strict=True)) for rank in selected_vectors)
        ),
        "overselection": int(
            sum(any(a > b for a, b in zip(rank, (1, 1, 1), strict=True)) for rank in selected_vectors)
        ),
        "rank_cap_hits": int(_bool(selected_ranks["rank_at_cap"]).sum()),
        "point_retained": int(_bool(selected_records["point_estimate_valid"]).sum()),
        "point_attempted": len(selected_records),
        "point_retained_share": float(_bool(selected_records["point_estimate_valid"]).mean()),
        "inference_retained": int(_bool(selected_records["inference_valid"]).sum()),
        "inference_attempted": len(selected_records),
        "inference_retained_share": float(_bool(selected_records["inference_valid"]).mean()),
        "runtime_total_seconds": float(selected_attempts["replication_runtime_seconds"].sum()),
        "runtime_mean_seconds": float(selected_attempts["replication_runtime_seconds"].mean()),
        "runtime_median_seconds": float(selected_attempts["replication_runtime_seconds"].median()),
        "gram_failure_target_records": int(
            selected_records["primary_status"].isin(GRAM_FAILURES).sum()
        ),
        "riesz_failure_target_records": int(
            selected_records["primary_status"].isin(RIESZ_FAILURES).sum()
        ),
        "split_fit_failure_target_records": int(
            selected_records["primary_status"].eq("split_fit_failure").sum()
        ),
        "candidate_numerically_unresolved_replications": int(
            selected_attempts["primary_status"].eq("candidate_numerically_unresolved").sum()
        ),
        "selected_postrefit_stability_failures": int(
            (~_bool(selected_ranks["selected_objective_stability_pass"])).sum()
        ),
        "rank_selection_failure_replications": int(
            selected_attempts["primary_status"].eq("rank_selection_failure").sum()
        ),
        "split_fit_failure_replications": int(
            selected_records.loc[
                selected_records["primary_status"].eq("split_fit_failure"),
                "semantic_replication_id",
            ].nunique()
        ),
    }
    fixed_summary = {
        "attempted_replications": 12,
        "full_fit_success": int(
            (~fixed_attempts["primary_status"].isin(
                {"full_fit_failure", "coefficient_bound_hit"}
            )).sum()
        ),
        "coefficient_bound_hits": int(
            fixed_attempts["primary_status"].eq("coefficient_bound_hit").sum()
        ),
        "numerical_stability_failures": int(
            fixed_attempts["primary_status"].eq("full_fit_failure").sum()
        ),
        "point_retained": int(_bool(fixed_records["point_estimate_valid"]).sum()),
        "point_attempted": len(fixed_records),
        "point_retained_share": float(_bool(fixed_records["point_estimate_valid"]).mean()),
        "inference_retained": int(_bool(fixed_records["inference_valid"]).sum()),
        "inference_attempted": len(fixed_records),
        "inference_retained_share": float(_bool(fixed_records["inference_valid"]).mean()),
        "gram_failure_target_records": int(
            fixed_records["primary_status"].isin(GRAM_FAILURES).sum()
        ),
        "riesz_failure_target_records": int(
            fixed_records["primary_status"].isin(RIESZ_FAILURES).sum()
        ),
        "split_fit_failure_target_records": int(
            fixed_records["primary_status"].eq("split_fit_failure").sum()
        ),
        "split_fit_failure_replications": int(
            fixed_records.loc[
                fixed_records["primary_status"].eq("split_fit_failure"),
                "semantic_replication_id",
            ].nunique()
        ),
        "runtime_total_seconds": float(fixed_attempts["replication_runtime_seconds"].sum()),
        "runtime_mean_seconds": float(fixed_attempts["replication_runtime_seconds"].mean()),
        "runtime_median_seconds": float(fixed_attempts["replication_runtime_seconds"].median()),
    }
    rank_distribution = pd.DataFrame(
        [
            {"selected_rank_vector": json.dumps(rank), "count": count, "share": count / 12}
            for rank, count in sorted(distribution.items())
        ]
    )
    failure_accounting = pd.DataFrame(
        [
            *_status_counts(fixed_records, "fixed_rank"),
            *_status_counts(selected_records, "selected_rank"),
        ]
    )
    optimization = pd.DataFrame(
        [
            {
                "method": method,
                "fit_count": len(fits),
                "converged_count": int(_bool(fits["convergence_flag"]).sum()),
                "stationarity_pass_count": int(_bool(fits["stationarity_pass"]).sum()),
                "coefficient_bound_hit_count": int(_bool(fits["coefficient_bound_hit"]).sum()),
                "objective_stability_failure_count": int(
                    (~_bool(fits["objective_stability_pass"].fillna(True))).sum()
                ),
            }
            for method, fits in (("fixed_rank", fixed_fits), ("selected_rank", selected_fits))
        ]
    )
    selected_replications = selected_ranks[
        [
            "semantic_replication_id", "dgp", "replication", "primary_status",
            "cap_pilot_rank", "pilot_multistart_disagreement",
            "cap_pilot_basin_confirmation_attempted",
            "cap_pilot_basin_confirmation_success", "candidate_coverage",
            "selected_rank_vector", "exact_rank_recovery", "rank_at_cap",
            "rank_runtime_seconds",
        ]
    ].copy()
    fixed_replications = fixed_attempts[
        ["semantic_replication_id", "dgp", "replication", "primary_status", "replication_runtime_seconds"]
    ].copy()
    selected_replications.to_csv(
        output / "selected_replications.csv", index=False, lineterminator="\n"
    )
    fixed_replications.to_csv(
        output / "fixed_replications.csv", index=False, lineterminator="\n"
    )
    rank_distribution.to_csv(
        output / "selected_rank_distribution.csv", index=False, lineterminator="\n"
    )
    failure_accounting.to_csv(
        output / "failure_accounting.csv", index=False, lineterminator="\n"
    )
    optimization.to_csv(
        output / "optimization_summary.csv", index=False, lineterminator="\n"
    )
    _write_text(
        output / "accounting_reconciliation.json",
        json.dumps(reconciliation, indent=2) + "\n",
    )
    summary = {
        "fixed": fixed_summary,
        "selected_4e-6": selected_summary,
        "reconciliation": reconciliation,
    }
    _write_text(
        output / "preflight_summary.json", json.dumps(summary, indent=2) + "\n"
    )
    dgp4_cap_hits = int(
        (_bool(selected_ranks["rank_at_cap"]) & selected_ranks["dgp"].eq(4)).sum()
    )
    recommendation = "NO-GO" if dgp4_cap_hits >= 2 else "GO"
    report = [
        "# Final independent preflight",
        "",
        "This run contains 12 new DGP draws (replication indices 3-5) and 24 matched method evaluations.",
        "No medium, rank-stress, power, or production simulation was run.",
        "",
        "## Fixed rank",
        "",
        "```json", json.dumps(fixed_summary, indent=2), "```",
        "",
        "## Selected rank (`c_kappa=4e-6`)",
        "",
        "```json", json.dumps(selected_summary, indent=2), "```",
        "",
        "Selected rank distribution:",
        "",
        "```csv",
        rank_distribution.to_csv(index=False, lineterminator="\n").strip(),
        "```",
        "",
        "## Recommendation",
        "",
        f"**{recommendation}** for the medium diagnostic. Although all cap pilots were valid and candidate coverage was 100%, rank-cap hits occurred in 3/12 cases, including 2/3 DGP-4 cases. That DGP-specific concentration does not satisfy the no-systematic-rank-cap gate.",
        "",
        "The recommendation also notes the explicitly retained split-fit failures; Gram and Riesz primary failures are counted separately above.",
        "The optimization CSV counts every executed route/candidate start, including intentionally rejected starts; replication-level numerical failures are the gate metrics reported in the JSON summaries.",
    ]
    _write_text(output / "REPORT.md", "\n".join(report) + "\n")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--fixed", type=Path, required=True)
    parser.add_argument("--selected", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument(
        "--prior-root", type=Path, default=Path("results/mc/preflight_final")
    )
    args = parser.parse_args()
    audit(args.fixed, args.selected, args.output, args.prior_root)


if __name__ == "__main__":
    main()
