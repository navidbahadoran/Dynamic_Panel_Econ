"""Offline diagnosis of the locked Revision-9 statistical-preflight no-go.

This script reads persisted preflight records only.  It does not generate a DGP,
fit a coefficient matrix, solve a Riesz system, or alter any locked procedure.
"""

from __future__ import annotations

import argparse
import json
import math
from collections import Counter
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd

from dynamic_panel_econ.ic_break_even import (
    ICPoint,
    optimality_interval,
    pairwise_break_even,
    rank_increment_dimension,
    select_ic_candidate,
)
from dynamic_panel_econ.rank_selection import model_dimension, revision8_kappa

TRUE_RANK = (1, 1, 1)
ZERO_RANK = (0, 0, 0)
ETA = 4.0
SPATIAL_DIMENSION = 1
BOUND = 10.0
STATIONARITY_TOL = 1e-6
BROAD_TARGETS = {
    "A_full_mean",
    "B_full_mean",
    "A_G1_time_average",
    "A_G2_time_average",
    "A_G2_minus_G1_time_average",
    "B_G1_time_average",
    "B_G2_time_average",
    "B_G2_minus_G1_time_average",
}
FAILED_FIXED_IDS = {
    "dgp1_N50_T50_r2026080902_truth1-1-1",
    "dgp3_N50_T50_r2026080902_truth1-1-1",
    "dgp3_N50_T50_r2026080903_truth1-1-1",
    "dgp4_N50_T50_r2026080901_truth1-1-1",
}


def _rank(value: Any) -> tuple[int, int, int]:
    parsed = json.loads(value) if isinstance(value, str) else value
    return tuple(int(item) for item in parsed)  # type: ignore[return-value]


def _rank_text(value: tuple[int, ...]) -> str:
    return "(" + ",".join(str(item) for item in value) + ")"


def _bool(value: Any) -> bool:
    if isinstance(value, str):
        return value.strip().lower() == "true"
    return bool(value)


def _read_chunks(root: Path) -> pd.DataFrame:
    files = sorted(root.glob("*.parquet"))
    if not files:
        raise FileNotFoundError(f"no parquet chunks under {root}")
    return pd.concat((pd.read_parquet(path) for path in files), ignore_index=True)


def _candidate_reason_map(rank_rows: pd.DataFrame) -> dict[tuple[str, tuple[int, ...]], list[str]]:
    result: dict[tuple[str, tuple[int, ...]], list[str]] = {}
    for saved in rank_rows.itertuples():
        diagnostics = json.loads(saved.rank_diagnostics_json)
        for record in diagnostics["candidate_records"]:
            result[(saved.semantic_replication_id, tuple(record["ranks"]))] = list(
                record.get("invalid_reasons", [])
            )
    return result


def _candidate_points(group: pd.DataFrame) -> list[ICPoint]:
    return [
        ICPoint(_rank(row.rank_vector), float(row.Q_hat), int(row.d_r))
        for row in group.itertuples()
        if _bool(row.valid)
    ]


def _penalty_scale() -> pd.DataFrame:
    rows: list[dict[str, Any]] = []
    for n in (50, 100, 200, 400):
        nt = n * n
        log_nt = math.log(nt)
        b_nt = nt ** (1.0 / (8.0 + ETA)) * log_nt
        kappa = revision8_kappa(
            n,
            n,
            eta_for_penalty=ETA,
            spatial_dimension=SPATIAL_DIMENSION,
        )
        increments = {
            f"penalty_increment_rank_{rank}_to_{rank + 1}": (
                kappa * rank_increment_dimension(rank, n, n) / nt
            )
            for rank in (0, 1, 2)
        }
        rows.append(
            {
                "N": n,
                "T": n,
                "NT": nt,
                "eta": ETA,
                "d_s": SPATIAL_DIMENSION,
                "log_NT": log_nt,
                "b_NT": b_nt,
                "kappa_base": kappa,
                **increments,
                "dimension_truth": model_dimension(TRUE_RANK, n, n),
                "dimension_zero": model_dimension(ZERO_RANK, n, n),
                "total_penalty_truth_minus_zero_c_1": (
                    kappa
                    * (
                        model_dimension(TRUE_RANK, n, n)
                        - model_dimension(ZERO_RANK, n, n)
                    )
                    / nt
                ),
            }
        )
    return pd.DataFrame(rows)


def _ic_outputs(
    candidates: pd.DataFrame,
    rank_rows: pd.DataFrame,
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame, float]:
    reason_map = _candidate_reason_map(rank_rows)
    validity_rows: list[dict[str, Any]] = []
    break_rows: list[dict[str, Any]] = []
    interval_rows: list[dict[str, Any]] = []
    max_identity_error = 0.0

    for semantic_id, group in candidates.groupby("semantic_replication_id", sort=True):
        group = group.copy()
        first = group.iloc[0]
        n, t = int(first["N"]), int(first["T"])
        kappa = revision8_kappa(n, t)
        base_per_dimension = kappa / (n * t)
        reconstructed = (
            group["log_Q_hat"].astype(float)
            + base_per_dimension * group["d_r"].astype(float)
        )
        valid_mask = group["valid"].map(_bool)
        errors = np.abs(reconstructed[valid_mask] - group.loc[valid_mask, "final_IC"].astype(float))
        if len(errors):
            max_identity_error = max(max_identity_error, float(errors.max()))

        true_rows = group.loc[
            group["rank_vector"].map(lambda value: _rank(value) == TRUE_RANK)
        ]
        if len(true_rows) != 1:
            raise AssertionError(f"expected one true-rank row for {semantic_id}")
        true_row = true_rows.iloc[0]
        truth_valid = _bool(true_row.valid)
        reasons = reason_map.get((semantic_id, TRUE_RANK), [])
        validity_rows.append(
            {
                "semantic_replication_id": semantic_id,
                "dgp": int(first.dgp),
                "N": n,
                "T": t,
                "replication": int(first.replication),
                "true_rank": _rank_text(TRUE_RANK),
                "true_rank_present_in_candidate_set": True,
                "true_rank_postrefit_numerically_valid": truth_valid,
                "true_rank_postrefit_numerically_invalid": not truth_valid,
                "invalid_reasons": "|".join(reasons),
                "Q_hat_true": float(true_row.Q_hat),
                "log_Q_hat_true": float(true_row.log_Q_hat),
                "d_true": int(true_row.d_r),
                "baseline_IC_true": float(true_row.final_IC),
            }
        )

        points = _candidate_points(group)
        by_rank = {point.ranks: point for point in points}
        if not truth_valid:
            interval_rows.append(
                {
                    "semantic_replication_id": semantic_id,
                    "dgp": int(first.dgp),
                    "N": n,
                    "T": t,
                    "replication": int(first.replication),
                    "eligible": False,
                    "c_lower": np.nan,
                    "c_upper": np.nan,
                    "lower_binding_candidate": "",
                    "upper_binding_candidate": "",
                    "nonempty": False,
                    "ineligibility_reason": "true_rank_postrefit_numerically_invalid",
                }
            )
            continue

        truth = by_rank[TRUE_RANK]
        interval = optimality_interval(
            truth,
            points,
            base_penalty_per_dimension=base_per_dimension,
        )
        interval_rows.append(
            {
                "semantic_replication_id": semantic_id,
                "dgp": int(first.dgp),
                "N": n,
                "T": t,
                "replication": int(first.replication),
                "eligible": True,
                "c_lower": interval.lower,
                "c_upper": interval.upper,
                "lower_binding_candidate": ""
                if interval.lower_binding_rank is None
                else _rank_text(interval.lower_binding_rank),
                "upper_binding_candidate": ""
                if interval.upper_binding_rank is None
                else _rank_text(interval.upper_binding_rank),
                "nonempty": not interval.empty and interval.upper >= max(0.0, interval.lower),
                "ineligibility_reason": "",
            }
        )
        for competitor in points:
            if competitor.ranks == TRUE_RANK:
                continue
            c_star = pairwise_break_even(
                truth,
                competitor,
                base_penalty_per_dimension=base_per_dimension,
            )
            break_rows.append(
                {
                    "semantic_replication_id": semantic_id,
                    "dgp": int(first.dgp),
                    "N": n,
                    "T": t,
                    "replication": int(first.replication),
                    "reference_rank": _rank_text(TRUE_RANK),
                    "competitor_rank": _rank_text(competitor.ranks),
                    "Q_hat_true": truth.qhat,
                    "Q_hat_competitor": competitor.qhat,
                    "log_Q_hat_true": truth.log_qhat,
                    "log_Q_hat_competitor": competitor.log_qhat,
                    "d_true": truth.dimension,
                    "d_competitor": competitor.dimension,
                    "P_true": base_per_dimension * truth.dimension,
                    "P_competitor": base_per_dimension * competitor.dimension,
                    "c_break_even": c_star,
                    "positive_finite_break_even": np.isfinite(c_star) and c_star > 0.0,
                    "especially_requested_comparison": competitor.ranks
                    in {
                        ZERO_RANK,
                        (1, 0, 0),
                        (0, 1, 0),
                        (0, 0, 1),
                        (1, 1, 0),
                        (1, 0, 1),
                        (0, 1, 1),
                    },
                }
            )

    if max_identity_error > 1e-10:
        raise AssertionError(f"locked IC reconstruction error {max_identity_error}")
    return (
        pd.DataFrame(validity_rows),
        pd.DataFrame(break_rows),
        pd.DataFrame(interval_rows),
        max_identity_error,
    )


def _intersection_summary(intervals: pd.DataFrame, label: str) -> dict[str, Any]:
    eligible = intervals.loc[intervals["eligible"].map(_bool)].copy()
    nonempty = eligible.loc[eligible["nonempty"].map(_bool)].copy()
    if len(nonempty) != len(eligible) or eligible.empty:
        common_lower, common_upper, common_nonempty = np.nan, np.nan, False
    else:
        common_lower = max(0.0, float(nonempty["c_lower"].max()))
        common_upper = float(nonempty["c_upper"].min())
        common_nonempty = common_upper >= common_lower and common_upper > 0.0

    candidates: set[float] = set()
    for row in nonempty.itertuples():
        lower = max(0.0, float(row.c_lower))
        upper = float(row.c_upper)
        if upper <= 0.0:
            continue
        candidates.add(max(lower, np.finfo(float).tiny))
        if np.isfinite(upper):
            candidates.add(upper)
            if lower > 0.0:
                candidates.add(math.sqrt(lower * upper))
            else:
                candidates.add(upper / 2.0)
    maximum_overlap = 0
    witness = np.nan
    for value in sorted(candidates):
        overlap = sum(
            max(0.0, float(row.c_lower)) <= value <= float(row.c_upper)
            for row in nonempty.itertuples()
        )
        if overlap > maximum_overlap:
            maximum_overlap, witness = overlap, value
    return {
        "scope": label,
        "eligible_replications": len(eligible),
        "nonempty_individual_intervals": len(nonempty),
        "common_c_lower": common_lower,
        "common_c_upper": common_upper,
        "common_intersection_nonempty": common_nonempty,
        "maximum_simultaneous_overlap": maximum_overlap,
        "maximum_overlap_c_witness": witness,
    }


def _grid_values(breaks: pd.DataFrame, intervals: pd.DataFrame) -> list[float]:
    observed = set(
        float(value)
        for value in breaks.loc[breaks["positive_finite_break_even"], "c_break_even"]
        if np.isfinite(value) and value > 0.0
    )
    for column in ("c_lower", "c_upper"):
        observed.update(
            float(value)
            for value in intervals[column].dropna()
            if np.isfinite(value) and value > 0.0
        )
    if not observed:
        return [1.0]
    ordered = sorted(observed)
    grid = {ordered[0] / 10.0, ordered[-1] * 10.0, 1.0, *ordered}
    grid.update(math.sqrt(left * right) for left, right in zip(ordered[:-1], ordered[1:], strict=True))
    return sorted(value for value in grid if np.isfinite(value) and value > 0.0)


def _selection_class(ranks: tuple[int, ...]) -> str:
    below = [rank < truth for rank, truth in zip(ranks, TRUE_RANK, strict=True)]
    above = [rank > truth for rank, truth in zip(ranks, TRUE_RANK, strict=True)]
    if ranks == TRUE_RANK:
        return "Exact"
    if any(below) and not any(above):
        return "Under only"
    if any(above) and not any(below):
        return "Over only"
    return "Mixed"


def _offline_grid(candidates: pd.DataFrame, grid_values: list[float]) -> pd.DataFrame:
    point_groups = {
        semantic_id: (
            int(group.iloc[0].N),
            _candidate_points(group),
        )
        for semantic_id, group in candidates.groupby("semantic_replication_id", sort=True)
    }
    rows: list[dict[str, Any]] = []
    for c_value in grid_values:
        for n in (50, 100):
            selections: list[tuple[int, ...]] = []
            classes: list[str] = []
            for _, (saved_n, points) in point_groups.items():
                if saved_n != n:
                    continue
                choice = select_ic_candidate(
                    points,
                    c_value,
                    base_penalty_per_dimension=revision8_kappa(n, n) / (n * n),
                )
                selections.append(choice.ranks)
                classes.append(_selection_class(choice.ranks))
            distribution = Counter(_rank_text(rank) for rank in selections)
            class_counts = Counter(classes)
            modal_count = max(distribution.values())
            modal_rank = sorted(rank for rank, count in distribution.items() if count == modal_count)[0]
            rows.append(
                {
                    "diagnostic_label": "OFFLINE EXISTING-CANDIDATE DIAGNOSTIC",
                    "c_kappa": c_value,
                    "N": n,
                    "eligible_saved_replications": len(selections),
                    "Exact": class_counts["Exact"],
                    "Under_only": class_counts["Under only"],
                    "Over_only": class_counts["Over only"],
                    "Mixed": class_counts["Mixed"],
                    "modal_selected_rank": modal_rank,
                    "modal_count": modal_count,
                    "selected_rank_distribution": json.dumps(dict(sorted(distribution.items()))),
                }
            )
    return pd.DataFrame(rows)


def _fixed_failure_details(fixed_rank_rows: pd.DataFrame, fits: pd.DataFrame) -> pd.DataFrame:
    rows: list[dict[str, Any]] = []
    failed = fixed_rank_rows.loc[
        fixed_rank_rows["semantic_replication_id"].isin(FAILED_FIXED_IDS)
    ]
    if set(failed["semantic_replication_id"]) != FAILED_FIXED_IDS:
        raise AssertionError("the four requested fixed-rank failures were not all found")
    fit_lookup = fits.loc[
        (fits["method"] == "fixed_rank")
        & (fits["fit_type"] == "full_fixed_rank")
        & fits["semantic_replication_id"].isin(FAILED_FIXED_IDS)
    ].set_index(["semantic_replication_id", "start_number"])
    for saved in failed.itertuples():
        diagnostics = json.loads(saved.fixed_rank_multistart_diagnostics)
        for start_number, record in enumerate(diagnostics["original_start_records"], 1):
            fit = fit_lookup.loc[(saved.semantic_replication_id, start_number)]
            invalid_reasons = list(record["invalid_reasons"])
            if "stationarity_high" in invalid_reasons:
                classification = "KKT/optimality failure"
            elif "numerical_rank_loss" in invalid_reasons:
                classification = "numerical rank loss"
            elif "constrained_feasibility_failure" in invalid_reasons:
                classification = "feasibility failure"
            elif not record["convergence"]:
                classification = "solver nonconvergence"
            elif not diagnostics["objective_stability_pass"]:
                classification = "start disagreement"
            else:
                classification = "other"
            rows.append(
                {
                    "semantic_replication_id": saved.semantic_replication_id,
                    "dgp": int(saved.dgp),
                    "N": int(saved.N),
                    "replication": int(saved.replication),
                    "start_number": start_number,
                    "start_id": record["start_id"],
                    "objective": float(record["final_objective"]),
                    "objective_gap_to_best": float(record["objective_gap_to_best"]),
                    "converged": bool(record["convergence"]),
                    "feasible": float(record["coefficient_envelope"]) <= BOUND,
                    "max_absolute_coefficient": float(record["coefficient_envelope"]),
                    "unconstrained_outside_box": _bool(fit.unconstrained_outside_box),
                    "constrained_fallback": _bool(fit.constrained_fallback_used),
                    "constrained_iterations": int(fit.constrained_iterations),
                    "stationarity_or_KKT_residual": float(record["stationarity"]),
                    "numerical_rank": _rank_text(tuple(record["numerical_rank"])),
                    "runtime_seconds": float(record["runtime"]),
                    "start_numerically_valid": bool(record["valid"]),
                    "exact_failure_reason": "|".join(invalid_reasons),
                    "replication_primary_classification": classification,
                    "objective_stability_pass": bool(diagnostics["objective_stability_pass"]),
                    "final_acceptance_basis": diagnostics["final_acceptance_basis"],
                }
            )
    return pd.DataFrame(rows).sort_values(
        ["dgp", "replication", "start_number"], kind="stable"
    )


def _split_boundary_details(fits: pd.DataFrame, fixed_raw: pd.DataFrame) -> pd.DataFrame:
    split = fits.loc[
        (fits["method"] == "fixed_rank")
        & (fits["N"] == 50)
        & fits["fit_type"].isin(
            ["time_split_0", "time_split_1", "unit_split_0", "unit_split_1"]
        )
    ].copy()
    boundary = split.loc[split["boundary_active"].map(_bool)].copy()
    if len(boundary) != 9:
        raise AssertionError(f"expected nine boundary-active split fits, found {len(boundary)}")
    rows: list[dict[str, Any]] = []
    for saved in boundary.itertuples():
        split_type, half_text = str(saved.fit_type).rsplit("_", 1)
        rows.append(
            {
                "record_scope": "boundary_fit",
                "semantic_replication_id": saved.semantic_replication_id,
                "dgp": int(saved.dgp),
                "N": int(saved.N),
                "replication": int(saved.replication),
                "split_type": split_type.replace("_split", ""),
                "half": int(half_text) + 1,
                "active_coefficient_matrix": "not_persisted_in_locked_records",
                "max_absolute_coefficient": float(saved.max_abs_coefficient),
                "number_of_active_entries": "not_persisted_in_locked_records",
                "unconstrained_max_absolute_coefficient": float(saved.unconstrained_max_abs),
                "unconstrained_overshoot_beyond_B": max(
                    0.0, float(saved.unconstrained_max_abs) - BOUND
                ),
                "constrained_objective": float(saved.constrained_objective),
                "KKT_residual": float(saved.constrained_KKT_residual),
                "constrained_iterations": int(saved.constrained_iterations),
                "constrained_solver_status": saved.constrained_solver_status,
                "boundary_fit_count_in_replication": np.nan,
                "boundary_count_category": "",
                "broad_target_status": "",
                "broad_targets_attempted": np.nan,
                "broad_targets_inference_valid": np.nan,
            }
        )

    broad = fixed_raw.loc[(fixed_raw["N"] == 50) & fixed_raw["target"].isin(BROAD_TARGETS)]
    for semantic_id, group in split.groupby("semantic_replication_id", sort=True):
        count = int(group["boundary_active"].map(_bool).sum())
        target_group = broad.loc[broad["semantic_replication_id"] == semantic_id]
        statuses = sorted(set(target_group["primary_status"].dropna().astype(str)))
        rows.append(
            {
                "record_scope": "replication_summary",
                "semantic_replication_id": semantic_id,
                "dgp": int(group.iloc[0].dgp),
                "N": 50,
                "replication": int(group.iloc[0].replication),
                "split_type": "all_four",
                "half": np.nan,
                "active_coefficient_matrix": "",
                "max_absolute_coefficient": np.nan,
                "number_of_active_entries": "",
                "unconstrained_max_absolute_coefficient": np.nan,
                "unconstrained_overshoot_beyond_B": np.nan,
                "constrained_objective": np.nan,
                "KKT_residual": np.nan,
                "constrained_iterations": np.nan,
                "constrained_solver_status": "",
                "boundary_fit_count_in_replication": count,
                "boundary_count_category": "0" if count == 0 else "1" if count == 1 else "multiple",
                "broad_target_status": "|".join(statuses),
                "broad_targets_attempted": len(target_group),
                "broad_targets_inference_valid": int(
                    target_group["primary_status"].eq("success").sum()
                ),
            }
        )
    return pd.DataFrame(rows)


def _theorem_text() -> str:
    return """# Theorem implication of a fixed positive multiplier

The maintained rate definitions are

```text
zeta_NT = b_NT^2 log(NT)^(d_s+2) (1/N + 1/T)
kappa_NT = c_kappa b_NT^2 log(NT)^(d_s+3).
```

Therefore

```text
zeta_NT / {kappa_NT (N+T)/(NT)} = 1 / {c_kappa log(NT)} -> 0
```

for every fixed `c_kappa > 0`. Multiplication by a fixed positive constant also does not alter
the maintained requirement `kappa_NT (N+T)/(NT) -> 0`. With `eta=4`, `d_s=1`, and balanced
`N=T=m`, that term is proportional to `c_kappa m^(-2/3) log(m^2)^6`, which converges to zero.

Thus `c_kappa=1` is a normalization, not a constant uniquely implied by the proof, and any fixed
positive multiplier preserves these asymptotic rate statements. That does not authorize choosing
a multiplier after inspecting production outcomes. If the paper is revised, the multiplier and
its rationale must be fixed before a fresh, independent validation experiment; the revised paper,
configuration, seeds, decision rule, and permitted sensitivities must all be locked in advance.
"""


def _common_markdown(summaries: list[dict[str, Any]], intervals: pd.DataFrame) -> str:
    lines = [
        "# Common fixed-c intersection",
        "",
        "Only replications with a numerically valid true-rank post-refit are eligible.",
        "",
        "| Scope | Eligible | Individual nonempty | Common lower | Common upper | Common? | Maximum overlap | Witness c |",
        "|:--|--:|--:|--:|--:|:--:|--:|--:|",
    ]
    for row in summaries:
        lines.append(
            "| {scope} | {eligible_replications} | {nonempty_individual_intervals} | "
            "{common_c_lower:.8g} | {common_c_upper:.8g} | {common_intersection_nonempty} | "
            "{maximum_simultaneous_overlap} | {maximum_overlap_c_witness:.8g} |".format(**row)
        )
    lines += [
        "",
        "The interval endpoints and binding saved candidates are in `truth_optimal_intervals.csv`.",
        "No multiplier is activated or recommended by this calculation.",
        "",
        f"Eligible N=50 intervals: {int(((intervals.N == 50) & intervals.eligible).sum())}/12.",
        f"Eligible N=100 intervals: {int(((intervals.N == 100) & intervals.eligible).sum())}/12.",
    ]
    return "\n".join(lines) + "\n"


def _classify_case(
    summaries: list[dict[str, Any]],
    grid: pd.DataFrame,
    validity: pd.DataFrame,
) -> tuple[str, str]:
    pooled = next(row for row in summaries if row["scope"] == "pooled_N50_N100")
    if pooled["common_intersection_nonempty"]:
        return "CASE A", "a nontrivial common fixed-c intersection exists"
    invalid_share = float((~validity["true_rank_postrefit_numerically_valid"]).mean())
    if invalid_share >= 0.5:
        return "CASE D", "true-rank numerical invalidity is too extensive"
    eligible_by_n = validity.groupby("N")[
        "true_rank_postrefit_numerically_valid"
    ].sum()
    pivot = grid.pivot(index="c_kappa", columns="N", values="Exact").dropna()
    substantial = (
        (pivot[50] >= max(1, math.ceil(float(eligible_by_n.loc[50]) / 2.0)))
        & (pivot[100] >= max(1, math.ceil(float(eligible_by_n.loc[100]) / 2.0)))
    ).any()
    if substantial:
        return "CASE B", "no common intersection, but one fixed-c tradeoff recovers truth substantially at both sizes"
    return "CASE C", "truth-optimal fixed-c requirements are incompatible across sizes"


def _report(
    validity: pd.DataFrame,
    penalty: pd.DataFrame,
    breaks: pd.DataFrame,
    intervals: pd.DataFrame,
    summaries: list[dict[str, Any]],
    grid: pd.DataFrame,
    failures: pd.DataFrame,
    split: pd.DataFrame,
    max_identity_error: float,
    case: str,
    case_reason: str,
) -> str:
    validity_summary = validity.groupby("N").agg(
        present=("true_rank_present_in_candidate_set", "sum"),
        valid=("true_rank_postrefit_numerically_valid", "sum"),
        invalid=("true_rank_postrefit_numerically_invalid", "sum"),
    )
    baseline = grid.loc[np.isclose(grid["c_kappa"], 1.0)]
    especially = breaks.loc[breaks["especially_requested_comparison"]].copy()
    break_summary = especially.groupby(["N", "competitor_rank"], sort=True).agg(
        available=("c_break_even", "size"),
        minimum=("c_break_even", "min"),
        median=("c_break_even", "median"),
        maximum=("c_break_even", "max"),
    )
    exact_pivot = grid.pivot(index="c_kappa", columns="N", values="Exact").dropna()
    exact_pivot["minimum_across_sizes"] = exact_pivot[[50, 100]].min(axis=1)
    exact_pivot["total_across_sizes"] = exact_pivot[50] + exact_pivot[100]
    best_tradeoff = exact_pivot.sort_values(
        ["minimum_across_sizes", "total_across_sizes"], ascending=False
    ).iloc[0]
    split_summary = split.loc[split["record_scope"] == "replication_summary"]
    lines = [
        "# Locked Revision-9 no-go diagnosis",
        "",
        "This is an offline analysis of already-computed preflight records. No simulation, fitting, rank selection, split fit, or Riesz solve was run.",
        "Accepted source preflight commit: `9f0778d2503a66d913274ac41ae2b6e94257a5b9`.",
        "",
        "## Candidate validity and locked IC",
        "",
    ]
    for n, row in validity_summary.iterrows():
        lines.append(
            f"- N={n}: truth present {int(row.present)}/12; valid truth post-refit {int(row.valid)}/12; invalid {int(row.invalid)}/12."
        )
    for row in validity.loc[
        validity["true_rank_postrefit_numerically_invalid"]
    ].itertuples():
        lines.append(
            f"  - {row.semantic_replication_id}: `{row.invalid_reasons}`."
        )
    lines += [
        f"- Maximum locked-IC reconstruction error over valid candidates: `{max_identity_error:.3g}`.",
        "- At c=1, both sizes select (0,0,0) in all 12 saved replications.",
        f"- Best simultaneous saved-grid tradeoff occurs at c={float(best_tradeoff.name):.8g}: Exact={int(best_tradeoff[50])}/12 at N=50 and {int(best_tradeoff[100])}/12 at N=100; this is not substantial at both sizes.",
        "",
        "## Especially requested pairwise break-even values",
        "",
        "| N | Competitor | Available | Minimum c | Median c | Maximum c |",
        "|--:|:--|--:|--:|--:|--:|",
    ]
    for (n, competitor), row in break_summary.iterrows():
        lines.append(
            f"| {n} | {competitor} | {int(row['available'])} | {row['minimum']:.8g} | {row['median']:.8g} | {row['maximum']:.8g} |"
        )
    lines += [
        "",
        "## Common-c result",
        "",
    ]
    for row in summaries:
        lines.append(
            f"- {row['scope']}: common [{row['common_c_lower']:.8g}, {row['common_c_upper']:.8g}], nonempty={row['common_intersection_nonempty']}; maximum overlap {row['maximum_simultaneous_overlap']}/{row['eligible_replications']}."
        )
    lines += [
        f"- Classification: **{case}** -- {case_reason}.",
        "",
        "## Baseline offline rank behavior",
        "",
    ]
    for row in baseline.itertuples():
        lines.append(
            f"- N={row.N}, c=1: Exact={row.Exact}, Under only={row.Under_only}, Over only={row.Over_only}, Mixed={row.Mixed}; distribution {row.selected_rank_distribution}."
        )
    lines += [
        "",
        "## Supplied-rank failures and split interiority",
        "",
        f"- The four N=50 full-fit failures contain {len(failures)} starts. Every start converged, remained feasible, retained rank (1,1,1), and failed only because its residual was slightly above the locked `1e-6` stationarity threshold. The primary classification is KKT/optimality failure, not start disagreement.",
        "",
        "| Failed semantic ID | Residual range | Maximum objective gap | Maximum coefficient | Classification |",
        "|:--|:--|--:|--:|:--|",
    ]
    for semantic_id, group in failures.groupby("semantic_replication_id", sort=True):
        lines.append(
            f"| {semantic_id} | {group.stationarity_or_KKT_residual.min():.8g} to {group.stationarity_or_KKT_residual.max():.8g} | {group.objective_gap_to_best.max():.8g} | {group.max_absolute_coefficient.max():.8g} | KKT/optimality failure |"
        )
    lines += [
        "",
        f"- Boundary-active split fits: {int((split.record_scope == 'boundary_fit').sum())}.",
        f"- Successful N=50 fixed-rank replications with zero/one/multiple active splits: {(split_summary.boundary_count_category == '0').sum()}/{(split_summary.boundary_count_category == '1').sum()}/{(split_summary.boundary_count_category == 'multiple').sum()}.",
        "- Five successful replications have at least one boundary-active split and all eight broad targets in each are suppressed by `boundary_interiority_failure`. The other three have no boundary-active split but fail the locked split numerical checks, producing `split_fit_failure`. Hence broad inference is 0/96.",
        "- The locked records do not persist the active matrix identity or entry-level active-set count. Those fields are explicitly marked unavailable; they are not inferred.",
        "",
        "## Fixed-positive-multiplier theorem implication",
        "",
        "For any fixed `c_kappa > 0`, the ratio of `zeta_NT` to `kappa_NT (N+T)/(NT)` is `1/{c_kappa log(NT)} -> 0`. A fixed positive multiplier also leaves the maintained vanishing-rate condition unchanged. Therefore c=1 is a normalization rather than a constant uniquely implied by the proof. Any revision must nevertheless freeze the multiplier before a fresh independent validation experiment; post-outcome selection would be inappropriate.",
        "",
        "## Statistical versus numerical decomposition",
        "",
        "- Supplied rank: the N=100 estimator/inference path is fully retained; N=50 has four numerical full-fit failures and universal broad-target split-interiority/numerical failure. Gram/Riesz diagnostics are reliable wherever reached.",
        "- Selected rank: candidate coverage is complete, but only 6/12 N=50 true-rank post-refits are numerically eligible. At N=100 all 12 are eligible. Conditional on eligibility, the locked IC penalty scale drives selection to zero rank at c=1, which removes tangent support and suppresses all selected-rank inference.",
        "- Paper-level issue: the finite-sample normalization of the paper IC. Numerical issue for the same estimator: N=50 stationarity acceptance and split interiority reliability. Neither is changed here.",
        "",
        "## Next paper-level decision",
        "",
        "The author should decide whether Revision 9 should retain its c=1 normalization or be revised to pre-specify a different fixed positive multiplier and a fresh validation protocol. This audit does not select a constant. Independently, the author should decide whether the N=50 numerical/interiority evidence is acceptable for the claimed design or requires a paper-level change in supported sample sizes or numerical assumptions before any code change.",
        "",
        "No medium or production run is authorized by this report.",
        "",
        "## Penalty scale snapshot",
        "",
        "| N | b_NT | kappa_base | 0->1 | 1->2 | 2->3 | truth-zero |",
        "|--:|--:|--:|--:|--:|--:|--:|",
    ]
    for row in penalty.itertuples():
        lines.append(
            f"| {row.N} | {row.b_NT:.8g} | {row.kappa_base:.8g} | {row.penalty_increment_rank_0_to_1:.8g} | {row.penalty_increment_rank_1_to_2:.8g} | {row.penalty_increment_rank_2_to_3:.8g} | {row.total_penalty_truth_minus_zero_c_1:.8g} |"
        )
    return "\n".join(lines) + "\n"


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--input-root",
        type=Path,
        default=Path("results/mc/preflight_revision9_locked"),
    )
    parser.add_argument(
        "--output-root",
        type=Path,
        default=Path("results/mc/audit/revision9_nogo_diagnosis"),
    )
    args = parser.parse_args()
    root = args.input_root
    fixed_root = root / "fixed_rank" / "f87b622f889053fe"
    selected_root = root / "selected_rank" / "aa152561964a7ec3"
    args.output_root.mkdir(parents=True, exist_ok=True)

    candidates = pd.read_parquet(root / "candidate_ic_diagnostics.parquet")
    fits = pd.read_parquet(root / "fit_diagnostics.parquet")
    selected_rank_rows = _read_chunks(selected_root / "rank")
    fixed_rank_rows = _read_chunks(fixed_root / "rank")
    fixed_raw = _read_chunks(fixed_root / "raw")

    if candidates["semantic_replication_id"].nunique() != 24:
        raise AssertionError("offline audit expected exactly 24 locked selected-rank realizations")
    validity, breaks, intervals, max_identity_error = _ic_outputs(
        candidates, selected_rank_rows
    )
    penalty = _penalty_scale()
    summaries = [
        _intersection_summary(intervals.loc[intervals.N == 50], "N50"),
        _intersection_summary(intervals.loc[intervals.N == 100], "N100"),
        _intersection_summary(intervals, "pooled_N50_N100"),
    ]
    grid = _offline_grid(candidates, _grid_values(breaks, intervals))
    failures = _fixed_failure_details(fixed_rank_rows, fits)
    split = _split_boundary_details(fits, fixed_raw)
    case, case_reason = _classify_case(summaries, grid, validity)

    penalty.to_csv(args.output_root / "ic_penalty_scale.csv", index=False)
    breaks.to_csv(args.output_root / "pairwise_break_even.csv", index=False)
    intervals.to_csv(args.output_root / "truth_optimal_intervals.csv", index=False)
    grid.to_csv(args.output_root / "offline_existing_candidate_grid.csv", index=False)
    validity.to_csv(args.output_root / "true_rank_postrefit_validity.csv", index=False)
    failures.to_csv(args.output_root / "fixed_n50_failure_details.csv", index=False)
    split.to_csv(args.output_root / "split_boundary_details.csv", index=False)
    (args.output_root / "common_c_intersection.md").write_text(
        _common_markdown(summaries, intervals), encoding="utf-8"
    )
    (args.output_root / "theorem_penalty_implications.md").write_text(
        _theorem_text(), encoding="utf-8"
    )
    (args.output_root / "revision9_nogo_diagnosis_report.md").write_text(
        _report(
            validity,
            penalty,
            breaks,
            intervals,
            summaries,
            grid,
            failures,
            split,
            max_identity_error,
            case,
            case_reason,
        ),
        encoding="utf-8",
    )
    print(
        json.dumps(
            {
                "case": case,
                "max_ic_reconstruction_error": max_identity_error,
                "true_rank_valid_by_N": validity.groupby("N")[
                    "true_rank_postrefit_numerically_valid"
                ].sum().astype(int).to_dict(),
                "boundary_split_fits": int(
                    (split["record_scope"] == "boundary_fit").sum()
                ),
                "output_root": str(args.output_root),
            },
            indent=2,
        )
    )


if __name__ == "__main__":
    main()
