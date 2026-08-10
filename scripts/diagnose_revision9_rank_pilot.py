"""Offline diagnosis of the locked Revision-9 thresholded cap pilot.

This script reads saved rank records only. It does not import or call DGP,
estimation, rank-selection, split-fit, Riesz, or inference routines.
"""

from __future__ import annotations

import argparse
import json
from collections import Counter
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd

EXPECTED_REALIZATIONS = 24
MATRIX_NAMES = ("A", "B", "H")
RECOVERY_CLASSES = ("Exact", "Under only", "Over only", "Mixed")


def _vector(value: Any) -> tuple[int, ...]:
    if isinstance(value, str):
        value = json.loads(value)
    return tuple(int(item) for item in value)


def _vector_text(value: tuple[int, ...]) -> str:
    return "(" + ",".join(str(item) for item in value) + ")"


def _classification(truth: tuple[int, ...], estimate: tuple[int, ...]) -> str:
    below = [left < right for left, right in zip(estimate, truth, strict=True)]
    above = [left > right for left, right in zip(estimate, truth, strict=True)]
    if estimate == truth:
        return "Exact"
    if any(below) and any(above):
        return "Mixed"
    if any(below):
        return "Under only"
    return "Over only"


def _one_coordinate_neighbor(left: tuple[int, ...], right: tuple[int, ...]) -> bool:
    return sum(abs(a - b) for a, b in zip(left, right, strict=True)) == 1


def _load_rank_rows(selected_root: Path) -> pd.DataFrame:
    files = sorted((selected_root / "rank").glob("*.parquet"))
    if not files:
        raise FileNotFoundError(f"no locked rank parquet files under {selected_root}")
    rows = pd.concat([pd.read_parquet(path) for path in files], ignore_index=True)
    if len(rows) != EXPECTED_REALIZATIONS:
        raise AssertionError(f"expected 24 locked rank rows, found {len(rows)}")
    if rows["semantic_replication_id"].nunique() != EXPECTED_REALIZATIONS:
        raise AssertionError("locked rank rows do not have 24 unique semantic IDs")
    if not rows["candidate_coverage"].astype(bool).all():
        raise AssertionError("locked candidate coverage is not 24/24")
    return rows.sort_values(["N", "dgp", "replication"]).reset_index(drop=True)


def _strata(rows: pd.DataFrame) -> list[tuple[str, str, pd.DataFrame]]:
    output = [("pooled", "all", rows)]
    output.extend(("N", str(value), group) for value, group in rows.groupby("N"))
    output.extend(("DGP", str(value), group) for value, group in rows.groupby("dgp"))
    return output


def _rank_summary(records: pd.DataFrame) -> pd.DataFrame:
    output: list[dict[str, Any]] = []
    for stratum_type, stratum_value, group in _strata(records):
        total = len(group)
        distribution = Counter(group["thresholded_cap_pilot_rank"])
        for rank, count in sorted(distribution.items()):
            output.append(
                {
                    "summary_type": "rank_distribution",
                    "stratum_type": stratum_type,
                    "stratum_value": stratum_value,
                    "category": rank,
                    "count": count,
                    "total": total,
                    "rate": count / total,
                }
            )
        recovery = Counter(group["recovery_classification"])
        for category in RECOVERY_CLASSES:
            count = recovery[category]
            output.append(
                {
                    "summary_type": "recovery_classification",
                    "stratum_type": stratum_type,
                    "stratum_value": stratum_value,
                    "category": category,
                    "count": count,
                    "total": total,
                    "rate": count / total,
                }
            )
    return pd.DataFrame(output)


def _candidate_source_record(row: pd.Series, diagnostics: dict[str, Any]) -> dict[str, Any]:
    truth = _vector(row["true_rank_vector"])
    pilot = _vector(row["rank_cap_thresholded_vector"])
    path = [_vector(value) for value in diagnostics["nuclear_path_rank_proposals"]]
    direct_cap = pilot == truth
    direct_nuclear = truth in path
    neighbor_cap = _one_coordinate_neighbor(pilot, truth)
    neighbor_nuclear = any(_one_coordinate_neighbor(proposal, truth) for proposal in path)
    mechanisms = []
    if direct_cap:
        mechanisms.append("thresholded_cap_pilot_direct")
    if direct_nuclear:
        mechanisms.append("nuclear_path_direct")
    if neighbor_cap:
        mechanisms.append("neighbor_of_cap_pilot")
    if neighbor_nuclear:
        mechanisms.append("neighbor_of_nuclear_proposal")
    attribution = mechanisms[0] if len(mechanisms) == 1 else "multiple_sources"
    return {
        "semantic_replication_id": row["semantic_replication_id"],
        "dgp": int(row["dgp"]),
        "N": int(row["N"]),
        "T": int(row["T"]),
        "replication": int(row["replication"]),
        "true_rank": _vector_text(truth),
        "candidate_coverage": bool(row["candidate_coverage"]),
        "thresholded_cap_pilot_direct": direct_cap,
        "nuclear_path_direct": direct_nuclear,
        "neighbor_of_cap_pilot": neighbor_cap,
        "neighbor_of_nuclear_proposal": neighbor_nuclear,
        "multiple_sources": len(mechanisms) > 1,
        "coverage_attribution": attribution,
        "all_identified_mechanisms": json.dumps(mechanisms),
        "saved_true_rank_sources": row["true_rank_sources"],
    }


def _nuclear_record(row: pd.Series, diagnostics: dict[str, Any]) -> dict[str, Any]:
    truth = _vector(row["true_rank_vector"])
    proposals = [_vector(value) for value in diagnostics["nuclear_path_rank_proposals"]]
    counts = Counter(proposals)
    modal_frequency = max(counts.values())
    modal = sorted(rank for rank, count in counts.items() if count == modal_frequency)
    return {
        "record_type": "replication",
        "stratum_type": "replication",
        "stratum_value": row["semantic_replication_id"],
        "semantic_replication_id": row["semantic_replication_id"],
        "dgp": int(row["dgp"]),
        "N": int(row["N"]),
        "T": int(row["T"]),
        "replication": int(row["replication"]),
        "true_rank": _vector_text(truth),
        "path_proposal_count": len(proposals),
        "distinct_path_rank_count": len(counts),
        "any_path_proposal_equals_truth": truth in proposals,
        "truth_proposal_frequency": counts[truth],
        "modal_path_rank": json.dumps([_vector_text(value) for value in modal]),
        "modal_path_rank_frequency": modal_frequency,
        "path_rank_proposals": json.dumps([_vector_text(value) for value in proposals]),
    }


def _nuclear_summary(replication_rows: pd.DataFrame) -> pd.DataFrame:
    output = replication_rows.to_dict("records")
    for stratum_type, stratum_value, group in _strata(replication_rows):
        modal_counter: Counter[str] = Counter()
        for proposals in group["path_rank_proposals"]:
            modal_counter.update(json.loads(proposals))
        maximum = max(modal_counter.values())
        modal = sorted(rank for rank, count in modal_counter.items() if count == maximum)
        output.append(
            {
                "record_type": "summary",
                "stratum_type": stratum_type,
                "stratum_value": stratum_value,
                "semantic_replication_id": "",
                "dgp": np.nan,
                "N": np.nan,
                "T": np.nan,
                "replication": np.nan,
                "true_rank": "(1,1,1)",
                "path_proposal_count": int(group["path_proposal_count"].sum()),
                "distinct_path_rank_count": len(modal_counter),
                "any_path_proposal_equals_truth": bool(
                    group["any_path_proposal_equals_truth"].all()
                ),
                "truth_proposal_frequency": int(group["truth_proposal_frequency"].sum()),
                "modal_path_rank": json.dumps(modal),
                "modal_path_rank_frequency": maximum,
                "path_rank_proposals": "",
            }
        )
    return pd.DataFrame(output)


def _theory_text() -> str:
    return """# Direct thresholded-pilot rank consistency

## Audit basis and conclusion

The maintained assumptions supplied for this audit are sufficient to imply consistency of the
thresholded pilot rank vector. This is an asymptotic implication; it does not make the observed
locked pilot ranks accurate at N=50 or N=100.

For each coefficient matrix M in the fixed collection
`{A^(1),...,A^(p),B^(1),...,B^(K),H}`, let M0 have fixed rank r_M and let M_tilde be the saved
rank-cap pilot conceptually covered by the maintained pilot lemma. The assumptions used are:

1. `max_M ||M_tilde-M0||_op = o_p(tau_NT)` (the maintained pilot operator-norm error condition);
2. `tau_NT=o(sqrt(NT))` and `tau_NT>0`;
3. for every positive-rank M0, `sigma_{r_M}(M0) >= c_M sqrt(NT)` with fixed `c_M>0` with
   probability tending to one, equivalently the maintained strong-factor singular-value condition;
4. p, K, the number of matrices, and all true ranks are fixed.

## Proof

Fix one matrix M. Weyl's singular-value perturbation inequality gives, for every j,

`|sigma_j(M_tilde)-sigma_j(M0)| <= ||M_tilde-M0||_op`.

Write `e_NT=||M_tilde-M0||_op`. From assumption 1,
`P(e_NT<tau_NT/2)->1`. If `r_M>0`, assumptions 2-3 imply
`sigma_{r_M}(M0)/tau_NT -> infinity` in probability, so
`P(sigma_{r_M}(M0)>2 tau_NT)->1`.

On the intersection of these events, for every `j<=r_M`,

`sigma_j(M_tilde) >= sigma_{r_M}(M0)-e_NT > 3 tau_NT/2 > tau_NT`.

Thus no true positive singular value is thresholded away. For every `j>r_M`,
`sigma_j(M0)=0`, and Weyl gives

`sigma_j(M_tilde) <= e_NT < tau_NT/2 < tau_NT`.

Thus no population-zero singular value exceeds the threshold. If `r_M=0`, only the second
argument is needed and the thresholded rank is zero. Hence the thresholded rank of M_tilde equals
`r_M` with probability tending to one.

The same argument applies separately to every A^(ell), B^(k), and H. Because their number is
fixed, a union bound over the finitely many failure events yields joint recovery of the complete
rank vector with probability tending to one.

## Manuscript implications

The existing pilot operator-error lemma and the manuscript's Weyl perturbation step supply the
entire mathematical argument above. The manuscript itself is not stored in this repository, so
this offline audit cannot responsibly assign exact lemma/proposition numbers or quote their
wording. Before editing, the author should map the maintained pilot lemma and strong-factor
assumption to their Revision-9 labels.

If the pilot became the final rank estimator, the algorithm definition, rank-selection
consistency theorem, and every downstream oracle/inference theorem that currently invokes the
IC-selected rank would need corresponding revision. IC-specific consistency arguments,
`kappa_NT` separation/rate material, candidate post-refit global-gap assumptions used only for IC
selection, and local-completion consistency material could potentially be removed. Candidate
fits might remain computational diagnostics, but they would no longer define the final rank.
No manuscript change is made by this audit.
"""


def _report_text(
    records: pd.DataFrame, sources: pd.DataFrame, nuclear: pd.DataFrame
) -> str:
    rank_counts = records["thresholded_cap_pilot_rank"].value_counts().sort_index()
    class_counts = records["recovery_classification"].value_counts()
    boundary = int(records["pilot_boundary_active"].sum())
    direct_cap = int(sources["thresholded_cap_pilot_direct"].sum())
    direct_nuclear = int(sources["nuclear_path_direct"].sum())
    neighbor_cap = int(sources["neighbor_of_cap_pilot"].sum())
    neighbor_nuclear = int(sources["neighbor_of_nuclear_proposal"].sum())
    multiple_sources = int(sources["multiple_sources"].sum())
    pooled_path_counts: Counter[str] = Counter()
    for proposals in nuclear["path_rank_proposals"]:
        pooled_path_counts.update(json.loads(proposals))
    pooled_modal_frequency = max(pooled_path_counts.values())
    pooled_modal = sorted(
        rank for rank, count in pooled_path_counts.items() if count == pooled_modal_frequency
    )
    distribution = ", ".join(f"{rank}: {count}" for rank, count in rank_counts.items())
    return f"""# Revision-9 thresholded rank-pilot diagnosis

This is an offline audit of 24 already-computed locked selected-rank records. It ran no DGP,
screening, fitting, IC selection, split fitting, Riesz solve, inference, or Monte Carlo.

## Cap-pilot ranks and numerical status

The pooled thresholded-rank distribution is: {distribution}. Exact recovery is 0/24; Over only
is {class_counts.get('Over only', 0)}/24; Under only is {class_counts.get('Under only', 0)}/24;
Mixed is {class_counts.get('Mixed', 0)}/24. At N=50 the counts are 0 Exact, 1 Under only, 11 Over
only, 0 Mixed. At N=100 they are 0 Exact, 0 Under only, 12 Over only, 0 Mixed. Every DGP has zero
exact recovery.

All 24 pilots are recorded as converged, objective-stable, and stationarity/KKT-valid under the
locked interior (`1e-6`) or constrained (`1e-4`) criterion. {boundary}/24 are boundary-active by
the saved coefficient-envelope diagnostic. This is therefore not a recorded pilot
numerical-validity failure; the accepted rank-adaptive pilot outputs are mostly too large.

## Candidate coverage source

Truth is in the candidates 24/24. Direct cap-pilot attribution is {direct_cap}/24. A nuclear-path
proposal equals truth in {direct_nuclear}/24. Truth is also a one-coordinate neighbor of the cap
pilot in {neighbor_cap}/24 and of at least one nuclear proposal in {neighbor_nuclear}/24. The
coverage attribution is nuclear-direct only in {24 - multiple_sources}/24 and multiple-source in
{multiple_sources}/24. The complete per-replication mechanisms are in
`candidate_coverage_source.csv`. The direct pilot did not account for any of the observed 24/24
coverage.

## Singular-value threshold margins

The locked rank rows save `tau_NT` and the thresholded vector but do not save the accepted pilot's
matrix-specific singular values. The fit-diagnostic table retains only aggregate `sigma_1` and
`sigma_r` for individual route fits; it cannot identify sigma_2 for each accepted A, B, and H.
Therefore sigma_1/tau, sigma_2/tau, signal/noise margins, and near-threshold counts are unavailable
without rerunning a prohibited fit. `cap_pilot_singular_value_margins.csv` records this limitation
for all 72 matrix-realization pairs rather than fabricating values.

## Pilot versus final IC

The pilot is incorrect and the final IC is incorrect in 24/24. The other three requested cells
(correct pilot/incorrect IC, incorrect pilot/correct IC, correct pilot/correct IC) are all zero.
Thus the IC did not destroy correct information already present in the pilot: the pilot overfit
while the locked IC selected `(0,0,0)` in every realization.

## Nuclear path

At least one path proposal equals `(1,1,1)` in 24/24. This makes the nuclear path an effective
screening source in these records, not a demonstrated final estimator. Modal proposals and full
saved proposal sequences are reported in `nuclear_path_rank_summary.csv`. Pooled over all 528
path positions, the modal rank is {json.dumps(pooled_modal)} with frequency
{pooled_modal_frequency}.

## Theory and evidence classification

The maintained `o_p(tau_NT)` pilot operator-error condition, strong-factor singular values of
order `sqrt(NT)`, `tau_NT=o(sqrt(NT))`, and fixed matrix count directly imply joint thresholded
rank consistency by Weyl's inequality. The complete argument is in
`direct_pilot_rank_theory.md`.

The evidence classification is **CASE P3**: thresholded cap-pilot recovery is poor (0/24 exact),
so removing the IC would not solve the locked finite-sample problem. This does not contradict the
asymptotic proof; it shows that its separation regime is not visible in these two saved sizes.
Missing singular spectra prevent a margin diagnosis but do not obscure the observed 0/24 rank
recovery.

## Paper-level options

1. Retaining the current IC preserves Revision 9 but must acknowledge its finite-sample NO-GO.
2. Using pilot ranks has a direct asymptotic proof, but the locked evidence rejects it as the
   current finite-sample solution; it should not be implemented on this record.
3. A new rank selector requires new theory and a separately locked validation design. It is the
   only option aimed at resolving both observed failures, but is a future paper-level project.

Recommended next decision: retain the NO-GO while the author decides whether to keep Option 1 or
undertake Option 3. Do not adopt Option 2 from these data and do not search for another fixed
`c_kappa`.

## Separate N=50 supplied-rank issue

The existing locked audit reports four narrow full-fit stationarity failures, nine boundary-active
split fits, broad-target inference failure at N=50, and clean N=100 supplied-rank performance.
That supports reporting N=50 as a small-sample stress design rather than a headline fixed-rank
design, without deleting it or changing B, tolerances, or the estimator. This issue is distinct
from selected-rank failure, which occurs at both N=50 and N=100.
"""


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--selected-root",
        type=Path,
        default=Path(
            "results/mc/preflight_revision9_locked/selected_rank/aa152561964a7ec3"
        ),
    )
    parser.add_argument(
        "--output-root",
        type=Path,
        default=Path("results/mc/audit/revision9_rank_pilot_diagnosis"),
    )
    args = parser.parse_args()
    rows = _load_rank_rows(args.selected_root)
    config = json.loads((args.selected_root / "resolved_config.json").read_text())[
        "config"
    ]
    stationarity_tol = float(config["estimation"]["stationarity_tol"])
    constrained_tol = float(config["estimation"]["constrained_kkt_tolerance"])
    constraint_tol = float(config["estimation"]["constraint_tolerance"])
    coefficient_bound = float(config["estimation"]["coefficient_bound"])

    record_rows: list[dict[str, Any]] = []
    source_rows: list[dict[str, Any]] = []
    margin_rows: list[dict[str, Any]] = []
    pilot_ic_rows: list[dict[str, Any]] = []
    nuclear_rows: list[dict[str, Any]] = []
    for _, row in rows.iterrows():
        diagnostics = json.loads(row["rank_diagnostics_json"])
        pilot = diagnostics["cap_pilot"]
        truth = _vector(row["true_rank_vector"])
        thresholded = _vector(row["rank_cap_thresholded_vector"])
        final_rank = _vector(row["selected_rank_vector"])
        boundary = float(row["cap_pilot_max_envelope_ratio"]) >= (
            1.0 - constraint_tol / coefficient_bound
        )
        stationarity_limit = constrained_tol if boundary else stationarity_tol
        stationarity_pass = (
            float(row["cap_pilot_stationarity_residual"]) <= stationarity_limit
        )
        valid = (
            bool(row["cap_pilot_converged"])
            and bool(row["cap_pilot_objective_stability_pass"])
            and stationarity_pass
        )
        recovery = _classification(truth, thresholded)
        record_rows.append(
            {
                "semantic_replication_id": row["semantic_replication_id"],
                "dgp": int(row["dgp"]),
                "N": int(row["N"]),
                "T": int(row["T"]),
                "replication": int(row["replication"]),
                "true_rank": _vector_text(truth),
                "cap_pilot_numerical_rank_before_thresholding": _vector_text(
                    _vector(row["cap_pilot_numerical_rank_before_thresholding"])
                ),
                "thresholded_cap_pilot_rank": _vector_text(thresholded),
                "recovery_classification": recovery,
                "pilot_numerically_valid": valid,
                "pilot_converged": bool(row["cap_pilot_converged"]),
                "pilot_objective_stability_pass": bool(
                    row["cap_pilot_objective_stability_pass"]
                ),
                "pilot_boundary_active": boundary,
                "pilot_max_envelope_ratio": float(row["cap_pilot_max_envelope_ratio"]),
                "pilot_stationarity_or_KKT_residual": float(
                    row["cap_pilot_stationarity_residual"]
                ),
                "pilot_stationarity_or_KKT_limit": stationarity_limit,
                "pilot_stationarity_or_KKT_pass": stationarity_pass,
                "pilot_objective": float(row["cap_pilot_best_valid_objective"]),
                "pilot_acceptance_basis": pilot["final_pilot_acceptance_basis"],
                "pilot_stable_thresholded_ranks_agree": bool(
                    pilot["stable_final_thresholded_ranks_agree"]
                ),
            }
        )
        source_rows.append(_candidate_source_record(row, diagnostics))
        nuclear_rows.append(_nuclear_record(row, diagnostics))
        pilot_exact = thresholded == truth
        ic_exact = final_rank == truth
        pilot_ic_rows.append(
            {
                "semantic_replication_id": row["semantic_replication_id"],
                "dgp": int(row["dgp"]),
                "N": int(row["N"]),
                "T": int(row["T"]),
                "replication": int(row["replication"]),
                "true_rank": _vector_text(truth),
                "thresholded_cap_pilot_rank": _vector_text(thresholded),
                "final_IC_selected_rank": _vector_text(final_rank),
                "pilot_exact": pilot_exact,
                "IC_exact": ic_exact,
                "pilot_vs_IC_classification": (
                    ("correct" if pilot_exact else "incorrect")
                    + " pilot -> "
                    + ("correct" if ic_exact else "incorrect")
                    + " IC"
                ),
            }
        )
        tau = float(diagnostics["threshold"])
        for matrix_index, matrix_name in enumerate(MATRIX_NAMES):
            margin_rows.append(
                {
                    "record_type": "matrix",
                    "semantic_replication_id": row["semantic_replication_id"],
                    "dgp": int(row["dgp"]),
                    "N": int(row["N"]),
                    "T": int(row["T"]),
                    "replication": int(row["replication"]),
                    "matrix": matrix_name,
                    "true_rank": truth[matrix_index],
                    "tau_NT": tau,
                    "thresholded_rank": thresholded[matrix_index],
                    "sigma_1_hat": np.nan,
                    "sigma_2_hat": np.nan,
                    "sigma_1_hat_over_tau_NT": np.nan,
                    "sigma_2_hat_over_tau_NT": np.nan,
                    "sigma_1_hat_minus_tau_NT": np.nan,
                    "tau_NT_minus_sigma_2_hat": np.nan,
                    "near_threshold_classification": "unavailable",
                    "singular_value_margin_available": False,
                    "minimum_signal_above_threshold_ratio": np.nan,
                    "median_signal_above_threshold_ratio": np.nan,
                    "minimum_noise_below_threshold_margin": np.nan,
                    "near_threshold_classification_count": np.nan,
                    "unavailability_reason": (
                        "accepted pilot matrix-specific singular values were not saved"
                    ),
                }
            )

    records = pd.DataFrame(record_rows)
    sources = pd.DataFrame(source_rows)
    nuclear_replications = pd.DataFrame(nuclear_rows)
    margins = pd.DataFrame(margin_rows)
    for n_value, group in margins.groupby("N"):
        margins.loc[len(margins)] = {
            "record_type": "N_summary",
            "semantic_replication_id": "",
            "dgp": np.nan,
            "N": n_value,
            "T": n_value,
            "replication": np.nan,
            "matrix": "all",
            "true_rank": 1,
            "tau_NT": float(group["tau_NT"].iloc[0]),
            "thresholded_rank": np.nan,
            "near_threshold_classification": "unavailable",
            "singular_value_margin_available": False,
            "minimum_signal_above_threshold_ratio": np.nan,
            "median_signal_above_threshold_ratio": np.nan,
            "minimum_noise_below_threshold_margin": np.nan,
            "near_threshold_classification_count": np.nan,
            "unavailability_reason": (
                "minimum/median margins and near-threshold counts cannot be computed"
            ),
        }

    args.output_root.mkdir(parents=True, exist_ok=True)
    records.to_csv(args.output_root / "cap_pilot_rank_records.csv", index=False)
    _rank_summary(records).to_csv(
        args.output_root / "cap_pilot_rank_summary.csv", index=False
    )
    sources.to_csv(args.output_root / "candidate_coverage_source.csv", index=False)
    margins.to_csv(
        args.output_root / "cap_pilot_singular_value_margins.csv", index=False
    )
    pd.DataFrame(pilot_ic_rows).to_csv(
        args.output_root / "pilot_vs_ic_rank.csv", index=False
    )
    _nuclear_summary(nuclear_replications).to_csv(
        args.output_root / "nuclear_path_rank_summary.csv", index=False
    )
    (args.output_root / "direct_pilot_rank_theory.md").write_text(
        _theory_text(), encoding="utf-8"
    )
    (args.output_root / "revision9_rank_pilot_diagnosis_report.md").write_text(
        _report_text(records, sources, nuclear_replications), encoding="utf-8"
    )

    if int((records["recovery_classification"] == "Exact").sum()) != 0:
        raise AssertionError("unexpected locked pilot exact recovery")
    if int(sources["thresholded_cap_pilot_direct"].sum()) != 0:
        raise AssertionError("unexpected direct cap-pilot candidate coverage")
    if not sources["nuclear_path_direct"].all():
        raise AssertionError("truth is not directly on every saved nuclear path")
    print(f"wrote offline rank-pilot diagnosis to {args.output_root}")


if __name__ == "__main__":
    main()
