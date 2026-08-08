"""Reconstruct Revision-8 IC break-even diagnostics from saved preflight fits only."""

from __future__ import annotations

import argparse
import json
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd

from dynamic_panel_econ.ic_break_even import (
    ICPoint,
    dimension_penalty,
    optimality_interval,
    pairwise_break_even,
    rank_increment_dimension,
    select_ic_candidate,
)
from dynamic_panel_econ.rank_selection import revision8_kappa

GRID = (1e-8, 3e-8, 1e-7, 3e-7, 1e-6, 3e-6, 1e-5, 3e-5, 1e-4, 3e-4, 1e-3)
TRUE_RANK = (1, 1, 1)


def _rank_text(ranks: tuple[int, ...] | list[int] | None) -> str:
    return "" if ranks is None else "(" + ",".join(str(int(x)) for x in ranks) + ")"


def _load_records(root: Path) -> pd.DataFrame:
    files = sorted(root.glob("*.parquet"))
    if not files:
        raise FileNotFoundError(f"no saved rank chunks found under {root}")
    frame = pd.concat((pd.read_parquet(path) for path in files), ignore_index=True)
    frame = frame.loc[frame["status"].eq("success")].copy()
    if frame.empty:
        raise ValueError("no successful saved preflight replications")
    return frame.sort_values(["dgp", "replication"], kind="stable")


def _points(records: list[dict[str, Any]]) -> list[ICPoint]:
    return [
        ICPoint(tuple(row["ranks"]), max(2.0 * float(row["objective"]), np.finfo(float).tiny), int(row["dimension"]))
        for row in records
        if bool(row["valid"])
    ]


def _reconstruct(frame: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    interval_rows: list[dict[str, Any]] = []
    grid_rows: list[dict[str, Any]] = []
    for _, saved in frame.iterrows():
        diagnostics = json.loads(saved["rank_diagnostics_json"])
        records = diagnostics["candidate_records"]
        required = {"ranks", "objective", "dimension", "valid", "sources", "ic"}
        incomplete = [index for index, row in enumerate(records) if not required <= row.keys()]
        if incomplete:
            raise ValueError(
                f"DGP {saved['dgp']} replication {saved['replication']} has incomplete "
                f"candidate records at indices {incomplete}"
            )
        if len(records) != int(saved["candidate_count_final"]):
            raise ValueError(
                f"DGP {saved['dgp']} replication {saved['replication']} did not persist "
                "the complete final candidate set"
            )
        valid = _points(records)
        truth_record = next((row for row in records if tuple(row["ranks"]) == TRUE_RANK), None)
        truth_valid = bool(truth_record is not None and truth_record["valid"])
        n, t = int(saved["N"]), int(saved["T"])
        base_per_dimension = revision8_kappa(n, t) / (n * t)
        zero = next((point for point in valid if point.ranks == (0, 0, 0)), None)
        interval = None
        truth = None
        truth_descriptive = None
        if truth_record is not None:
            truth_descriptive = ICPoint(
                TRUE_RANK,
                max(2.0 * float(truth_record["objective"]), np.finfo(float).tiny),
                int(truth_record["dimension"]),
            )
        if truth_valid:
            truth = next(point for point in valid if point.ranks == TRUE_RANK)
            interval = optimality_interval(
                truth, valid, base_penalty_per_dimension=base_per_dimension
            )
        record_by_rank = {tuple(row["ranks"]): row for row in records}
        baseline_choice = select_ic_candidate(
            valid, 1.0, base_penalty_per_dimension=base_per_dimension
        )
        baseline_reconstructed_ic = (
            baseline_choice.log_qhat
            + base_per_dimension * baseline_choice.dimension
        )
        interval_rows.append(
            {
                "dgp": int(saved["dgp"]),
                "N": n,
                "T": t,
                "replication": int(saved["replication"]),
                "analysis_status": "eligible" if truth_valid else "true_rank_candidate_invalid",
                "candidate_count_saved": len(records),
                "valid_candidate_count": len(valid),
                "candidate_data_complete": True,
                "final_candidate_set_persisted": True,
                "true_candidate_valid": truth_valid,
                "true_candidate_invalid_reasons": "|".join(
                    [] if truth_record is None else truth_record["invalid_reasons"]
                ),
                "true_candidate_sources": "|".join(
                    [] if truth_record is None else truth_record["sources"]
                ),
                "Qhat_true": np.nan if truth_record is None else max(2.0 * float(truth_record["objective"]), np.finfo(float).tiny),
                "log_Qhat_true": np.nan
                if truth_record is None
                else np.log(max(2.0 * float(truth_record["objective"]), np.finfo(float).tiny)),
                "d_true": int(truth_record["dimension"])
                if truth_record is not None
                else int(sum(n + t - 1 for _ in TRUE_RANK)),
                "P_true": dimension_penalty(TRUE_RANK, n, t),
                "Qhat_zero": np.nan if zero is None else zero.qhat,
                "log_Qhat_zero": np.nan if zero is None else zero.log_qhat,
                "log_loss_improvement_truth_over_zero": np.nan
                if zero is None or truth_record is None
                else zero.log_qhat
                - np.log(max(2.0 * float(truth_record["objective"]), np.finfo(float).tiny)),
                "d_zero": np.nan if zero is None else zero.dimension,
                "P_zero": dimension_penalty((0, 0, 0), n, t),
                "P_true_minus_zero": dimension_penalty(TRUE_RANK, n, t)
                - dimension_penalty((0, 0, 0), n, t),
                "c_star_zero_vs_true": np.nan
                if not truth_valid or zero is None
                else pairwise_break_even(
                    truth, zero, base_penalty_per_dimension=base_per_dimension
                ),
                "c_lower": np.nan if interval is None else interval.lower,
                "c_upper": np.nan if interval is None else interval.upper,
                "interval_empty": pd.NA if interval is None else interval.empty,
                "lower_binding_rank": "" if interval is None else _rank_text(interval.lower_binding_rank),
                "upper_binding_rank": "" if interval is None else _rank_text(interval.upper_binding_rank),
                "lower_binding_sources": ""
                if interval is None or interval.lower_binding_rank is None
                else "|".join(record_by_rank[interval.lower_binding_rank]["sources"]),
                "upper_binding_sources": ""
                if interval is None or interval.upper_binding_rank is None
                else "|".join(record_by_rank[interval.upper_binding_rank]["sources"]),
                "baseline_selected_rank": _rank_text(json.loads(saved["selected_rank_vector"])),
                "baseline_selected_ic": float(saved["selected_ic"]),
                "baseline_true_ic": float(saved["true_rank_ic"]),
                "baseline_reconstructed_rank": _rank_text(baseline_choice.ranks),
                "baseline_reconstructed_ic": baseline_reconstructed_ic,
                "baseline_reconstruction_pass": (
                    _rank_text(baseline_choice.ranks)
                    == _rank_text(json.loads(saved["selected_rank_vector"]))
                    and np.isclose(
                        baseline_reconstructed_ic,
                        float(saved["selected_ic"]),
                        rtol=1e-12,
                        atol=1e-12,
                    )
                ),
            }
        )
        for multiplier in GRID:
            choice = select_ic_candidate(
                valid, multiplier, base_penalty_per_dimension=base_per_dimension
            )
            selected_ic = choice.log_qhat + multiplier * base_per_dimension * choice.dimension
            true_ic = (
                np.nan
                if truth_descriptive is None
                else truth_descriptive.log_qhat
                + multiplier * base_per_dimension * truth_descriptive.dimension
            )
            grid_rows.append(
                {
                    "dgp": int(saved["dgp"]),
                    "N": n,
                    "T": t,
                    "replication": int(saved["replication"]),
                    "c_kappa": multiplier,
                    "selected_rank": _rank_text(choice.ranks),
                    "selected_ic": selected_ic,
                    "true_ic": true_ic,
                    "true_candidate_valid": truth_valid,
                    "true_rank_selected": truth_valid and choice.ranks == TRUE_RANK,
                    "underfit": any(rank < target for rank, target in zip(choice.ranks, TRUE_RANK, strict=True)),
                    "overfit": any(rank > target for rank, target in zip(choice.ranks, TRUE_RANK, strict=True)),
                }
            )
    return pd.DataFrame(interval_rows), pd.DataFrame(grid_rows)


def _penalty_table() -> pd.DataFrame:
    rows = []
    for n in (50, 100, 200, 400):
        nt = n * n
        log_nt = np.log(nt)
        b_nt = nt ** (1.0 / 12.0) * log_nt
        base_kappa = revision8_kappa(n, n)
        for rank in (0, 1, 2):
            delta_d = rank_increment_dimension(rank, n, n)
            row: dict[str, Any] = {
                "N": n,
                "T": n,
                "NT": nt,
                "rank_r": rank,
                "log_NT": log_nt,
                "b_NT": b_nt,
                "base_kappa": base_kappa,
                "delta_d": delta_d,
                "base_increment": base_kappa * delta_d / nt,
                "delta_d_1": 2 * n - 1,
                "base_increment_1": base_kappa * (2 * n - 1) / nt,
            }
            for multiplier in GRID:
                row[f"increment_c_{multiplier:.0e}"] = multiplier * row["base_increment"]
                row[f"increment_1_c_{multiplier:.0e}"] = (
                    multiplier * row["base_increment_1"]
                )
            rows.append(row)
    return pd.DataFrame(rows)


def _latex_escape(value: Any) -> str:
    if pd.isna(value):
        return ""
    if isinstance(value, (float, np.floating)):
        if value == np.inf:
            return r"$\infty$"
        if value == -np.inf:
            return r"$-\infty$"
        return f"{value:.6g}"
    text = str(value)
    return (
        text.replace("\\", r"\textbackslash{}")
        .replace("_", r"\_")
        .replace("%", r"\%")
        .replace("&", r"\&")
    )


def _tabular(
    rows: list[list[Any]], headers: list[str], caption: str, label: str, notes: str
) -> str:
    alignment = "l" + "r" * (len(headers) - 1)
    lines = [
        r"\begin{table}[!htbp]",
        r"\centering",
        rf"\caption{{{caption}}}",
        rf"\label{{{label}}}",
        rf"\begin{{tabular}}{{{alignment}}}",
        r"\toprule",
        " & ".join(headers) + r" \\",
        r"\midrule",
    ]
    lines.extend(" & ".join(_latex_escape(value) for value in row) + r" \\" for row in rows)
    lines += [
        r"\bottomrule",
        r"\end{tabular}",
        rf"\begin{{minipage}}{{0.96\linewidth}}\footnotesize \textit{{Notes:}} {notes}\end{{minipage}}",
        r"\end{table}",
        "",
    ]
    return "\n".join(lines)


def _write_interval_latex(frame: pd.DataFrame, path: Path) -> None:
    rows = [
        [
            row.dgp,
            row.replication,
            "yes" if row.true_candidate_valid else "no",
            row.Qhat_true,
            row.Qhat_zero,
            row.c_star_zero_vs_true,
            row.c_lower,
            row.c_upper,
            row.lower_binding_rank,
            row.upper_binding_rank,
        ]
        for row in frame.itertuples()
    ]
    path.write_text(
        _tabular(
            rows,
            ["DGP", "Rep.", "Valid truth", r"$\widehat Q_1$", r"$\widehat Q_0$", "$c^*_{0,1}$", "$c_L$", "$c_U$", "Lower binder", "Upper binder"],
            "Preflight IC break-even intervals",
            "tab:ic-break-even",
            "Intervals intersect comparisons against every valid saved candidate. An empty entry means the true-rank post-refit was numerically invalid.",
        ),
        encoding="utf-8",
    )


def _write_grid_latex(frame: pd.DataFrame, path: Path) -> None:
    summary = frame.groupby("c_kappa", sort=True).agg(
        true_rank=("true_rank_selected", "sum"),
        underfit=("underfit", "sum"),
        overfit=("overfit", "sum"),
    )
    by_dgp = frame.pivot_table(
        index="c_kappa", columns="dgp", values="true_rank_selected", aggfunc="sum"
    )
    rows = []
    for multiplier, row in summary.iterrows():
        rows.append(
            [multiplier, int(row.true_rank), int(row.underfit), int(row.overfit)]
            + [int(by_dgp.loc[multiplier].get(dgp, 0)) for dgp in (1, 2, 3, 4)]
        )
    path.write_text(
        _tabular(
            rows,
            [r"$c_\kappa$", "True", "Under", "Over", "DGP 1", "DGP 2", "DGP 3", "DGP 4"],
            "Preflight IC multiplier grid",
            "tab:ic-grid",
            "True-rank selections are shown overall and by DGP. The numerically invalid DGP 3, replication 0 true-rank fit cannot count as selected. Under- and overfit are separate componentwise indicators across all eight replications.",
        ),
        encoding="utf-8",
    )


def _write_penalty_latex(frame: pd.DataFrame, path: Path) -> None:
    parts = []
    for rank in (0, 1, 2):
        subset = frame.loc[frame["rank_r"].eq(rank)]
        rows = []
        for multiplier in GRID:
            rows.append(
                [multiplier]
                + [
                    subset.loc[subset["N"].eq(n), f"increment_c_{multiplier:.0e}"].iloc[0]
                    for n in (50, 100, 200, 400)
                ]
            )
        parts.append(
            _tabular(
                rows,
                [r"$c_\kappa$", "$N=50$", "$N=100$", "$N=200$", "$N=400$"],
                f"Revision-8 penalty increments: rank {rank} to {rank + 1}",
                f"tab:ic-penalty-r{rank}",
                rf"Each entry is $c_\kappa\kappa_{{NT}}\Delta d/(NT)$ for one matrix, with $\Delta d=N+T-2({rank})-1$, $\eta=4$, and $d_s=1$.",
            )
        )
    path.write_text("\n".join(parts), encoding="utf-8")


def _report(intervals: pd.DataFrame, grid: pd.DataFrame, penalty: pd.DataFrame) -> str:
    eligible = intervals.loc[intervals["true_candidate_valid"]]
    finite_upper = eligible.loc[np.isfinite(eligible["c_upper"]), "c_upper"]
    nonempty = eligible.loc[~eligible["interval_empty"].astype(bool)]
    zero_breaks = eligible["c_star_zero_vs_true"].dropna()
    common_lower = eligible["c_lower"].max()
    common_upper = eligible["c_upper"].min()
    common_nonempty = common_lower <= common_upper
    lower_row = eligible.loc[eligible["c_lower"].idxmax()]
    upper_row = eligible.loc[eligible["c_upper"].idxmin()]
    endpoints = sorted(set(eligible["c_lower"].tolist() + eligible["c_upper"].tolist()))
    probes = endpoints + [
        (left + right) / 2
        for left, right in zip(endpoints[:-1], endpoints[1:], strict=True)
        if np.isfinite(right)
    ]
    max_overlap = max(
        sum(row.c_lower <= value <= row.c_upper for row in eligible.itertuples())
        for value in probes
    )
    grid_summary = grid.groupby("c_kappa", sort=True).agg(
        exact_true=("true_rank_selected", "sum"),
        underfit=("underfit", "sum"),
        overfit=("overfit", "sum"),
    )
    recommended = grid_summary.sort_values(
        ["exact_true", "underfit", "overfit"], ascending=[False, True, True]
    ).head(1)
    first_increment = penalty.loc[penalty["rank_r"].eq(0)].set_index("N")
    lines = [
        "# Offline IC break-even diagnostic",
        "",
        "This report reconstructs the Revision-8 IC from saved preflight candidate fits. It performs no fitting, simulation, or estimator changes.",
        "",
        "## Saved-candidate audit",
        "",
        f"- Successful replications: {len(intervals)}.",
        "- Every persisted candidate record contains rank, objective, normalized residual variance, model dimension, validity, source attribution, and baseline IC.",
        "- Candidate-record counts equal the saved final locally completed candidate counts in all eight replications.",
        f"- Replications with a valid true-rank post-refit: {len(eligible)}.",
        f"- Replications with an invalid true-rank post-refit: {len(intervals) - len(eligible)}.",
        "- DGP 3, replication 0 has the true rank in the screened candidate set, but it is excluded from IC competition because objective stability failed.",
        "",
        "## Exact interval summary",
        "",
        f"- Nonempty true-rank intervals among eligible replications: {len(nonempty)}/{len(eligible)}.",
        f"- Median lower endpoint: {eligible['c_lower'].median():.6g}.",
        f"- Median finite upper endpoint: {finite_upper.median():.6g}." if len(finite_upper) else "- No finite upper endpoints.",
        f"- Smallest finite upper endpoint: {finite_upper.min():.6g}." if len(finite_upper) else "- No finite upper endpoints.",
        f"- Maximum lower endpoint: {common_lower:.6g} (DGP {int(lower_row.dgp)}, replication {int(lower_row.replication)}).",
        f"- Minimum upper endpoint: {common_upper:.6g} (DGP {int(upper_row.dgp)}, replication {int(upper_row.replication)}).",
        f"- Common intersection: [{common_lower:.6g}, {common_upper:.6g}]."
        if common_nonempty
        else f"- Common intersection: empty; at most {max_overlap}/{len(eligible)} eligible intervals overlap. DGP {int(lower_row.dgp)}, replication {int(lower_row.replication)} sets the largest lower bound, while DGP {int(upper_row.dgp)}, replication {int(upper_row.replication)} sets the smallest upper bound.",
        "",
        "The lower endpoint is imposed by higher-dimensional competitors; the upper endpoint is imposed by lower-dimensional competitors. Bounds are the intersection of all valid saved candidates, not only the all-zero model.",
        "",
        "## Diagnostic grid",
        "",
        "| c_kappa | true rank | underfit | overfit | selected-rank distribution |",
        "|---:|---:|---:|---:|:---|",
        *[
            f"| {multiplier:.0e} | {int(row.exact_true)}/7 ({row.exact_true / 7:.1%}) | {int(row.underfit)}/8 ({row.underfit / 8:.1%}) | {int(row.overfit)}/8 ({row.overfit / 8:.1%}) | {', '.join(f'{rank}: {count}' for rank, count in grid.loc[grid['c_kappa'].eq(multiplier), 'selected_rank'].value_counts().sort_index().items())} |"
            for multiplier, row in grid_summary.iterrows()
        ],
        "",
        "The true-rank count uses only numerically valid true-rank post-refits; its maximum is therefore seven, not eight. Underfit and overfit are separate componentwise indicators, so a mixed rank vector can satisfy both.",
        "",
        "## Shortlist for later testing",
        "",
    ]
    for multiplier, row in recommended.iterrows():
        lines.append(
            f"- `c_kappa = {multiplier:.0e}`: true rank in {int(row['exact_true'])}/7 eligible replications; underfit in {int(row['underfit'])}/8 and overfit in {int(row['overfit'])}/8 saved replications."
        )
    lines += [
        "",
        "These are diagnostic shortlist values only. The active configuration remains unchanged; a later lower-cap preflight must be approved and executed before choosing a production constant.",
        "",
        "## Literal Revision-8 penalty",
        "",
        "The base rate is `kappa_NT = b_NT^2 log(NT)^(d_s+3)`, with `b_NT = (NT)^(1/(8+eta)) log(NT)`, `eta=4`, and `d_s=1`. The penalty table reports increments for raising one matrix from rank r to r+1 at N=T in {50,100,200,400}.",
        "",
        "| N=T | NT | log(NT) | b_NT | base kappa | delta d1 | base increment 1 |",
        "|---:|---:|---:|---:|---:|---:|---:|",
        *[
            f"| {n} | {int(row.NT)} | {row.log_NT:.6g} | {row.b_NT:.6g} | {row.base_kappa:.6g} | {int(row.delta_d_1)} | {row.base_increment_1:.6g} |"
            for n, row in first_increment.iterrows()
        ],
        "",
        "First-rank IC increments under each diagnostic multiplier are tabulated in `tab_ic_penalty_magnitude.csv`; the LaTeX artifact additionally provides separate panels for rank increases 0-to-1, 1-to-2, and 2-to-3.",
        "",
        "## Zero-versus-truth break-even distribution",
        "",
        f"Among the {len(zero_breaks)} replications with valid true-rank candidates: minimum {zero_breaks.min():.6g}; 10th percentile {zero_breaks.quantile(0.1):.6g}; median {zero_breaks.median():.6g}; 90th percentile {zero_breaks.quantile(0.9):.6g}; maximum {zero_breaks.max():.6g}.",
    ]
    return "\n".join(lines) + "\n"


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--input",
        type=Path,
        default=Path("results/mc/cap_pilot_preflight/538330907b21af7f/rank"),
    )
    parser.add_argument(
        "--output", type=Path, default=Path("results/mc/diagnostics/ic_break_even")
    )
    args = parser.parse_args()
    args.output.mkdir(parents=True, exist_ok=True)
    intervals, grid = _reconstruct(_load_records(args.input))
    penalty = _penalty_table()

    intervals.to_csv(args.output / "tab_ic_break_even_preflight.csv", index=False)
    intervals.to_parquet(args.output / "tab_ic_break_even_preflight.parquet", index=False)
    grid.to_csv(args.output / "tab_ic_grid_preflight.csv", index=False)
    grid.to_parquet(args.output / "tab_ic_grid_preflight.parquet", index=False)
    penalty.to_csv(args.output / "tab_ic_penalty_magnitude.csv", index=False)
    _write_interval_latex(intervals, args.output / "tab_ic_break_even_preflight.tex")
    _write_grid_latex(grid, args.output / "tab_ic_grid_preflight.tex")
    _write_penalty_latex(penalty, args.output / "tab_ic_penalty_magnitude.tex")
    (args.output / "ic_break_even_report.md").write_text(
        _report(intervals, grid, penalty), encoding="utf-8"
    )


if __name__ == "__main__":
    main()
