"""Aggregation and journal tables for the frozen medium diagnostics."""

from __future__ import annotations

import json
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd

from .reporting import TARGET_TITLES

TRUE_BASELINE = "[1, 1, 1]"
UNSUPPORTED = {
    "target_unsupported_selected_rank",
    "split_target_unsupported_selected_rank",
}
GRAM_FAILURES = {
    "tangent_gram_eigensolver_failure",
    "tangent_gram_nearly_singular",
    "split_tangent_gram_eigensolver_failure",
    "split_tangent_gram_nearly_singular",
}


def _planned_rank(config: dict[str, Any], true_ranks: list[str]) -> pd.DataFrame:
    return pd.DataFrame(
        [
            {
                "dgp": int(dgp),
                "N": int(n),
                "T": int(t),
                "replication": replication,
                "true_rank_vector": true_rank,
            }
            for dgp in config["run"]["dgps"]
            for n, t in config["run"]["cells"]
            for true_rank in true_ranks
            for replication in range(int(config["run"]["replications"]))
        ]
    )


def _failure_status(raw_rows: pd.DataFrame) -> pd.DataFrame:
    if raw_rows.empty or "record_type" not in raw_rows:
        return pd.DataFrame(
            columns=["dgp", "N", "T", "replication", "true_rank_vector", "failure_status"]
        )
    failures = raw_rows.loc[raw_rows["record_type"] == "failure"].copy()
    if failures.empty:
        return pd.DataFrame(
            columns=["dgp", "N", "T", "replication", "true_rank_vector", "failure_status"]
        )
    keys = ["dgp", "N", "T", "replication", "true_rank_vector"]
    return failures.groupby(keys, as_index=False)["status"].first().rename(
        columns={"status": "failure_status"}
    )


def complete_rank_replications(
    rank_rows: pd.DataFrame,
    raw_rows: pd.DataFrame,
    config: dict[str, Any],
    true_ranks: list[str],
) -> pd.DataFrame:
    """Complete the requested rank grid and retain replication failure reasons."""

    keys = ["dgp", "N", "T", "replication", "true_rank_vector"]
    complete = _planned_rank(config, true_ranks).merge(
        rank_rows, how="left", on=keys, suffixes=("", "_rank")
    )
    complete = complete.merge(_failure_status(raw_rows), how="left", on=keys)
    complete["numerical_instability"] = complete["selected_rank_vector"].isna()
    complete["candidate_covered"] = complete["true_rank_in_candidates"].fillna(False).astype(bool)
    complete["exact"] = complete["exact_rank_recovery"].fillna(False).astype(bool)
    complete["selected_all_zero"] = complete["selected_rank_vector"].eq("[0, 0, 0]")
    return complete


def aggregate_rank_medium(replications: pd.DataFrame) -> pd.DataFrame:
    """Aggregate baseline rank recovery and its failure decomposition."""

    rows: list[dict[str, Any]] = []
    group_keys = ["dgp", "N", "T", "true_rank_vector"]
    for keys, group in replications.groupby(group_keys, sort=True):
        available = group.loc[~group["numerical_instability"]]
        covered = group.loc[group["candidate_covered"]]
        failure_distribution = (
            group["failure_status"].dropna().value_counts().sort_index().to_dict()
        )
        rows.append(
            {
                **dict(zip(group_keys, keys, strict=True)),
                "requested_replications": len(group),
                "rank_records": len(available),
                "candidate_coverage": float(group["candidate_covered"].mean()),
                "candidate_absence_rate": float((~group["candidate_covered"]).mean()),
                "exact_rank_recovery": float(group["exact"].mean()),
                "ic_selection_error_conditional_coverage": (
                    float((~covered["exact"]).mean()) if len(covered) else np.nan
                ),
                "A_underselection": float(available["A_underselected"].mean()),
                "A_overselection": float(available["A_overselected"].mean()),
                "B_underselection": float(available["B_underselected"].mean()),
                "B_overselection": float(available["B_overselected"].mean()),
                "H_underselection": float(available["H_underselected"].mean()),
                "H_overselection": float(available["H_overselected"].mean()),
                "rank_zero_recovery": (
                    float(available["zero_rank_recovery"].mean())
                    if " 0" in keys[3]
                    else np.nan
                ),
                "selected_all_zero_rate": float(group["selected_all_zero"].mean()),
                "mean_candidate_count": float(available["candidate_count_final"].mean()),
                "mean_true_rank_ic_when_present": float(available["true_rank_ic"].mean()),
                "mean_selected_rank_ic": float(available["selected_ic"].mean()),
                "mean_selected_minus_true_ic": float(
                    available["selected_minus_true_ic"].mean()
                ),
                "candidate_coverage_failure_rate": float(
                    ((~group["candidate_covered"]) & ~group["numerical_instability"]).mean()
                ),
                "ic_choice_failure_rate": float(
                    (group["candidate_covered"] & ~group["exact"]).mean()
                ),
                "numerical_instability_rate": float(group["numerical_instability"].mean()),
                "cap_hit_rate": float(available["rank_at_cap"].fillna(False).mean()),
                "rank_cap_thresholded_distribution": json.dumps(
                    available["rank_cap_thresholded_vector"]
                    .value_counts(normalize=True)
                    .sort_index()
                    .to_dict(),
                    sort_keys=True,
                ),
                "selected_rank_distribution": json.dumps(
                    available["selected_rank_vector"]
                    .value_counts(normalize=True)
                    .sort_index()
                    .to_dict(),
                    sort_keys=True,
                ),
                "failure_status_distribution": json.dumps(
                    failure_distribution, sort_keys=True
                ),
            }
        )
    return pd.DataFrame(rows)


def aggregate_rank_sensitivity_medium(replications: pd.DataFrame) -> pd.DataFrame:
    """Aggregate candidate coverage and exact recovery for frozen multipliers."""

    specifications = [("baseline", "selected_rank_vector", "true_rank_in_candidates")]
    for family in ("ic", "threshold"):
        for multiplier in ("0.5", "1.0", "2.0"):
            specifications.append(
                (
                    f"{family}_multiplier_{multiplier}",
                    f"{family}_multiplier_{multiplier}_selected_rank",
                    f"{family}_multiplier_{multiplier}_true_rank_in_candidates",
                )
            )
    rows: list[dict[str, Any]] = []
    for keys, group in replications.groupby(["dgp", "N", "T"], sort=True):
        for specification, selected_column, coverage_column in specifications:
            selected = (
                group[selected_column]
                if selected_column in group
                else pd.Series(index=group.index, dtype=object)
            )
            coverage = (
                group[coverage_column].fillna(False).astype(bool)
                if coverage_column in group
                else pd.Series(False, index=group.index, dtype=bool)
            )
            exact = selected.eq(TRUE_BASELINE)
            covered = coverage & selected.notna()
            rows.append(
                {
                    "dgp": keys[0],
                    "N": keys[1],
                    "T": keys[2],
                    "specification": specification,
                    "requested_replications": len(group),
                    "available_replications": int(selected.notna().sum()),
                    "candidate_coverage": float(coverage.mean()),
                    "exact_rank_recovery": float(exact.fillna(False).mean()),
                    "ic_selection_error_conditional_coverage": (
                        float((~exact[covered]).mean()) if covered.any() else np.nan
                    ),
                    "selected_rank_distribution": json.dumps(
                        selected.dropna().value_counts(normalize=True).sort_index().to_dict(),
                        sort_keys=True,
                    ),
                }
            )
    return pd.DataFrame(rows)


def _headline_target(name: str, dgp: int) -> bool:
    if name == "B_entry":
        return False
    return not ("G2_minus_G1_fixed_time" in name and dgp < 4)


def _planned_targets(config: dict[str, Any]) -> pd.DataFrame:
    targets = config["inference"]["targets"]
    return pd.DataFrame(
        [
            {
                "dgp": int(dgp),
                "N": int(n),
                "T": int(t),
                "replication": replication,
                "target": target,
            }
            for dgp in config["run"]["dgps"]
            for n, t in config["run"]["cells"]
            for replication in range(int(config["run"]["replications"]))
            for target in targets
            if _headline_target(target, int(dgp))
        ]
    )


def complete_target_replications(
    raw_rows: pd.DataFrame, config: dict[str, Any]
) -> pd.DataFrame:
    """Complete the theorem-covered target grid."""

    targets = raw_rows.loc[raw_rows["record_type"] == "target"].copy()
    keys = ["dgp", "N", "T", "replication", "target"]
    complete = _planned_targets(config).merge(targets, how="left", on=keys)
    complete["target_record"] = complete["status"].notna()
    complete["full_target_supported"] = (
        complete["target_record"]
        & complete["riesz_target_tangent_norm"].fillna(0.0).gt(
            float(config["inference"]["target_support_tolerance"])
        )
    )
    complete["primary_target_supported"] = (
        complete["full_target_supported"] & ~complete["status"].isin(UNSUPPORTED)
    )
    return complete


def _split_solver_failure(value: Any) -> bool:
    if not isinstance(value, str) or not value:
        return False
    return any(
        item.get("target_supported", False) and not item.get("riesz_converged", True)
        for item in json.loads(value)
    )


def _split_target_instability(value: Any) -> bool:
    if not isinstance(value, str) or not value:
        return False
    return any(
        item.get("target_supported", False)
        and item.get("riesz_converged", False)
        and not item.get("riesz_target_stable", True)
        for item in json.loads(value)
    )


def aggregate_target_support_medium(targets: pd.DataFrame) -> pd.DataFrame:
    """Aggregate support and conditional numerical-failure rates."""

    rows: list[dict[str, Any]] = []
    group_keys = ["dgp", "N", "T", "target"]
    for keys, group in targets.groupby(group_keys, sort=True):
        supported = group["primary_target_supported"]
        denominator = int(supported.sum())
        split_solver = group["split_diagnostics_json"].map(_split_solver_failure)
        split_instability = group["split_diagnostics_json"].map(_split_target_instability)
        solver_failure = group["status"].eq("riesz_not_converged") | split_solver
        target_instability = group["status"].eq("riesz_target_instability") | split_instability
        gram_failure = group["status"].isin(GRAM_FAILURES)
        rows.append(
            {
                **dict(zip(group_keys, keys, strict=True)),
                "requested_replications": len(group),
                "target_records": int(group["target_record"].sum()),
                "selected_rank_target_support_rate": float(supported.mean()),
                "full_target_support_rate": float(group["full_target_supported"].mean()),
                "target_unsupported_selected_rank_rate": float(
                    group["status"].isin(UNSUPPORTED).mean()
                ),
                "riesz_solver_failure_rate_conditional_support": (
                    float(solver_failure[supported].mean()) if denominator else np.nan
                ),
                "riesz_target_instability_rate_conditional_support": (
                    float(target_instability[supported].mean()) if denominator else np.nan
                ),
                "tangent_gram_failure_rate_conditional_support": (
                    float(gram_failure[supported].mean()) if denominator else np.nan
                ),
                "successful_inference_rate": float(group["status"].eq("success").mean()),
                "mean_target_tangent_norm": float(group["riesz_target_tangent_norm"].mean()),
                "minimum_target_tangent_norm": float(group["riesz_target_tangent_norm"].min()),
                "mean_true_projection_ratio": float(group["true_target_projection_ratio"].mean()),
                "minimum_true_projection_ratio": float(group["true_target_projection_ratio"].min()),
                "mean_true_entry_unit_leverage_scaled": float(
                    group["true_entry_unit_leverage_scaled"].mean()
                ),
                "minimum_true_entry_unit_leverage_scaled": float(
                    group["true_entry_unit_leverage_scaled"].min()
                ),
                "mean_true_entry_time_leverage_scaled": float(
                    group["true_entry_time_leverage_scaled"].mean()
                ),
                "minimum_true_entry_time_leverage_scaled": float(
                    group["true_entry_time_leverage_scaled"].min()
                ),
            }
        )
    return pd.DataFrame(rows)


def _system_rows(targets: pd.DataFrame) -> pd.DataFrame:
    rows: list[dict[str, Any]] = []
    for row in targets.loc[targets["target_record"]].itertuples(index=False):
        base = {"dgp": row.dgp, "N": row.N, "T": row.T, "target": row.target}
        rows.append(
            {
                **base,
                "system": "full",
                "target_supported": row.full_target_supported,
                "solver_called": row.status not in GRAM_FAILURES
                and row.status not in UNSUPPORTED
                and row.full_target_supported,
                "solver_failure": row.status == "riesz_not_converged",
                "rayleigh": row.riesz_target_rayleigh_quotient,
                "gram_min": row.tangent_gram_smallest_eigenvalue,
                "gram_condition": row.tangent_gram_condition_number,
            }
        )
        split_text = getattr(row, "split_diagnostics_json", None)
        if isinstance(split_text, str) and split_text:
            for split in json.loads(split_text):
                rows.append(
                    {
                        **base,
                        "system": f"{split['kind']}_{split['part']}",
                        "target_supported": split.get("target_supported", False),
                        "solver_called": "riesz_converged" in split,
                        "solver_failure": not split.get("riesz_converged", True),
                        "rayleigh": split.get("riesz_target_rayleigh_quotient", np.nan),
                        "gram_min": split.get("tangent_gram_smallest_eigenvalue", np.nan),
                        "gram_condition": split.get("tangent_gram_condition_number", np.nan),
                    }
                )
    return pd.DataFrame(rows)


def aggregate_riesz_medium(targets: pd.DataFrame) -> pd.DataFrame:
    """Summarize full and split Riesz/Gram distributions on supported systems."""

    systems = _system_rows(targets)
    rows: list[dict[str, Any]] = []
    for keys, group in systems.groupby(["dgp", "N", "T", "target"], sort=True):
        supported = group.loc[group["target_supported"]]
        solver_attempts = supported.loc[supported["solver_called"]]
        rows.append(
            {
                "dgp": keys[0],
                "N": keys[1],
                "T": keys[2],
                "target": keys[3],
                "supported_systems": len(supported),
                "riesz_solver_failure_rate": (
                    float(solver_attempts["solver_failure"].mean())
                    if len(solver_attempts)
                    else np.nan
                ),
                "minimum_tangent_gram_eigenvalue": float(supported["gram_min"].min()),
                "p10_tangent_gram_eigenvalue": float(supported["gram_min"].quantile(0.10)),
                "median_tangent_gram_eigenvalue": float(supported["gram_min"].median()),
                "median_tangent_gram_condition_number": float(
                    supported["gram_condition"].median()
                ),
                "p90_tangent_gram_condition_number": float(
                    supported["gram_condition"].quantile(0.90)
                ),
                "maximum_tangent_gram_condition_number": float(
                    supported["gram_condition"].max()
                ),
                "minimum_riesz_target_rayleigh_quotient": float(
                    supported["rayleigh"].min()
                ),
                "p10_riesz_target_rayleigh_quotient": float(
                    supported["rayleigh"].quantile(0.10)
                ),
                "median_riesz_target_rayleigh_quotient": float(
                    supported["rayleigh"].median()
                ),
            }
        )
    return pd.DataFrame(rows)


def attach_rank_stress_calibration(
    summary: pd.DataFrame, calibration: pd.DataFrame
) -> pd.DataFrame:
    """Attach achieved R2 and deterministic envelopes to rank-stress results."""

    calibration = calibration.rename(columns={"n": "N", "t": "T"})
    columns = [
        "dgp",
        "N",
        "T",
        "true_rank_vector",
        "target_r2",
        "requested_r2",
        "achieved_r2",
        "r2_scale_identified",
        "c_h",
        "c_xi",
        "theoretical_coefficient_envelope",
    ]
    result = summary.merge(calibration[columns], how="left", on=["dgp", "N", "T", "true_rank_vector"])
    has_zero = result["true_rank_vector"].str.contains(" 0")
    result["rank_zero_recovery"] = np.where(
        has_zero, result["rank_zero_recovery"], np.nan
    )
    return result


def _escape(value: Any) -> str:
    return str(value).replace("_", r"\_").replace("%", r"\%").replace("&", r"\&")


def _format(value: Any, digits: int = 3) -> str:
    if value is None or pd.isna(value):
        return "--"
    if isinstance(value, (int, np.integer)):
        return str(int(value))
    return f"{float(value):.{digits}f}"


def _write_table(
    frame: pd.DataFrame,
    output: Path,
    stem: str,
    *,
    caption: str,
    label: str,
    panel_column: str,
    columns: list[tuple[str, str]],
    notes: str,
) -> list[Path]:
    output.mkdir(parents=True, exist_ok=True)
    csv_path = output / f"{stem}.csv"
    parquet_path = output / f"{stem}.parquet"
    tex_path = output / f"{stem}.tex"
    frame.to_csv(csv_path, index=False)
    frame.to_parquet(parquet_path, index=False)
    newline = r"\\"
    alignment = "l" + "r" * len(columns)
    lines = [
        r"\begingroup",
        r"\small",
        rf"\begin{{longtable}}{{{alignment}}}",
        rf"\caption{{{caption}}}\label{{{label}}}" + newline,
        r"\toprule",
        "Panel & " + " & ".join(title for _, title in columns) + " " + newline,
        r"\midrule",
        r"\endfirsthead",
        r"\toprule",
        "Panel & " + " & ".join(title for _, title in columns) + " " + newline,
        r"\midrule",
        r"\endhead",
    ]
    panels = frame[panel_column].drop_duplicates().tolist()
    for index, panel in enumerate(panels):
        title = TARGET_TITLES.get(panel, str(panel))
        lines.append(
            rf"\multicolumn{{{len(columns) + 1}}}{{l}}{{\textit{{Panel {chr(65 + index)}. {_escape(title)}}}}} "
            + newline
        )
        for row in frame.loc[frame[panel_column] == panel].to_dict("records"):
            values = [_format(row.get(key)) for key, _ in columns]
            lines.append(" & " + " & ".join(values) + " " + newline)
        lines.append(r"\addlinespace")
    lines.extend(
        [
            r"\bottomrule",
            rf"\multicolumn{{{len(columns) + 1}}}{{p{{0.96\linewidth}}}}{{\footnotesize \textit{{Notes:}} {notes}}} "
            + newline,
            r"\end{longtable}",
            r"\endgroup",
            "",
        ]
    )
    tex_path.write_text("\n".join(lines), encoding="utf-8")
    return [csv_path, parquet_path, tex_path]


def write_medium_tables(
    rank: pd.DataFrame,
    support: pd.DataFrame,
    riesz: pd.DataFrame,
    sensitivity: pd.DataFrame,
    stress: pd.DataFrame,
    output: Path,
) -> list[Path]:
    """Write the five prespecified CSV/Parquet/LaTeX table families."""

    paths: list[Path] = []
    paths += _write_table(
        rank,
        output,
        "tab_mc_rank_medium",
        caption="Medium pre-production rank diagnostic",
        label="tab:mc-rank-medium",
        panel_column="dgp",
        columns=[
            ("N", "$N$"), ("T", "$T$"), ("candidate_coverage", "Coverage"),
            ("exact_rank_recovery", "Exact"),
            ("ic_selection_error_conditional_coverage", "IC error$|$covered"),
            ("A_underselection", "$A<$"), ("B_underselection", "$B<$"),
            ("H_underselection", "$H<$"), ("selected_all_zero_rate", "All zero"),
            ("mean_candidate_count", "Candidates"),
            ("numerical_instability_rate", "Numerical"), ("cap_hit_rate", "Cap hit"),
        ],
        notes="Coverage absence, IC error conditional on coverage, numerical instability, and cap hits are reported separately.",
    )
    paths += _write_table(
        support,
        output,
        "tab_mc_target_support",
        caption="Medium target-support diagnostic",
        label="tab:mc-target-support",
        panel_column="target",
        columns=[
            ("dgp", "DGP"), ("N", "$N$"), ("T", "$T$"),
            ("requested_replications", "Requested"),
            ("selected_rank_target_support_rate", "Supported"),
            ("target_unsupported_selected_rank_rate", "Unsupported"),
            ("riesz_solver_failure_rate_conditional_support", "Solver fail$|$support"),
            ("riesz_target_instability_rate_conditional_support", "Target fail$|$support"),
            ("tangent_gram_failure_rate_conditional_support", "Gram fail$|$support"),
            ("successful_inference_rate", "Success"),
            ("minimum_true_projection_ratio", r"Min $\|P_0D\|/\|D\|$"),
        ],
        notes="Only theorem-covered targets enter this table. B-entry and DGP 1--3 fixed-time group contrasts are excluded from headline rates.",
    )
    paths += _write_table(
        riesz,
        output,
        "tab_mc_riesz_medium",
        caption="Medium Riesz and tangent-Gram diagnostic",
        label="tab:mc-riesz-medium",
        panel_column="target",
        columns=[
            ("dgp", "DGP"), ("N", "$N$"), ("T", "$T$"),
            ("supported_systems", "Systems"),
            ("riesz_solver_failure_rate", "Solver fail"),
            ("minimum_tangent_gram_eigenvalue", "Min Gram eig."),
            ("p10_tangent_gram_eigenvalue", "P10 Gram eig."),
            ("median_tangent_gram_condition_number", "Med. cond."),
            ("p90_tangent_gram_condition_number", "P90 cond."),
            ("minimum_riesz_target_rayleigh_quotient", "Min target RQ"),
            ("median_riesz_target_rayleigh_quotient", "Med. target RQ"),
        ],
        notes="Gram eigenvalues use cached nonredundant tangent coordinates; target Rayleigh quotients remain target-specific diagnostics.",
    )
    paths += _write_table(
        sensitivity,
        output,
        "tab_mc_rank_sensitivity_medium",
        caption="Medium rank-selection sensitivity diagnostic",
        label="tab:mc-rank-sensitivity-medium",
        panel_column="specification",
        columns=[
            ("dgp", "DGP"), ("N", "$N$"), ("T", "$T$"),
            ("candidate_coverage", "Coverage"),
            ("exact_rank_recovery", "Exact"),
            ("ic_selection_error_conditional_coverage", "Error$|$covered"),
            ("available_replications", "Available"),
        ],
        notes="Each row applies the same rank-selection algorithm with only the named multiplier changed; no multiplier is selected automatically.",
    )
    paths += _write_table(
        stress,
        output,
        "tab_mc_rank_stress_medium",
        caption="Medium rank-stress diagnostic",
        label="tab:mc-rank-stress-medium",
        panel_column="true_rank_vector",
        columns=[
            ("dgp", "DGP"), ("N", "$N$"), ("T", "$T$"),
            ("candidate_coverage", "Coverage"), ("exact_rank_recovery", "Exact"),
            ("A_underselection", "$A<$"), ("A_overselection", "$A>$"),
            ("B_underselection", "$B<$"), ("B_overselection", "$B>$"),
            ("H_underselection", "$H<$"), ("H_overselection", "$H>$"),
            ("rank_zero_recovery", "Zero recovery"),
            ("achieved_r2", "$R^2$"),
            ("theoretical_coefficient_envelope", "Envelope"),
        ],
        notes="For true rank (1,0,2), c_xi=1 is a normalization and the displayed pooled R-squared is induced, not a failed calibration target.",
    )
    return paths
