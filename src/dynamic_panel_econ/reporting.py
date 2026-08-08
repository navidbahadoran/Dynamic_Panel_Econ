"""Journal-style aggregation and deterministic panel-oriented table output."""

from __future__ import annotations

import json
import string
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd

from .monte_carlo import FAILURE_CODES

TARGET_TITLES = {
    "A_entry": "Autoregressive coefficient, individual entry",
    "B_entry": "Slope coefficient, individual entry",
    "A_fixed_time_mean": "Autoregressive coefficient, fixed-time mean",
    "B_fixed_time_mean": "Slope coefficient, fixed-time mean",
    "A_full_mean": "Autoregressive coefficient, full-panel mean---corrected",
    "B_full_mean": "Slope coefficient, full-panel mean---corrected",
    "A_G1_fixed_time": "Autoregressive coefficient, fixed-time group 1 mean",
    "A_G2_fixed_time": "Autoregressive coefficient, fixed-time group 2 mean",
    "A_G2_minus_G1_fixed_time": "Autoregressive coefficient, fixed-time group contrast",
    "B_G1_fixed_time": "Slope coefficient, fixed-time group 1 mean",
    "B_G2_fixed_time": "Slope coefficient, fixed-time group 2 mean",
    "B_G2_minus_G1_fixed_time": "Slope coefficient, fixed-time group contrast",
    "A_G1_time_average": "Autoregressive coefficient, time-averaged group 1 mean---corrected",
    "A_G2_time_average": "Autoregressive coefficient, time-averaged group 2 mean---corrected",
    "A_G2_minus_G1_time_average": "Autoregressive coefficient, time-averaged group contrast---corrected",
    "B_G1_time_average": "Slope coefficient, time-averaged group 1 mean---corrected",
    "B_G2_time_average": "Slope coefficient, time-averaged group 2 mean---corrected",
    "B_G2_minus_G1_time_average": "Slope coefficient, time-averaged group contrast---corrected",
}


def _read_chunks(directory: Path) -> pd.DataFrame:
    files = sorted(directory.glob("*.parquet"))
    return pd.concat((pd.read_parquet(path) for path in files), ignore_index=True) if files else pd.DataFrame()


def _manifest(root: Path) -> dict[str, Any]:
    return json.loads((root / "run_manifest.json").read_text(encoding="utf-8"))


def aggregate_run(root: str | Path) -> pd.DataFrame:
    """Aggregate primary target performance and a separate true-rank-conditioned summary."""

    root = Path(root)
    raw = _read_chunks(root / "raw")
    rank_raw = _read_chunks(root / "rank")
    if raw.empty and rank_raw.empty:
        raise FileNotFoundError(f"no raw chunks under {root / 'raw'}")
    manifest = _manifest(root)
    requested = int(manifest["requested_replications_per_cell"])
    successes = (
        raw.loc[(raw["record_type"] == "target") & (raw["status"] == "success")].copy()
        if not raw.empty
        else pd.DataFrame(columns=["dgp", "N", "T", "target", "true_rank_vector"])
    )
    rows = []
    for keys, group in successes.groupby(["dgp", "N", "T", "target", "true_rank_vector"], dropna=False):
        errors = group["estimate"] - group["truth"]
        plugin_errors = group["plugin_estimate"] - group["truth"]
        corrected_errors = group["corrected_estimate"] - group["truth"]
        mc_sd = float(group["estimate"].std(ddof=1)) if len(group) > 1 else float("nan")
        mean_se = float(group["standard_error"].mean())
        size = float(group["centered_reject_5pct"].mean())
        coverage = float(group["covered_95pct"].mean())
        record = {
            "dgp": keys[0],
            "N": keys[1],
            "T": keys[2],
            "target": keys[3],
            "true_rank_vector": keys[4],
            "mean_truth": float(group["truth"].mean()),
            "mean_estimate": float(group["estimate"].mean()),
            "bias": float(errors.mean()),
            "rmse": float(np.sqrt(np.mean(errors**2))),
            "mc_sd": mc_sd,
            "mean_estimated_se": mean_se,
            "se_to_mc_sd": mean_se / mc_sd if mc_sd > 0 else float("nan"),
            "size_5pct": size,
            "coverage_95pct": coverage,
            "power_against_zero": float(group["reject_zero_5pct"].mean()),
            "mcse_size": float(np.sqrt(size * (1.0 - size) / len(group))),
            "mcse_coverage": float(np.sqrt(coverage * (1.0 - coverage) / len(group))),
            "successful_replications": int(group["replication"].nunique()),
            "requested_replications": requested,
            "plugin_mean_estimate": float(group["plugin_estimate"].mean()),
            "plugin_bias": float(plugin_errors.mean()),
            "plugin_rmse": float(np.sqrt(np.mean(plugin_errors**2))),
            "plugin_mean_se": float(group["plugin_standard_error"].mean()),
            "corrected_mean_estimate": float(group["corrected_estimate"].mean()),
            "corrected_bias": float(corrected_errors.mean()),
            "corrected_rmse": float(np.sqrt(np.mean(corrected_errors**2))),
            "corrected_mean_se": float(group["corrected_standard_error"].mean()),
            "mean_rank_runtime_seconds": float(group["rank_runtime_seconds"].mean()),
            "mean_inference_runtime_seconds": float(group["inference_runtime_seconds"].mean()),
            "mean_replication_runtime_seconds": float(group["replication_runtime_seconds"].mean()),
            "mean_best_neighbor_target_change": float(group["best_neighbor_target_change"].mean()),
            "spatial_c": float(group["spatial_c"].mean()) if "spatial_c" in group else np.nan,
        }
        relevant = raw.loc[
            (raw["dgp"] == keys[0])
            & (raw["N"] == keys[1])
            & (raw["T"] == keys[2])
            & (raw["true_rank_vector"] == keys[4])
            & raw["target"].isin([keys[3], "__replication__"])
        ]
        failures = relevant.drop_duplicates(["replication", "status"])["status"].value_counts()
        for code in FAILURE_CODES:
            record[f"failure_{code}"] = int(failures.get(code, 0)) if code != "success" else 0
        for truth_column in [column for column in group if column.endswith("_true")]:
            record[f"mean_{truth_column}"] = float(group[truth_column].mean())
            record[f"sd_{truth_column}"] = float(group[truth_column].std(ddof=1)) if len(group) > 1 else 0.0
        rows.append(record)
    summary = pd.DataFrame(rows)
    if not summary.empty:
        summary = summary.sort_values(["dgp", "N", "T", "true_rank_vector", "target"]).reset_index(drop=True)
    else:
        summary = pd.DataFrame(
            columns=[
                "dgp", "N", "T", "target", "true_rank_vector", "mean_truth",
                "mean_estimate", "bias", "rmse", "mc_sd", "mean_estimated_se",
                "se_to_mc_sd", "size_5pct", "coverage_95pct", "power_against_zero",
                "successful_replications", "requested_replications", "plugin_bias",
                "corrected_bias", "plugin_rmse", "corrected_rmse", "corrected_mean_se",
                "mean_rank_runtime_seconds", "mean_inference_runtime_seconds",
                "mean_replication_runtime_seconds", "spatial_c",
            ]
        )

    rank_rows = []
    if not rank_raw.empty:
        for keys, group in rank_raw.groupby(["dgp", "N", "T", "true_rank_vector"], dropna=False):
            successful = group.loc[group["status"] == "success"]
            denominator = len(group)
            rank_rows.append(
                {
                    "dgp": keys[0],
                    "N": keys[1],
                    "T": keys[2],
                    "true_rank_vector": keys[3],
                    "candidate_coverage": float(group["true_rank_in_candidates"].mean()),
                    "exact_rank_recovery": float(group["exact_rank_recovery"].mean()),
                    "A_underselection": float(group["A_underselected"].mean()),
                    "A_overselection": float(group["A_overselected"].mean()),
                    "B_underselection": float(group["B_underselected"].mean()),
                    "B_overselection": float(group["B_overselected"].mean()),
                    "H_underselection": float(group["H_underselected"].mean()),
                    "H_overselection": float(group["H_overselected"].mean()),
                    "zero_rank_recovery": float(group["zero_rank_recovery"].mean()),
                    "cap_hit_rate": float(group["rank_at_cap"].mean()),
                    "mean_candidate_count": float(group["candidate_count_final"].mean()),
                    "mean_ic_gap": float(group["ic_gap"].mean()),
                    "mean_rank_runtime_seconds": float(group["rank_runtime_seconds"].mean()),
                    "successful_replications": int(len(successful)),
                    "requested_replications": requested,
                    "rank_records": denominator,
                }
            )
    rank_summary = pd.DataFrame(rank_rows)
    if not rank_summary.empty:
        rank_summary = rank_summary.sort_values(["true_rank_vector", "dgp", "N", "T"])

    destination = root / "summary"
    destination.mkdir(exist_ok=True)
    for name, frame in (("mc_summary", summary), ("rank_summary", rank_summary)):
        frame.to_parquet(destination / f"{name}.parquet", index=False)
        frame.to_csv(destination / f"{name}.csv", index=False)
    return summary


def _escape(value: object) -> str:
    if value is None or (isinstance(value, float) and not np.isfinite(value)):
        return "--"
    text = str(value)
    for old, new in (("%", "\\%"), ("&", "\\&"), ("_", "\\_")):
        text = text.replace(old, new)
    return text


def _number(value: object, digits: int = 3) -> str:
    if value is None or pd.isna(value):
        return "--"
    if isinstance(value, (int, np.integer)):
        return str(int(value))
    return f"{float(value):.{digits}f}"


def _longtable(
    panels: list[tuple[str, pd.DataFrame]],
    columns: list[tuple[str, str]],
    *,
    caption: str,
    label: str,
    notes: str,
    first_text_column: str | None = None,
) -> str:
    alignment = "rr" + ("l" if first_text_column else "") + "r" * (len(columns) - 2 - bool(first_text_column))
    lines = [
        "\\begingroup",
        "\\small",
        f"\\begin{{longtable}}{{{alignment}}}",
        f"\\caption{{{caption}}}\\label{{{label}}}\\\\",
        "\\toprule",
        " & ".join(title for _, title in columns) + " \\\\",
        "\\midrule",
        "\\endfirsthead",
        "\\toprule",
        " & ".join(title for _, title in columns) + " \\\\",
        "\\midrule",
        "\\endhead",
    ]
    for title, frame in panels:
        lines.append(f"\\multicolumn{{{len(columns)}}}{{l}}{{\\textit{{{title}}}}} \\\\")
        for row in frame.to_dict("records"):
            values = []
            for key, _ in columns:
                values.append(_escape(row.get(key)) if key == first_text_column else _number(row.get(key)))
            lines.append(" & ".join(values) + " \\\\")
        lines.append("\\addlinespace")
    lines.extend(
        [
            "\\bottomrule",
            f"\\multicolumn{{{len(columns)}}}{{p{{0.96\\linewidth}}}}{{\\footnotesize \\textit{{Notes:}} {notes}}} \\\\",
            "\\end{longtable}",
            "\\endgroup",
            "",
        ]
    )
    return "\n".join(lines)


def _write_artifact(base: Path, data: pd.DataFrame, latex: str) -> list[Path]:
    paths = [base.with_suffix(".csv"), base.with_suffix(".parquet"), base.with_suffix(".tex")]
    data.to_csv(paths[0], index=False)
    data.to_parquet(paths[1], index=False)
    paths[2].write_text(latex, encoding="utf-8")
    return paths


def _performance_table(summary: pd.DataFrame, dgp: int) -> tuple[pd.DataFrame, str]:
    data = summary.loc[summary["dgp"] == dgp].copy()
    panels = []
    for letter, target in zip(string.ascii_uppercase, TARGET_TITLES, strict=False):
        subset = data.loc[data["target"] == target]
        if not subset.empty:
            panels.append((f"Panel {letter}. {TARGET_TITLES[target]}", subset))
    columns = [
        ("N", "$N$"), ("T", "$T$"), ("mean_truth", "True"),
        ("mean_estimate", "Mean est."), ("bias", "Bias"), ("rmse", "RMSE"),
        ("mc_sd", "MC sd"), ("mean_estimated_se", "Mean s.e."),
        ("se_to_mc_sd", "s.e./MC sd"), ("size_5pct", "Size"),
        ("coverage_95pct", "Coverage"), ("power_against_zero", "Power"),
        ("successful_replications", "Success"), ("requested_replications", "Requested"),
    ]
    latex = _longtable(
        panels,
        columns,
        caption=f"Monte Carlo performance for DGP {dgp}",
        label=f"tab:mc-dgp{dgp}",
        notes=(
            "Rates are on the 0--1 scale. Broad targets report the prescribed two-way corrected "
            "estimator; local and fixed-time targets report the full-panel plug-in. DGPs 2--4 "
            "use the logarithmic Bartlett spatial cutoff."
        ),
    )
    return data, latex


def make_tables(root: str | Path) -> list[Path]:
    root = Path(root)
    summary_path = root / "summary" / "mc_summary.parquet"
    summary = pd.read_parquet(summary_path) if summary_path.exists() else aggregate_run(root)
    rank_summary = pd.read_parquet(root / "summary" / "rank_summary.parquet")
    table_root = root / "tables"
    table_root.mkdir(exist_ok=True)
    outputs: list[Path] = []

    for dgp in range(1, 5):
        data, latex = _performance_table(summary, dgp)
        outputs.extend(_write_artifact(table_root / f"tab_mc_dgp{dgp}", data, latex))
    main_tex = "\n".join(f"\\input{{tab_mc_dgp{dgp}.tex}}" for dgp in range(1, 5)) + "\n"
    main_data = summary.copy()
    outputs.extend(_write_artifact(table_root / "tab_mc_main_summary", main_data, main_tex))

    calibration = pd.read_parquet(root / "calibration.parquet").copy()
    calibration["requested_replications"] = _manifest(root)["requested_replications_per_cell"]
    coeff_columns = [
        ("dgp", "DGP"), ("n", "$N$"), ("t", "$T$"), ("mean_a", "Mean A"),
        ("sd_a", "sd A"), ("min_a", "Min A"), ("max_a", "Max A"),
        ("mean_b", "Mean B"), ("sd_b", "sd B"), ("min_b", "Min B"),
        ("max_b", "Max B"), ("achieved_r2", "$R^2$"), ("c_h", "$c_H$"),
        ("c_xi", "$c_\\xi$"), ("mean_c_a", "Mean $c_a$"),
        ("requested_replications", "Reps"),
    ]
    coeff_tex = _longtable(
        [("Panel A. Calibration and coefficient summaries", calibration)],
        coeff_columns,
        caption="Calibration and coefficient summaries",
        label="tab:mc-coeff-summary",
        notes="Calibration draws are deterministic and separate from reported replications.",
    )
    outputs.extend(_write_artifact(table_root / "tab_mc_coeff_summary", calibration, coeff_tex))

    broad = summary.loc[summary["target"].map(TARGET_TITLES).notna() & summary["target"].str.contains("full_mean|time_average")].copy()
    bias_panels = [
        (f"Panel {letter}. {TARGET_TITLES[target]}", broad.loc[broad["target"] == target])
        for letter, target in zip(string.ascii_uppercase, broad["target"].drop_duplicates(), strict=False)
    ]
    bias_columns = [
        ("N", "$N$"), ("T", "$T$"), ("mean_truth", "True"),
        ("plugin_bias", "Plug-in bias"), ("corrected_bias", "Corrected bias"),
        ("plugin_rmse", "Plug-in RMSE"), ("corrected_rmse", "Corrected RMSE"),
        ("corrected_mean_se", "Corrected s.e."), ("coverage_95pct", "Coverage"),
        ("successful_replications", "Success"), ("requested_replications", "Requested"),
    ]
    bias_tex = _longtable(
        bias_panels,
        bias_columns,
        caption="Effect of the two-way split-panel bias correction",
        label="tab:mc-bias-correction",
        notes="The primary estimator for every broad target is the corrected estimator.",
    )
    outputs.extend(_write_artifact(table_root / "tab_mc_bias_correction", broad, bias_tex))

    group_specs = [
        ("A fixed-time group truths", ["A_G1_fixed_time_true", "A_G2_fixed_time_true", "A_G2_minus_G1_fixed_time_true"]),
        ("A time-averaged group truths", ["A_G1_time_average_true", "A_G2_time_average_true", "A_G2_minus_G1_time_average_true"]),
        ("B fixed-time group truths", ["B_G1_fixed_time_true", "B_G2_fixed_time_true", "B_G2_minus_G1_fixed_time_true"]),
        ("B time-averaged group truths", ["B_G1_time_average_true", "B_G2_time_average_true", "B_G2_minus_G1_time_average_true"]),
    ]
    group_labels = {
        "A_G1_fixed_time_true": "$A$: group 1 mean",
        "A_G2_fixed_time_true": "$A$: group 2 mean",
        "A_G2_minus_G1_fixed_time_true": "$A$: group 2 minus group 1",
        "A_G1_time_average_true": "$A$: group 1 mean",
        "A_G2_time_average_true": "$A$: group 2 mean",
        "A_G2_minus_G1_time_average_true": "$A$: group 2 minus group 1",
        "B_G1_fixed_time_true": "$B$: group 1 mean",
        "B_G2_fixed_time_true": "$B$: group 2 mean",
        "B_G2_minus_G1_fixed_time_true": "$B$: group 2 minus group 1",
        "B_G1_time_average_true": "$B$: group 1 mean",
        "B_G2_time_average_true": "$B$: group 2 mean",
        "B_G2_minus_G1_time_average_true": "$B$: group 2 minus group 1",
    }
    dgp4_base = summary.loc[summary["dgp"] == 4].drop_duplicates(["N", "T"])
    group_rows, group_panels = [], []
    for letter, (title, names) in zip(string.ascii_uppercase, group_specs, strict=False):
        panel_rows = []
        for _, row in dgp4_base.iterrows():
            for name in names:
                item = {
                    "N": row["N"],
                    "T": row["T"],
                    "quantity": group_labels[name],
                    "truth_field": name,
                    "true": row.get(f"mean_{name}"),
                }
                panel_rows.append(item)
                group_rows.append(item)
        group_panels.append((f"Panel {letter}. {title}", pd.DataFrame(panel_rows)))
    group_tex = _longtable(
        group_panels,
        [("N", "$N$"), ("T", "$T$"), ("quantity", "Quantity"), ("true", "True")],
        caption="DGP 4 post-scaling group truths",
        label="tab:mc-dgp4-groups",
        notes="A truths are after the common stability rescaling; the raw A contrasts and $c_a$ remain in the accompanying data file.",
        first_text_column="quantity",
    )
    outputs.extend(_write_artifact(table_root / "tab_mc_dgp4_groups", pd.DataFrame(group_rows), group_tex))

    rank_panels = []
    for letter, true_rank in zip(string.ascii_uppercase, rank_summary["true_rank_vector"].drop_duplicates(), strict=False):
        rank_panels.append((f"Panel {letter}. True rank vector {true_rank}", rank_summary.loc[rank_summary["true_rank_vector"] == true_rank]))
    rank_columns = [
        ("dgp", "DGP"), ("N", "$N$"), ("T", "$T$"),
        ("candidate_coverage", "Candidate cover"), ("exact_rank_recovery", "Exact recover"),
        ("A_underselection", "A under"), ("A_overselection", "A over"),
        ("B_underselection", "B under"), ("B_overselection", "B over"),
        ("H_underselection", "H under"), ("H_overselection", "H over"),
        ("zero_rank_recovery", "Zero recover"), ("cap_hit_rate", "Cap hit"),
        ("mean_candidate_count", "Candidates"), ("mean_ic_gap", "IC gap"),
        ("successful_replications", "Success"), ("requested_replications", "Requested"),
    ]
    rank_tex = _longtable(
        rank_panels,
        rank_columns,
        caption="Rank-selection performance",
        label="tab:mc-rank",
        notes="Rows are unique by DGP, panel cell, and true-rank design; no inference-target duplication is present.",
    )
    outputs.extend(_write_artifact(table_root / "tab_mc_rank", rank_summary, rank_tex))

    computation = rank_summary[
        ["dgp", "N", "T", "true_rank_vector", "mean_rank_runtime_seconds", "mean_candidate_count", "mean_ic_gap", "successful_replications", "requested_replications"]
    ].copy()
    comp_panels = [("Panel A. Replication-level rank computation", computation)]
    comp_tex = _longtable(
        comp_panels,
        [("dgp", "DGP"), ("N", "$N$"), ("T", "$T$"), ("true_rank_vector", "True rank"),
         ("mean_rank_runtime_seconds", "Rank sec."), ("mean_candidate_count", "Candidates"),
         ("mean_ic_gap", "IC gap"), ("successful_replications", "Success"),
         ("requested_replications", "Requested")],
        caption="Computation and numerical diagnostics",
        label="tab:mc-computation",
        notes="Quantities are replication-level and are not duplicated by inference target.",
        first_text_column="true_rank_vector",
    )
    outputs.extend(_write_artifact(table_root / "tab_mc_computation", computation, comp_tex))

    spatial = summary.loc[summary["dgp"].isin([2, 3, 4])].copy()
    spatial_tex = _longtable(
        [("Panel A. Available spatial-cutoff specifications", spatial)],
        [("N", "$N$"), ("T", "$T$"), ("spatial_c", "$c_{sp}$"),
         ("mean_truth", "True"), ("bias", "Bias"), ("coverage_95pct", "Coverage"),
         ("successful_replications", "Success"), ("requested_replications", "Requested")],
        caption="Spatial-HAC sensitivity",
        label="tab:mc-spatial-sensitivity",
        notes="Rows appear only for sensitivity specifications present in the aggregated run; the Bartlett cutoff is logarithmic.",
    )
    outputs.extend(_write_artifact(table_root / "tab_mc_spatial_sensitivity", spatial, spatial_tex))
    return outputs
