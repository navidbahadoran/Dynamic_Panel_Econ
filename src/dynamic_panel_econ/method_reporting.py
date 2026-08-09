"""Method-comparison tables and figures from reconciled Monte Carlo records."""

from __future__ import annotations

import json
from collections.abc import Sequence
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd

from .mc_accounting import (
    PRIMARY_STATUSES,
    reconcile_fit_rows,
    reconcile_matched_draws,
    reconcile_summary,
    summarize_accounting,
    summarize_power,
)

TABLE_NAMES = (
    "dgp1", "dgp2", "dgp3", "dgp4", "bias_correction", "rank",
    "optimization", "failure_accounting", "target_regularities", "power", "runtime",
)
FIGURE_NAMES = (
    "bias", "rmse", "coverage", "failure_rate", "retained_share", "runtime",
    "interval_length", "power_A", "power_B",
)


def _read_many(roots: Sequence[Path], filename: str) -> pd.DataFrame:
    files = []
    for root in roots:
        direct = root / filename
        files.extend([direct] if direct.exists() else root.glob(f"**/{filename}"))
    frames = [pd.read_parquet(path) for path in dict.fromkeys(files)]
    return pd.concat(frames, ignore_index=True) if frames else pd.DataFrame()


def _write_exclusive(path: Path, content: str | bytes, *, overwrite: bool) -> None:
    if path.exists() and not overwrite:
        raise FileExistsError(f"refusing to overwrite existing reporting artifact: {path}")
    path.parent.mkdir(parents=True, exist_ok=True)
    if isinstance(content, bytes):
        path.write_bytes(content)
    else:
        path.write_text(content, encoding="utf-8")


def _latex_table(frame: pd.DataFrame, *, caption: str, label: str, panel: str) -> str:
    columns = list(frame.columns)
    alignment = "l" * len(columns)
    lines = [
        r"\begingroup", r"\small", rf"\begin{{longtable}}{{{alignment}}}",
        rf"\caption{{{caption}}}\label{{{label}}}\\", r"\toprule",
        " & ".join(column.replace("_", r"\_") for column in columns) + r" \\",
        r"\midrule", r"\endfirsthead", r"\toprule",
        " & ".join(column.replace("_", r"\_") for column in columns) + r" \\",
        r"\midrule", r"\endhead", rf"\multicolumn{{{len(columns)}}}{{l}}{{\textit{{{panel}}}}} \\",
    ]
    for row in frame.to_dict("records"):
        values = []
        for column in columns:
            value = row[column]
            if value is None or (isinstance(value, float) and not np.isfinite(value)):
                values.append("--")
            elif isinstance(value, float):
                values.append(f"{value:.4g}")
            else:
                values.append(str(value).replace("_", r"\_"))
        lines.append(" & ".join(values) + r" \\")
    lines += [r"\bottomrule", r"\end{longtable}", r"\endgroup", ""]
    return "\n".join(lines)


def _performance_columns(summary: pd.DataFrame) -> pd.DataFrame:
    wanted = [
        "method", "target", "N", "T", "mean_truth", "mean_estimate", "bias", "rmse",
        "mc_sd", "mean_se", "coverage", "rejection_probability", "R_attempted", "R_point",
        "R_inference", "point_retained_share", "inference_retained_share",
    ]
    return summary[[column for column in wanted if column in summary]].copy()


def _table_data(
    name: str,
    summary: pd.DataFrame,
    replication: pd.DataFrame,
    fit: pd.DataFrame,
    inference: pd.DataFrame,
    target_records: pd.DataFrame,
) -> pd.DataFrame:
    if name.startswith("dgp"):
        return _performance_columns(summary.loc[summary["dgp"].eq(int(name[-1]))])
    if name == "failure_accounting":
        columns = ["dgp", "N", "T", "method", "target", "R_attempted", "R_point", "R_inference", "numerical_failure_rate", "target_support_failure_rate", "total_inference_failure_rate"]
        columns += [f"failure_{status}" for status in PRIMARY_STATUSES if f"failure_{status}" in summary]
        return summary[[column for column in columns if column in summary]]
    if name == "runtime":
        columns = ["dgp", "N", "T", "method", "target", "runtime_median", "runtime_mean", "runtime_p10", "runtime_p90", "runtime_p95"]
        return summary[[column for column in columns if column in summary]]
    if name == "optimization":
        if fit.empty:
            return fit
        keys = ["dgp", "N", "T", "method", "fit_type"]
        aggregations = {
            "fits": ("fit_type", "size"),
            "convergence_rate": ("convergence_flag", "mean"),
            "stationarity_pass_rate": ("stationarity_pass", "mean"),
            "median_runtime": ("runtime_seconds", "median"),
        }
        if "boundary_active" in fit:
            aggregations["boundary_active_rate"] = ("boundary_active", "mean")
        if "constrained_fallback_used" in fit:
            aggregations["constrained_fallback_rate"] = (
                "constrained_fallback_used",
                "mean",
            )
        if "constrained_solver_status" in fit:
            failures = fit["constrained_solver_status"].isin(
                {
                    "constrained_solver_failure",
                    "constrained_feasibility_failure",
                    "constrained_optimality_failure",
                    "nonfinite_constrained_solution",
                }
            )
            fit = fit.assign(_constrained_solver_failure=failures)
            aggregations["constrained_solver_failure_rate"] = (
                "_constrained_solver_failure",
                "mean",
            )
        if "coefficient_bound_hit" in fit:
            aggregations["legacy_coefficient_bound_hit_rate"] = (
                "coefficient_bound_hit",
                "mean",
            )
        return fit.groupby(keys, dropna=False).agg(**aggregations).reset_index()
    if name == "rank":
        columns = ["dgp", "N", "T", "method", "supplied_rank_vector", "selected_rank_vector", "cap_pilot_rank", "candidate_coverage", "primary_status"]
        return replication[[column for column in columns if column in replication]]
    if name == "target_regularities":
        columns = ["dgp", "N", "T", "method", "target", "target_tangent_norm", "target_supported", "tangent_gram_min_eigenvalue", "tangent_gram_condition_number"]
        return inference[[column for column in columns if column in inference]]
    if name == "power":
        return (
            summarize_power(target_records)
            if "nominal_delta" in target_records
            else pd.DataFrame()
        )
    if name == "bias_correction":
        columns = ["dgp", "N", "T", "method", "target", "phi_full", "phi_time_sum", "phi_unit_sum", "phi_corrected"]
        return inference[[column for column in columns if column in inference]]
    raise ValueError(f"unknown table: {name}")


def _plot_metric(summary: pd.DataFrame, name: str, destination: Path, *, overwrite: bool) -> None:
    try:
        import matplotlib.pyplot as plt
    except ImportError as exc:  # pragma: no cover - dependency error is environment-specific
        raise RuntimeError("matplotlib is required to render Monte Carlo figures") from exc
    metric = {
        "bias": "bias", "rmse": "rmse", "coverage": "coverage",
        "failure_rate": "total_inference_failure_rate",
        "retained_share": "inference_retained_share", "runtime": "runtime_median",
        "interval_length": "mean_interval_length",
    }[name]
    figure, axis = plt.subplots(figsize=(6.5, 4.0))
    keys = ["dgp", "method", "target"]
    for label, group in summary.groupby(keys, dropna=False, sort=True):
        ordered = group.sort_values("N")
        axis.plot(ordered["N"], ordered[metric], marker="o", label=" / ".join(map(str, label)))
        if name == "runtime" and {"runtime_p10", "runtime_p90"} <= set(group):
            axis.fill_between(ordered["N"], ordered["runtime_p10"], ordered["runtime_p90"], alpha=0.12)
    if name == "coverage":
        axis.axhline(0.95, color="black", linestyle="--", linewidth=0.8)
    axis.set_xlabel("N = T")
    axis.set_ylabel(metric.replace("_", " "))
    axis.legend(fontsize=6, ncol=2)
    figure.tight_layout()
    for extension in ("png", "pdf"):
        path = destination / f"{name}.{extension}"
        if path.exists() and not overwrite:
            plt.close(figure)
            raise FileExistsError(f"refusing to overwrite existing reporting artifact: {path}")
        path.parent.mkdir(parents=True, exist_ok=True)
        figure.savefig(path, dpi=300)
    plt.close(figure)


def _plot_power(power: pd.DataFrame, block: str, destination: Path, *, overwrite: bool) -> None:
    try:
        import matplotlib.pyplot as plt
    except ImportError as exc:  # pragma: no cover
        raise RuntimeError("matplotlib is required to render Monte Carlo figures") from exc
    data = power.loc[power["target"].astype(str).str.startswith(block)]
    figure, axis = plt.subplots(figsize=(6.5, 4.0))
    for method, group in data.groupby("method", dropna=False):
        ordered = group.sort_values("realized_true_contrast")
        axis.errorbar(
            ordered["realized_true_contrast"], ordered["rejection_probability"],
            yerr=1.96 * ordered["mc_binomial_se"], marker="o", label=method,
        )
    axis.axhline(0.05, color="black", linestyle="--", linewidth=0.8)
    axis.set_xlabel("Realized true post-scaling group contrast")
    axis.set_ylabel("Rejection probability")
    axis.legend()
    figure.tight_layout()
    for extension in ("png", "pdf"):
        path = destination / f"power_{block}.{extension}"
        if path.exists() and not overwrite:
            plt.close(figure)
            raise FileExistsError(f"refusing to overwrite existing reporting artifact: {path}")
        path.parent.mkdir(parents=True, exist_ok=True)
        figure.savefig(path, dpi=300)
    plt.close(figure)


def report_method_comparison(
    input_roots: Sequence[str | Path],
    output_dir: str | Path,
    *,
    tables: Sequence[str],
    figures: Sequence[str],
    overwrite: bool = False,
) -> dict[str, Any]:
    roots = [Path(root) for root in input_roots]
    output = Path(output_dir)
    replication = _read_many(roots, "replication_records.parquet")
    attempts = _read_many(roots, "attempted_replications.parquet")
    targets = _read_many(roots, "raw/*.parquet")
    if targets.empty:
        raw_files = [path for root in roots for path in root.glob("raw/*.parquet")]
        targets = pd.concat([pd.read_parquet(path) for path in raw_files], ignore_index=True) if raw_files else pd.DataFrame()
    target_records = targets.loc[targets.get("record_type", pd.Series(index=targets.index, dtype=object)).eq("target")].copy()
    fit = _read_many(roots, "fit_diagnostics.parquet")
    inference = _read_many(roots, "inference_diagnostics.parquet")
    if replication.empty:
        raise FileNotFoundError("no replication_records.parquet found")
    attempts = attempts if not attempts.empty else replication
    reconcile_matched_draws(attempts)
    if "expected_fit_count" in attempts and not fit.empty:
        keys = ["run_id", "dgp", "N", "T", "replication", "method", "expected_fit_count"]
        reconcile_fit_rows(fit, attempts[keys].drop_duplicates())
    target_names = sorted(
        set(target_records.get("target", pd.Series(dtype=str)).dropna().astype(str))
        | set(
            replication.loc[
                replication.get("target", pd.Series(index=replication.index, dtype=str)).ne(
                    "__replication__"
                ),
                "target",
            ].dropna().astype(str)
        )
    )
    summary = summarize_accounting(attempts, target_records, targets=target_names)
    reconcile_summary(summary)
    selected_tables = TABLE_NAMES if "all" in tables else tuple(tables)
    selected_figures = FIGURE_NAMES if "all" in figures else tuple(figures)
    for name in selected_tables:
        data = _table_data(name, summary, replication, fit, inference, target_records)
        latex = _latex_table(
            data,
            caption=f"Monte Carlo {name.replace('_', ' ')}",
            label=f"tab:mc-{name.replace('_', '-')}",
            panel="Panel A. Fixed-rank and selected-rank methods",
        )
        _write_exclusive(output / "tables" / f"tab_mc_{name}.tex", latex, overwrite=overwrite)
    power = summarize_power(target_records) if "nominal_delta" in target_records else pd.DataFrame()
    for name in selected_figures:
        if name.startswith("power_"):
            _plot_power(power, name[-1], output / "figures", overwrite=overwrite)
        else:
            _plot_metric(summary, name, output / "figures", overwrite=overwrite)
    summary_path = output / "summary.csv"
    if summary_path.exists() and not overwrite:
        raise FileExistsError(f"refusing to overwrite existing reporting artifact: {summary_path}")
    output.mkdir(parents=True, exist_ok=True)
    summary.to_csv(summary_path, index=False)
    summary.to_parquet(output / "summary.parquet", index=False)
    captions = "\n".join(f"- `{name}`: Monte Carlo {name.replace('_', ' ')} by method and sample size." for name in selected_figures) + "\n"
    _write_exclusive(output / "captions" / "figure_captions.md", captions, overwrite=overwrite)
    tex_captions = "\n".join(rf"\newcommand{{\caption{name.replace('_', '')}}}{{Monte Carlo {name.replace('_', ' ')}}}" for name in selected_figures) + "\n"
    _write_exclusive(output / "captions" / "figure_captions.tex", tex_captions, overwrite=overwrite)
    reconciliation = {
        "summary_rows": len(summary),
        "attempted": int(summary["R_attempted"].sum()),
        "point": int(summary["R_point"].sum()),
        "inference": int(summary["R_inference"].sum()),
        "status": "passed",
    }
    _write_exclusive(
        output / "audit" / "reconciliation_report.md",
        "# Reconciliation report\n\n```json\n" + json.dumps(reconciliation, indent=2) + "\n```\n",
        overwrite=overwrite,
    )
    return reconciliation
