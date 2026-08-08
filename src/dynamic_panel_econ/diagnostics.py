"""Pre-production diagnostics that leave maintained estimators and production config unchanged."""

from __future__ import annotations

from collections.abc import Callable
from dataclasses import dataclass
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd
from scipy.optimize import brentq

from .dgp import DGPParameters, _draw_raw
from .seeds import rng_for


@dataclass(frozen=True, slots=True)
class R2Curve:
    """A deterministic common-random-number pooled-R2 calibration curve."""

    c_h: float
    r2: Callable[[float], float]
    large_c_xi_floor: float


def build_r2_curve(
    dgp: int,
    n: int,
    t: int,
    master_seed: int,
    *,
    params: DGPParameters,
    pi_h: float,
    draws: int,
) -> R2Curve:
    """Construct the exact affine-in-``c_xi`` calibration curve for one cell."""

    raws = [
        _draw_raw(dgp, n, t, rng_for(master_seed, "calibration", dgp, n, t, j), params)
        for j in range(draws)
    ]
    observed = slice(params.burn_in, params.burn_in + t)
    var_u = float(np.mean([np.var(raw.u_tilde[:, observed]) for raw in raws]))
    var_h = float(np.mean([np.var(raw.h_raw[:, observed]) for raw in raws]))
    if var_u <= 0.0 or var_h <= 0.0:
        raise RuntimeError("nonpositive calibration variance")
    c_h = float(np.sqrt((pi_h / (1.0 - pi_h)) * var_u / var_h))

    # For each draw, the pooled total sum of squares of
    # y_base + c_xi*y_scale is exactly a quadratic in c_xi.  Retaining only
    # those four scalar coefficients makes the dense feasibility grid exact
    # without rescanning large arrays at every grid point.
    quadratic_components: list[tuple[float, float, float, float]] = []
    for raw in raws:
        length = raw.x.shape[1]
        y_base = np.empty((n, length), dtype=np.float64)
        y_scale = np.empty((n, length), dtype=np.float64)
        previous_base = np.zeros(n, dtype=np.float64)
        previous_scale = np.zeros(n, dtype=np.float64)
        shock_scale = c_h * raw.h_raw + raw.u_tilde
        for column in range(length):
            current_base = raw.a[:, column] * previous_base + raw.beta[:, column] * raw.x[:, column]
            current_scale = raw.a[:, column] * previous_scale + shock_scale[:, column]
            y_base[:, column] = current_base
            y_scale[:, column] = current_scale
            previous_base, previous_scale = current_base, current_scale
        base = y_base[:, observed]
        scale = y_scale[:, observed]
        primitive_u = raw.u_tilde[:, observed]
        base_centered = base - base.mean()
        scale_centered = scale - scale.mean()
        quadratic_components.append(
            (
                float(np.sum(base_centered**2)),
                float(2.0 * np.sum(base_centered * scale_centered)),
                float(np.sum(scale_centered**2)),
                float(np.sum(primitive_u**2)),
            )
        )

    def average_r2(c_xi: float) -> float:
        values = []
        for base_ss, cross, scale_ss, u_ss in quadratic_components:
            denominator = base_ss + c_xi * cross + c_xi**2 * scale_ss
            values.append(1.0 - c_xi**2 * u_ss / denominator)
        return float(np.mean(values))

    floors = [
        1.0 - u_ss / scale_ss
        for _, _, scale_ss, u_ss in quadratic_components
    ]
    return R2Curve(c_h, average_r2, float(np.mean(floors)))


def r2_feasibility_rows(
    config: dict[str, Any],
    *,
    targets: tuple[float, ...] = (0.60, 0.65, 0.70),
) -> pd.DataFrame:
    """Evaluate all configured DGP/cell calibration curves and target roots."""

    from .monte_carlo import _params

    params = _params(config)
    grid = np.geomspace(1e-8, 1e8, 321)
    rows: list[dict[str, Any]] = []
    for dgp in config["run"]["dgps"]:
        for n, t in config["run"]["cells"]:
            curve = build_r2_curve(
                int(dgp),
                int(n),
                int(t),
                int(config["run"]["master_seed"]),
                params=params,
                pi_h=float(config["dgp"]["pi_h"]),
                draws=int(config["dgp"]["calibration_draws"]),
            )
            r2_values = np.array([curve.r2(float(value)) for value in grid])
            row: dict[str, Any] = {
                "dgp": int(dgp),
                "N": int(n),
                "T": int(t),
                "calibration_draws": int(config["dgp"]["calibration_draws"]),
                "c_xi_lower": float(grid[0]),
                "r2_at_positive_lower_end": float(r2_values[0]),
                "large_c_xi_r2_floor": curve.large_c_xi_floor,
                "minimum_grid_r2": float(r2_values.min()),
                "minimum_grid_c_xi": float(grid[int(r2_values.argmin())]),
                "c_h": curve.c_h,
            }
            for target in targets:
                label = f"{target:.2f}".replace(".", "_")
                differences = r2_values - target
                brackets = [
                    (float(grid[index]), float(grid[index + 1]))
                    for index in range(len(grid) - 1)
                    if differences[index] == 0.0
                    or differences[index] * differences[index + 1] < 0.0
                ]
                feasible = bool(brackets)
                root = float("nan")
                residual = float("nan")
                if feasible:
                    left, right = brackets[0]
                    root = float(
                        brentq(
                            lambda value, curve=curve, target=target: curve.r2(value) - target,
                            left,
                            right,
                            xtol=1e-12,
                            rtol=1e-12,
                        )
                    )
                    residual = float(curve.r2(root) - target)
                row[f"target_{label}_feasible"] = feasible
                row[f"target_{label}_c_xi"] = root
                row[f"target_{label}_root_residual"] = residual
            rows.append(row)
    return pd.DataFrame(rows).sort_values(["dgp", "N", "T"]).reset_index(drop=True)


def write_r2_feasibility(frame: pd.DataFrame, output: Path) -> list[Path]:
    """Write machine-readable data and a compact journal-style longtable."""

    output.mkdir(parents=True, exist_ok=True)
    csv_path = output / "r2_feasibility.csv"
    parquet_path = output / "r2_feasibility.parquet"
    tex_path = output / "tab_mc_r2_feasibility.tex"
    frame.to_csv(csv_path, index=False)
    frame.to_parquet(parquet_path, index=False)
    lines = [
        "\\begingroup",
        "\\small",
        "\\begin{longtable}{rrrrrrrrrrrrrr}",
        "\\caption{Pooled-$R^2$ calibration feasibility}\\label{tab:mc-r2-feasibility}\\\\",
        "\\toprule",
        "DGP & $N$ & $T$ & Lower $R^2$ & Floor $R^2$ & Grid min. & "
        "$0.60$ & $c_\\xi$ & $0.65$ & $c_\\xi$ & $0.70$ & $c_\\xi$ & Draws & $c_H$ \\\\",
        "\\midrule",
        "\\endfirsthead",
        "\\toprule",
        "DGP & $N$ & $T$ & Lower $R^2$ & Floor $R^2$ & Grid min. & "
        "$0.60$ & $c_\\xi$ & $0.65$ & $c_\\xi$ & $0.70$ & $c_\\xi$ & Draws & $c_H$ \\\\",
        "\\midrule",
        "\\endhead",
    ]
    for row in frame.itertuples(index=False):
        cells = [
            str(row.dgp),
            str(row.N),
            str(row.T),
            f"{row.r2_at_positive_lower_end:.4f}",
            f"{row.large_c_xi_r2_floor:.4f}",
            f"{row.minimum_grid_r2:.4f}",
        ]
        for label in ("0_60", "0_65", "0_70"):
            feasible = getattr(row, f"target_{label}_feasible")
            root = getattr(row, f"target_{label}_c_xi")
            cells.extend(["Yes" if feasible else "No", f"{root:.4f}" if feasible else "--"])
        cells.extend([str(row.calibration_draws), f"{row.c_h:.4f}"])
        lines.append(" & ".join(cells) + r" \\")
    lines.extend(
        [
            "\\bottomrule",
            "\\multicolumn{14}{p{0.96\\linewidth}}{\\footnotesize \\textit{Notes:} "
            "The floor is the numerical affine large-$c_\\xi$ limit. Roots use deterministic "
            r"calibration draws independent of Monte Carlo replications.} \\",
            "\\end{longtable}",
            "\\endgroup",
        ]
    )
    tex_path.write_text("\n".join(lines) + "\n", encoding="utf-8")
    return [csv_path, parquet_path, tex_path]


def aggregate_rank_pilot(
    rank_rows: pd.DataFrame,
    *,
    dgps: list[int],
    cells: list[list[int]],
    true_rank_vectors: list[list[int]],
    replications: int,
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    """Create complete replication, baseline, and sensitivity rank-pilot data."""

    import json

    planned = pd.DataFrame(
        [
            {
                "dgp": int(dgp),
                "N": int(n),
                "T": int(t),
                "replication": replication,
                "true_rank_vector": json.dumps(tuple(int(value) for value in true_rank)),
            }
            for dgp in dgps
            for n, t in cells
            for true_rank in true_rank_vectors
            for replication in range(replications)
        ]
    )
    replication_frame = planned.merge(
        rank_rows,
        how="left",
        on=["dgp", "N", "T", "replication", "true_rank_vector"],
        suffixes=("", "_rank"),
    )
    replication_frame["numerical_fit_failure"] = replication_frame[
        "selected_rank_vector"
    ].isna()
    coverage = replication_frame["true_rank_in_candidates"].fillna(False).astype(bool)
    exact = replication_frame["exact_rank_recovery"].fillna(False).astype(bool)
    replication_frame["rank_failure_class"] = np.select(
        [
            replication_frame["numerical_fit_failure"],
            exact,
            ~coverage,
        ],
        ["numerical_fit_failure", "exact_recovery", "candidate_coverage_failure"],
        default="ic_choice_failure",
    )

    rows: list[dict[str, Any]] = []
    sensitivity_rows: list[dict[str, Any]] = []
    group_columns = ["dgp", "N", "T", "true_rank_vector"]
    for keys, group in replication_frame.groupby(group_columns, sort=True, dropna=False):
        available = group[~group["numerical_fit_failure"]]
        incorrect_covered = available[
            available["true_rank_in_candidates"].fillna(False)
            & ~available["exact_rank_recovery"].fillna(False)
        ]
        distribution = (
            available["selected_rank_vector"].value_counts(normalize=True).sort_index().to_dict()
        )
        row: dict[str, Any] = dict(zip(group_columns, keys, strict=True))
        row.update(
            {
                "requested_replications": len(group),
                "rank_records": len(available),
                "candidate_coverage": float(group["true_rank_in_candidates"].fillna(False).mean()),
                "exact_rank_recovery": float(group["exact_rank_recovery"].fillna(False).mean()),
                "A_underselection": float(available["A_underselected"].mean()),
                "A_overselection": float(available["A_overselected"].mean()),
                "B_underselection": float(available["B_underselected"].mean()),
                "B_overselection": float(available["B_overselected"].mean()),
                "H_underselection": float(available["H_underselected"].mean()),
                "H_overselection": float(available["H_overselected"].mean()),
                "rank_zero_recovery": float(available["zero_rank_recovery"].mean()),
                "candidate_coverage_failure_rate": float(
                    (group["rank_failure_class"] == "candidate_coverage_failure").mean()
                ),
                "ic_choice_failure_rate": float(
                    (group["rank_failure_class"] == "ic_choice_failure").mean()
                ),
                "numerical_fit_failure_rate": float(group["numerical_fit_failure"].mean()),
                "mean_selected_minus_true_ic_conditional_incorrect": float(
                    incorrect_covered["selected_minus_true_ic"].mean()
                ),
                "selected_rank_distribution": json.dumps(distribution, sort_keys=True),
            }
        )
        rows.append(row)
        specifications = {"baseline": "selected_rank_vector"}
        specifications.update(
            {
                f"ic_multiplier_{multiplier}": f"ic_multiplier_{multiplier}_selected_rank"
                for multiplier in ("0.5", "1.0", "2.0")
            }
        )
        specifications.update(
            {
                f"threshold_multiplier_{multiplier}": (
                    f"threshold_multiplier_{multiplier}_selected_rank"
                )
                for multiplier in ("0.5", "1.0", "2.0")
            }
        )
        for specification, column in specifications.items():
            selected_values = group[column] if column in group else pd.Series(index=group.index, dtype=object)
            is_exact = selected_values == group["true_rank_vector"]
            sensitivity_rows.append(
                {
                    **dict(zip(group_columns, keys, strict=True)),
                    "specification": specification,
                    "exact_rank_recovery": float(is_exact.fillna(False).mean()),
                    "available_replications": int(selected_values.notna().sum()),
                    "requested_replications": len(group),
                    "selected_rank_distribution": json.dumps(
                        selected_values.dropna().value_counts(normalize=True).sort_index().to_dict(),
                        sort_keys=True,
                    ),
                }
            )
    return (
        replication_frame,
        pd.DataFrame(rows),
        pd.DataFrame(sensitivity_rows),
    )


def write_rank_pilot(
    replications: pd.DataFrame,
    summary: pd.DataFrame,
    sensitivity: pd.DataFrame,
    output: Path,
) -> list[Path]:
    """Write rank-pilot data and a panel-oriented journal table."""

    output.mkdir(parents=True, exist_ok=True)
    paths: list[Path] = []
    for stem, frame in (
        ("rank_pilot_replications", replications),
        ("rank_pilot_summary", summary),
        ("rank_pilot_sensitivity", sensitivity),
    ):
        csv_path = output / f"{stem}.csv"
        parquet_path = output / f"{stem}.parquet"
        frame.to_csv(csv_path, index=False)
        frame.to_parquet(parquet_path, index=False)
        paths.extend([csv_path, parquet_path])
    tex_path = output / "tab_mc_rank_pilot.tex"
    newline = r"\\"
    lines = [
        r"\begingroup",
        r"\small",
        r"\begin{longtable}{rrlrrrrrrrrrrrr}",
        r"\caption{Medium-sample rank-selection diagnostic}\label{tab:mc-rank-pilot}" + newline,
        r"\toprule",
        r"$N$ & $T$ & Specification & Coverage & Exact & $A<$ & $A>$ & $B<$ & $B>$ & $H<$ & $H>$ & Zero & Cov. fail & IC fail & Num. fail " + newline,
        r"\midrule",
        r"\endfirsthead",
        r"\toprule",
        r"$N$ & $T$ & Specification & Coverage & Exact & $A<$ & $A>$ & $B<$ & $B>$ & $H<$ & $H>$ & Zero & Cov. fail & IC fail & Num. fail " + newline,
        r"\midrule",
        r"\endhead",
    ]
    panel_letters = "ABCDEFGHIJKLMNOPQRSTUVWXYZ"
    for letter, true_rank in zip(
        panel_letters, summary["true_rank_vector"].drop_duplicates(), strict=False
    ):
        lines.append(
            rf"\multicolumn{{15}}{{l}}{{\textit{{Panel {letter}. True rank {true_rank}}}}} "
            + newline
        )
        for row in summary[summary["true_rank_vector"] == true_rank].itertuples(index=False):
            values = [
                str(row.N), str(row.T), "Baseline", f"{row.candidate_coverage:.3f}",
                f"{row.exact_rank_recovery:.3f}", f"{row.A_underselection:.3f}",
                f"{row.A_overselection:.3f}", f"{row.B_underselection:.3f}",
                f"{row.B_overselection:.3f}", f"{row.H_underselection:.3f}",
                f"{row.H_overselection:.3f}", f"{row.rank_zero_recovery:.3f}",
                f"{row.candidate_coverage_failure_rate:.3f}",
                f"{row.ic_choice_failure_rate:.3f}", f"{row.numerical_fit_failure_rate:.3f}",
            ]
            lines.append(" & ".join(values) + " " + newline)
            cell_sensitivity = sensitivity[
                (sensitivity["true_rank_vector"] == true_rank)
                & (sensitivity["N"] == row.N)
                & (sensitivity["T"] == row.T)
                & (sensitivity["specification"] != "baseline")
            ]
            for sens in cell_sensitivity.itertuples(index=False):
                values = [
                    "", "", str(sens.specification).replace("_", r"\_"), "--",
                    f"{sens.exact_rank_recovery:.3f}", *(["--"] * 10),
                ]
                lines.append(" & ".join(values) + " " + newline)
        lines.append(r"\addlinespace")
    lines.extend(
        [
            r"\bottomrule",
            r"\multicolumn{15}{p{0.96\linewidth}}{\footnotesize \textit{Notes:} Rates use requested replications; under/over-selection rates condition on an available rank fit. Sensitivity rows report exact recovery only.} " + newline,
            r"\end{longtable}",
            r"\endgroup",
        ]
    )
    tex_path.write_text("\n".join(lines) + "\n", encoding="utf-8")
    paths.append(tex_path)
    return paths


def aggregate_riesz_diagnostic(
    raw_rows: pd.DataFrame,
    *,
    dgps: list[int],
    cells: list[list[int]],
    targets: list[str],
    replications: int,
    target_rayleigh_floor: float,
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    """Complete the target grid and extract target-specific instability events."""

    import json

    target_rows = raw_rows[raw_rows.get("record_type") == "target"].copy()
    failures = raw_rows[raw_rows.get("record_type") == "failure"].copy()
    failure_status = (
        failures.groupby(["dgp", "N", "T", "replication"], as_index=False)["status"].first()
        if not failures.empty
        else pd.DataFrame(columns=["dgp", "N", "T", "replication", "status"])
    )
    planned = pd.DataFrame(
        [
            {
                "dgp": int(dgp),
                "N": int(n),
                "T": int(t),
                "replication": replication,
                "target": target,
            }
            for dgp in dgps
            for n, t in cells
            for replication in range(replications)
            for target in targets
        ]
    )
    complete = planned.merge(
        target_rows,
        how="left",
        on=["dgp", "N", "T", "replication", "target"],
        suffixes=("", "_result"),
    ).merge(
        failure_status.rename(columns={"status": "replication_failure_status"}),
        how="left",
        on=["dgp", "N", "T", "replication"],
    )
    complete["status"] = complete["status"].fillna(
        complete["replication_failure_status"].fillna("missing_target_record")
    )
    events: list[dict[str, Any]] = []
    for row in target_rows.itertuples(index=False):
        base = {
            "dgp": row.dgp,
            "replication": row.replication,
            "N": row.N,
            "T": row.T,
            "target": row.target,
            "selected_rank": row.selected_rank,
            "selected_ic_gap": row.ic_gap,
            "cap_status": row.rank_at_cap,
            "coefficient_bound_status": row.coefficient_bound_active,
        }
        rayleigh = row.riesz_target_rayleigh_quotient
        if row.status == "riesz_target_instability" and (
            not np.isfinite(rayleigh) or rayleigh < target_rayleigh_floor
        ):
            events.append(
                {
                    **base,
                    "system": "full_panel",
                    "event_type": "riesz_target_instability",
                    "riesz_target_rayleigh_quotient": rayleigh,
                    "tangent_gram_smallest_eigenvalue": getattr(
                        row, "tangent_gram_smallest_eigenvalue", np.nan
                    ),
                    "tangent_gram_largest_eigenvalue": getattr(
                        row, "tangent_gram_largest_eigenvalue", np.nan
                    ),
                    "tangent_gram_condition_number": getattr(
                        row, "tangent_gram_condition_number", np.nan
                    ),
                    "tangent_gram_eigensolver_status": getattr(
                        row, "tangent_gram_eigensolver_status", "not_requested"
                    ),
                    "riesz_equation_residual": row.riesz_equation_residual,
                    "target_tangent_norm": row.riesz_target_tangent_norm,
                    "target_status": row.status,
                }
            )
        split_text = getattr(row, "split_diagnostics_json", None)
        if isinstance(split_text, str) and split_text:
            for split in json.loads(split_text):
                split_rayleigh = split.get("riesz_target_rayleigh_quotient", np.nan)
                if row.status == "split_riesz_target_instability" and (
                    not np.isfinite(split_rayleigh)
                    or split_rayleigh < target_rayleigh_floor
                ):
                    events.append(
                        {
                            **base,
                            "system": f"{split['kind']}_split_{split['part']}",
                            "event_type": "riesz_target_instability",
                            "riesz_target_rayleigh_quotient": split_rayleigh,
                            "tangent_gram_smallest_eigenvalue": split.get(
                                "tangent_gram_smallest_eigenvalue", np.nan
                            ),
                            "tangent_gram_largest_eigenvalue": split.get(
                                "tangent_gram_largest_eigenvalue", np.nan
                            ),
                            "tangent_gram_condition_number": split.get(
                                "tangent_gram_condition_number", np.nan
                            ),
                            "riesz_equation_residual": split.get("riesz_equation_residual"),
                            "target_tangent_norm": split.get("target_support_norm"),
                            "target_status": row.status,
                        }
                    )
        if row.status in {
            "tangent_gram_eigensolver_failure",
            "tangent_gram_nearly_singular",
        }:
            events.append(
                {
                    **base,
                    "system": "full_panel",
                    "event_type": row.status,
                    "riesz_target_rayleigh_quotient": rayleigh,
                    "tangent_gram_smallest_eigenvalue": getattr(
                        row, "tangent_gram_smallest_eigenvalue", np.nan
                    ),
                    "tangent_gram_largest_eigenvalue": getattr(
                        row, "tangent_gram_largest_eigenvalue", np.nan
                    ),
                    "tangent_gram_condition_number": getattr(
                        row, "tangent_gram_condition_number", np.nan
                    ),
                    "tangent_gram_eigensolver_status": getattr(
                        row, "tangent_gram_eigensolver_status", "failed"
                    ),
                    "riesz_equation_residual": row.riesz_equation_residual,
                    "target_tangent_norm": row.riesz_target_tangent_norm,
                    "target_status": row.status,
                }
            )
        if row.status in {
            "split_tangent_gram_eigensolver_failure",
            "split_tangent_gram_nearly_singular",
        }:
            for split in json.loads(split_text or "[]"):
                events.append(
                    {
                        **base,
                        "system": f"{split['kind']}_split_{split['part']}",
                        "event_type": row.status,
                        "riesz_target_rayleigh_quotient": split.get(
                            "riesz_target_rayleigh_quotient", np.nan
                        ),
                        "tangent_gram_smallest_eigenvalue": split.get(
                            "tangent_gram_smallest_eigenvalue", np.nan
                        ),
                        "tangent_gram_largest_eigenvalue": split.get(
                            "tangent_gram_largest_eigenvalue", np.nan
                        ),
                        "tangent_gram_condition_number": split.get(
                            "tangent_gram_condition_number", np.nan
                        ),
                        "riesz_equation_residual": split.get(
                            "riesz_equation_residual", np.nan
                        ),
                        "target_tangent_norm": split.get("target_support_norm", np.nan),
                        "target_status": row.status,
                    }
                )
    event_frame = pd.DataFrame(events)
    event_keys = (
        event_frame[["dgp", "N", "T", "replication", "target"]]
        .drop_duplicates()
        .assign(any_target_instability_event=True)
        if not event_frame.empty
        else pd.DataFrame(
            columns=[
                "dgp",
                "N",
                "T",
                "replication",
                "target",
                "any_target_instability_event",
            ]
        )
    )
    complete = complete.merge(
        event_keys,
        how="left",
        on=["dgp", "N", "T", "replication", "target"],
    )
    complete["any_target_instability_event"] = complete[
        "any_target_instability_event"
    ].fillna(False)
    summary_rows = []
    for keys, group in complete.groupby(["dgp", "N", "T", "target"], sort=True):
        statuses = group["status"]
        summary_rows.append(
            {
                "dgp": keys[0],
                "N": keys[1],
                "T": keys[2],
                "target": keys[3],
                "requested_replications": len(group),
                "success_rate": float((statuses == "success").mean()),
                "any_failure_rate": float((statuses != "success").mean()),
                "full_riesz_target_instability_rate": float(
                    (statuses == "riesz_target_instability").mean()
                ),
                "split_riesz_target_instability_rate": float(
                    (statuses == "split_riesz_target_instability").mean()
                ),
                "target_unsupported_rate": float(
                    statuses.isin(
                        [
                            "target_unsupported_selected_rank",
                            "split_target_unsupported_selected_rank",
                        ]
                    ).mean()
                ),
                "tangent_gram_failure_rate": float(
                    statuses.isin(
                        [
                            "tangent_gram_eigensolver_failure",
                            "tangent_gram_nearly_singular",
                            "split_tangent_gram_eigensolver_failure",
                            "split_tangent_gram_nearly_singular",
                        ]
                    ).mean()
                ),
                "any_target_instability_event_rate": float(
                    group["any_target_instability_event"].mean()
                ),
                "minimum_tangent_gram_eigenvalue": float(
                    group["tangent_gram_smallest_eigenvalue"].min()
                ),
                "median_tangent_gram_largest_eigenvalue": float(
                    group["tangent_gram_largest_eigenvalue"].median()
                ),
                "median_tangent_gram_condition_number": float(
                    group["tangent_gram_condition_number"].median()
                ),
                "tangent_gram_eigensolver_failure_rate": float(
                    (~group["tangent_gram_eigensolver_converged"].fillna(False).astype(bool)).mean()
                ),
                "rank_selection_failure_rate": float(
                    statuses.isin(
                        ["rank_at_cap", "rank_pilot_failure", "rank_selection_failure"]
                    ).mean()
                ),
                "other_failure_rate": float(
                    (
                        (statuses != "success")
                        & ~statuses.isin(
                            [
                                "riesz_target_instability",
                                "split_riesz_target_instability",
                                "target_unsupported_selected_rank",
                                "split_target_unsupported_selected_rank",
                                "tangent_gram_eigensolver_failure",
                                "tangent_gram_nearly_singular",
                                "split_tangent_gram_eigensolver_failure",
                                "split_tangent_gram_nearly_singular",
                                "rank_at_cap",
                                "rank_pilot_failure",
                                "rank_selection_failure",
                            ]
                        )
                    ).mean()
                ),
            }
        )
    return complete, event_frame, pd.DataFrame(summary_rows)


def write_riesz_diagnostic(
    replications: pd.DataFrame,
    events: pd.DataFrame,
    summary: pd.DataFrame,
    output: Path,
) -> list[Path]:
    """Write Riesz diagnostic data and named-target panel table."""

    from .reporting import TARGET_TITLES

    output.mkdir(parents=True, exist_ok=True)
    paths: list[Path] = []
    for stem, frame in (
        ("riesz_diagnostic_replications", replications),
        ("riesz_target_instability_events", events),
        ("riesz_diagnostic_summary", summary),
    ):
        csv_path = output / f"{stem}.csv"
        parquet_path = output / f"{stem}.parquet"
        frame.to_csv(csv_path, index=False)
        frame.to_parquet(parquet_path, index=False)
        paths.extend([csv_path, parquet_path])
    tex_path = output / "tab_mc_riesz_diagnostics.tex"
    newline = r"\\"
    lines = [
        r"\begingroup",
        r"\small",
        r"\begin{longtable}{rrrrrrrrrrrrrrr}",
        r"\caption{Riesz target-stability and tangent-Gram diagnostics}\label{tab:mc-riesz-diagnostics}" + newline,
        r"\toprule",
        r"DGP & $N$ & $T$ & Success & Failure & Unsupported & Full target & Split target & Gram fail & Any event & Min Gram eig. & Med. cond. & Eig. fail & Rank fail & Other " + newline,
        r"\midrule",
        r"\endfirsthead",
        r"\toprule",
        r"DGP & $N$ & $T$ & Success & Failure & Unsupported & Full target & Split target & Gram fail & Any event & Min Gram eig. & Med. cond. & Eig. fail & Rank fail & Other " + newline,
        r"\midrule",
        r"\endhead",
    ]
    for letter, target in zip(
        "ABCDEFGHIJKLMNOPQRSTUVWXYZ", summary["target"].drop_duplicates(), strict=False
    ):
        title = TARGET_TITLES.get(target, target)
        lines.append(rf"\multicolumn{{15}}{{l}}{{\textit{{Panel {letter}. {title}}}}} " + newline)
        for row in summary[summary["target"] == target].itertuples(index=False):
            values = [
                str(row.dgp), str(row.N), str(row.T), f"{row.success_rate:.3f}",
                f"{row.any_failure_rate:.3f}", f"{row.target_unsupported_rate:.3f}",
                f"{row.full_riesz_target_instability_rate:.3f}",
                f"{row.split_riesz_target_instability_rate:.3f}",
                f"{row.tangent_gram_failure_rate:.3f}",
                f"{row.any_target_instability_event_rate:.3f}",
                f"{row.minimum_tangent_gram_eigenvalue:.3g}",
                f"{row.median_tangent_gram_condition_number:.3g}",
                f"{row.tangent_gram_eigensolver_failure_rate:.3f}",
                f"{row.rank_selection_failure_rate:.3f}", f"{row.other_failure_rate:.3f}",
            ]
            lines.append(" & ".join(values) + " " + newline)
        lines.append(r"\addlinespace")
    lines.extend(
        [
            r"\bottomrule",
            r"\multicolumn{15}{p{0.96\linewidth}}{\footnotesize \textit{Notes:} Unsupported targets are detected before Riesz iteration. Target events use the target-specific Rayleigh quotient; Gram failures use the separately cached matrix-free Lanczos spectrum on nonredundant tangent coordinates. These are numerical diagnostics, not a proof of the maintained global objective-gap or conditioning assumptions.} " + newline,
            r"\end{longtable}",
            r"\endgroup",
        ]
    )
    tex_path.write_text("\n".join(lines) + "\n", encoding="utf-8")
    paths.append(tex_path)
    return paths
