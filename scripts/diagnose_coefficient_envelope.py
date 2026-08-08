"""Calibrate every actual rank-stress design and derive one deterministic common B."""

from __future__ import annotations

import argparse
import json
from dataclasses import asdict
from pathlib import Path

import numpy as np
import pandas as pd

from dynamic_panel_econ.calibration import calibrate_rank_stress_cell
from dynamic_panel_econ.config import load_config
from dynamic_panel_econ.monte_carlo import _params


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--config", default="configs/mc/rank_pilot_diagnostic.toml")
    parser.add_argument("--output", default="results/mc/calibration")
    parser.add_argument("--margin", type=float, default=1.0)
    args = parser.parse_args()
    config = load_config(args.config)
    rows = []
    params = _params(config)
    for dgp in config["run"]["dgps"]:
        for n, t in config["run"]["cells"]:
            for vector in config["rank_stress"]["true_rank_vectors"]:
                true_rank = tuple(int(value) for value in vector)
                result = calibrate_rank_stress_cell(
                    int(dgp),
                    int(n),
                    int(t),
                    true_rank,
                    int(config["run"]["master_seed"]),
                    component_strengths=tuple(config["rank_stress"]["component_strengths"]),
                    params=params,
                    pi_h=float(config["dgp"]["pi_h"]),
                    target_r2=float(config["dgp"]["target_r2"]),
                    draws=int(config["dgp"]["calibration_draws"]),
                    allow_infeasible_diagnostic=True,
                )
                row = asdict(result)
                row.update(row.pop("diagnostics"))
                row["true_rank_vector"] = json.dumps(true_rank)
                rows.append(row)
    frame = pd.DataFrame(rows).sort_values(["dgp", "n", "t", "true_rank_vector"])
    maximum = float(frame["theoretical_coefficient_envelope"].max())
    common_bound = float(np.ceil(maximum + args.margin))
    frame["recommended_common_B"] = common_bound
    frame["required_simulation_interior_margin"] = args.margin
    frame["deterministic_interior_margin"] = common_bound - frame[
        "theoretical_coefficient_envelope"
    ]
    frame["envelope_condition_pass"] = (
        frame["theoretical_coefficient_envelope"] <= common_bound - args.margin
    )
    output = Path(args.output)
    output.mkdir(parents=True, exist_ok=True)
    csv_path = output / "rank_stress_calibration.csv"
    parquet_path = output / "rank_stress_calibration.parquet"
    tex_path = output / "tab_mc_rank_stress_calibration.tex"
    frame.to_csv(csv_path, index=False)
    frame.to_parquet(parquet_path, index=False)
    newline = r"\\"
    lines = [
        r"\begingroup",
        r"\small",
        r"\begin{longtable}{rrlrrrrrrrrr}",
        r"\caption{Rank-stress calibration and coefficient envelopes}\label{tab:mc-rank-stress-calibration}" + newline,
        r"\toprule",
        r"DGP & $N=T$ & True rank & $c_H$ & $c_\xi$ & $\pi_H$ & $R^2$ & Feasible & Env. A & Env. B & Env. H & Margin " + newline,
        r"\midrule",
        r"\endfirsthead",
        r"\toprule",
        r"DGP & $N=T$ & True rank & $c_H$ & $c_\xi$ & $\pi_H$ & $R^2$ & Feasible & Env. A & Env. B & Env. H & Margin " + newline,
        r"\midrule",
        r"\endhead",
    ]
    for row in frame.itertuples(index=False):
        values = [
            str(row.dgp), str(row.n), str(row.true_rank_vector), f"{row.c_h:.4f}",
            f"{row.c_xi:.4f}", f"{row.achieved_h_share:.4f}", f"{row.achieved_r2:.4f}",
            "yes" if row.target_r2_feasible else "no",
            f"{row.theoretical_max_abs_A:.4f}", f"{row.theoretical_max_abs_B:.4f}",
            f"{row.theoretical_max_abs_H:.4f}", f"{row.deterministic_interior_margin:.4f}",
        ]
        lines.append(" & ".join(values) + " " + newline)
    lines.extend(
        [
            r"\bottomrule",
            r"\multicolumn{12}{p{0.96\linewidth}}{\footnotesize \textit{Notes:} In infeasible rows, $c_\xi=1$ is a reporting normalization, not a successful calibration. The strict Monte Carlo runner rejects those cells.}",
            r"\end{longtable}",
            r"\endgroup",
        ]
    )
    tex_path.write_text("\n".join(lines) + "\n", encoding="utf-8")
    print(frame.to_string(index=False))
    print(f"maximum_deterministic_envelope={maximum:.10f}")
    print(f"recommended_common_B={common_bound:.1f}")
    print(f"required_margin={args.margin:.1f}")
    print(csv_path, parquet_path, tex_path, sep="\n")


if __name__ == "__main__":
    main()
