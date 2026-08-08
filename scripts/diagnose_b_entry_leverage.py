"""Compare the maintained and proposed slope-factor leverage specifications."""

from __future__ import annotations

import argparse
from dataclasses import asdict, replace
from pathlib import Path

import numpy as np
import pandas as pd

from dynamic_panel_econ.calibration import calibrate_cell
from dynamic_panel_econ.config import load_config
from dynamic_panel_econ.dgp import bounded_ar_envelope, generate_panel
from dynamic_panel_econ.monte_carlo import _params
from dynamic_panel_econ.seeds import seed_sequence


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--config", default="configs/mc/rank_pilot_diagnostic.toml")
    parser.add_argument("--output", default="results/mc/calibration")
    parser.add_argument("--draws", type=int, default=None)
    args = parser.parse_args()
    config = load_config(args.config)
    draws = int(args.draws or config["dgp"]["calibration_draws"])
    base = _params(config)
    rows = []
    for kappa in (0.20, 0.15):
        params = replace(base, kappa_f_b=kappa)
        factor_lower = params.mu_f_b - abs(kappa) * bounded_ar_envelope(params.rho_g)
        factor_upper = params.mu_f_b + abs(kappa) * bounded_ar_envelope(params.rho_g)
        for dgp in config["run"]["dgps"]:
            for n, t in config["run"]["cells"]:
                calibration = calibrate_cell(
                    int(dgp),
                    int(n),
                    int(t),
                    int(config["run"]["master_seed"]),
                    params=params,
                    pi_h=float(config["dgp"]["pi_h"]),
                    target_r2=0.65,
                    draws=draws,
                )
                record = asdict(calibration)
                diagnostics = record.pop("diagnostics")
                record.update(
                    {
                        "kappa_f_b": kappa,
                        "theoretical_factor_b_lower": factor_lower,
                        "theoretical_factor_b_upper": factor_upper,
                        "mean_b": diagnostics["mean_b"],
                        "sd_b": diagnostics["sd_b"],
                        "min_b": diagnostics["min_b"],
                        "max_b": diagnostics["max_b"],
                        "theoretical_coefficient_envelope": diagnostics[
                            "theoretical_coefficient_envelope"
                        ],
                        "B_G1_fixed_time_true": np.nan,
                        "B_G2_fixed_time_true": np.nan,
                        "B_G2_minus_G1_fixed_time_true": np.nan,
                        "B_G1_time_average_true": np.nan,
                        "B_G2_time_average_true": np.nan,
                        "B_G2_minus_G1_time_average_true": np.nan,
                    }
                )
                if int(dgp) == 4:
                    truths = []
                    for replication in range(draws):
                        panel = generate_panel(
                            4,
                            int(n),
                            int(t),
                            seed_sequence(
                                config["run"]["master_seed"],
                                "b_entry_comparison",
                                kappa,
                                n,
                                t,
                                replication,
                            ),
                            c_h=calibration.c_h,
                            c_xi=calibration.c_xi,
                            params=params,
                            coefficient_bound=9.0,
                            simulation_interior_margin=1.0,
                        )
                        truths.append(panel.truths)
                    for name in (
                        "B_G1_fixed_time_true",
                        "B_G2_fixed_time_true",
                        "B_G2_minus_G1_fixed_time_true",
                        "B_G1_time_average_true",
                        "B_G2_time_average_true",
                        "B_G2_minus_G1_time_average_true",
                    ):
                        record[name] = float(np.mean([truth[name] for truth in truths]))
                rows.append(record)
    frame = pd.DataFrame(rows).sort_values(["dgp", "n", "t", "kappa_f_b"])
    output = Path(args.output)
    output.mkdir(parents=True, exist_ok=True)
    csv_path = output / "b_entry_leverage_comparison.csv"
    parquet_path = output / "b_entry_leverage_comparison.parquet"
    tex_path = output / "tab_mc_b_entry_leverage_comparison.tex"
    frame.to_csv(csv_path, index=False)
    frame.to_parquet(parquet_path, index=False)
    newline = r"\\"
    lines = [
        r"\begingroup",
        r"\small",
        r"\begin{longtable}{rrrrrrrrrrrr}",
        r"\caption{Deterministic B-entry leverage comparison}\label{tab:mc-b-entry-leverage}" + newline,
        r"\toprule",
        r"DGP & $N=T$ & $\kappa_{f,b}$ & $\inf f_b$ & Mean B & sd B & Min B & Max B & $c_H$ & $c_\xi$ & $R^2$ & Envelope " + newline,
        r"\midrule",
        r"\endfirsthead",
        r"\toprule",
        r"DGP & $N=T$ & $\kappa_{f,b}$ & $\inf f_b$ & Mean B & sd B & Min B & Max B & $c_H$ & $c_\xi$ & $R^2$ & Envelope " + newline,
        r"\midrule",
        r"\endhead",
    ]
    for row in frame.itertuples(index=False):
        values = [
            str(row.dgp),
            str(row.n),
            f"{row.kappa_f_b:.2f}",
            f"{row.theoretical_factor_b_lower:.3f}",
            f"{row.mean_b:.4f}",
            f"{row.sd_b:.4f}",
            f"{row.min_b:.4f}",
            f"{row.max_b:.4f}",
            f"{row.c_h:.4f}",
            f"{row.c_xi:.4f}",
            f"{row.achieved_r2:.4f}",
            f"{row.theoretical_coefficient_envelope:.4f}",
        ]
        lines.append(" & ".join(values) + " " + newline)
    lines.extend(
        [
            r"\bottomrule",
            r"\multicolumn{12}{p{0.96\linewidth}}{\footnotesize \textit{Notes:} The maintained design uses $\kappa_{f,b}=0.20$ and has deterministic lower support zero. The proposed comparison uses 0.15 and has lower support 0.15. No production DGP is changed by this table.} " + newline,
            r"\end{longtable}",
            r"\endgroup",
        ]
    )
    tex_path.write_text("\n".join(lines) + "\n", encoding="utf-8")
    print(frame.to_string(index=False))
    print(csv_path, parquet_path, tex_path, sep="\n")


if __name__ == "__main__":
    main()
