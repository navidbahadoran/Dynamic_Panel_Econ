"""Build the ex-ante calibration candidates without running Monte Carlo replications."""

from __future__ import annotations

import argparse
import math
from pathlib import Path

import pandas as pd

from dynamic_panel_econ.calibration import (
    calibrate_cell,
    calibrate_rank_stress_cell,
)
from dynamic_panel_econ.dgp import DGPParameters

DGPS = (1, 2, 3, 4)
BASELINE_CELLS = ((50, 50), (100, 100), (200, 200), (400, 400))
STRESS_CELLS = ((50, 50), (100, 100), (200, 200))
STRESS_RANKS = ((1, 1, 1), (2, 1, 1), (1, 0, 2))
COMPONENT_STRENGTHS = (1.0, 1.0)


def _row(result, *, rank_stress: bool, ranks: tuple[int, int, int]) -> dict[str, object]:
    diagnostics = result.diagnostics
    return {
        "design": "rank_stress" if rank_stress else "baseline",
        "rank_stress": rank_stress,
        "dgp": result.dgp,
        "n": result.n,
        "t": result.t,
        "true_rank_vector": ",".join(map(str, ranks)),
        "c_h": result.c_h,
        "c_xi": result.c_xi,
        "intended_r2": 0.65,
        "achieved_r2": result.achieved_r2,
        "r2_scale_identified": ranks[1] > 0,
        "pi_h": result.achieved_h_share,
        "realized_calibration_h_share": diagnostics["realized_calibration_h_share"],
        "population_var_u_tilde": diagnostics["population_var_u_tilde"],
        "population_var_h_raw": diagnostics["population_var_h_raw"],
        "C_A": diagnostics["theoretical_max_abs_A"],
        "C_beta": diagnostics["theoretical_max_abs_B"],
        "C_H": diagnostics["theoretical_max_abs_H"],
        "C_Theta": diagnostics["theoretical_coefficient_envelope"],
        "calibration_seed": 20260807,
        "calibration_draws": diagnostics["calibration_draws"],
    }


def build(draws: int) -> pd.DataFrame:
    params = DGPParameters()
    rows: list[dict[str, object]] = []
    for dgp in DGPS:
        for n, t in BASELINE_CELLS:
            result = calibrate_cell(
                dgp, n, t, 20260807, params=params, target_r2=0.65, draws=draws
            )
            rows.append(_row(result, rank_stress=False, ranks=(1, 1, 1)))
        for n, t in STRESS_CELLS:
            for ranks in STRESS_RANKS:
                result = calibrate_rank_stress_cell(
                    dgp,
                    n,
                    t,
                    ranks,
                    20260807,
                    params=params,
                    target_r2=0.65,
                    draws=draws,
                    component_strengths=COMPONENT_STRENGTHS,
                )
                rows.append(_row(result, rank_stress=True, ranks=ranks))
    return pd.DataFrame(rows)


def _toml_bool(value: bool) -> str:
    return "true" if value else "false"


def write_toml(frame: pd.DataFrame, destination: Path, draws: int) -> None:
    lines = [
        "schema_version = 1",
        'status = "candidate_not_activated"',
        'method = "population_c_h_and_independent_frozen_cell_c_xi"',
        "calibration_seed = 20260807",
        f"calibration_draws = {draws}",
        "",
    ]
    for row in frame.itertuples(index=False):
        ranks = [int(value) for value in row.true_rank_vector.split(",")]
        lines.extend(
            [
                "[[calibration]]",
                f'rank_stress = {_toml_bool(bool(row.rank_stress))}',
                f'design = "{row.design}"',
                f"dgp = {int(row.dgp)}",
                f"n = {int(row.n)}",
                f"t = {int(row.t)}",
                f"true_rank_vector = {ranks}",
                f"c_h = {float(row.c_h):.17g}",
                f"c_xi = {float(row.c_xi):.17g}",
                "intended_r2 = 0.65",
                f"achieved_r2 = {float(row.achieved_r2):.17g}",
                f'r2_scale_identified = {_toml_bool(bool(row.r2_scale_identified))}',
                f"pi_h = {float(row.pi_h):.17g}",
                f"population_var_u_tilde = {float(row.population_var_u_tilde):.17g}",
                f"population_var_h_raw = {float(row.population_var_h_raw):.17g}",
                f"C_A = {float(row.C_A):.17g}",
                f"C_beta = {float(row.C_beta):.17g}",
                f"C_H = {float(row.C_H):.17g}",
                f"C_Theta = {float(row.C_Theta):.17g}",
                "",
            ]
        )
    destination.write_text("\n".join(lines), encoding="utf-8")


def write_outputs(frame: pd.DataFrame, output: Path) -> None:
    output.mkdir(parents=True, exist_ok=True)
    frame.to_csv(output / "frozen_calibration_candidates.csv", index=False)
    frame[
        [
            "design",
            "dgp",
            "n",
            "t",
            "true_rank_vector",
            "population_var_u_tilde",
            "population_var_h_raw",
            "c_h",
            "pi_h",
            "realized_calibration_h_share",
        ]
    ].to_csv(output / "population_variance_components.csv", index=False)
    frame[
        [
            "design",
            "dgp",
            "n",
            "t",
            "true_rank_vector",
            "C_A",
            "C_beta",
            "C_H",
            "C_Theta",
        ]
    ].to_csv(output / "theoretical_coefficient_envelopes.csv", index=False)
    maximum = float(frame["C_Theta"].max())
    retained = maximum <= 8.0
    proposed_b = 9.0 if retained else float(math.ceil(maximum) + 1)
    pd.DataFrame(
        [
            {
                "C_Theta_max": maximum,
                "current_B": 9.0,
                "c_B": 1.0,
                "current_B_minus_c_B": 8.0,
                "B9_verified": retained,
                "proposed_common_B": proposed_b,
                "proposed_c_B": 1.0,
                "proposed_interior_margin": proposed_b - maximum,
                "status": "CANDIDATE_NOT_ACTIVATED",
            }
        ]
    ).to_csv(output / "proposed_B_summary.csv", index=False)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--draws", type=int, default=50)
    parser.add_argument(
        "--config-output", default="configs/mc/frozen_dgp_calibration.toml"
    )
    parser.add_argument(
        "--audit-output", default="results/mc/audit/dgp_theorem_alignment"
    )
    args = parser.parse_args()
    frame = build(args.draws)
    config_output = Path(args.config_output)
    config_output.parent.mkdir(parents=True, exist_ok=True)
    write_toml(frame, config_output, args.draws)
    write_outputs(frame, Path(args.audit_output))


if __name__ == "__main__":
    main()
