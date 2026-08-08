"""Evaluate pooled-R2 feasibility without changing the production configuration."""

from __future__ import annotations

import argparse
from pathlib import Path

from dynamic_panel_econ.config import load_config
from dynamic_panel_econ.diagnostics import r2_feasibility_rows, write_r2_feasibility


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--config", default="configs/mc/production.toml")
    parser.add_argument("--output", default="results/mc/calibration")
    args = parser.parse_args()
    frame = r2_feasibility_rows(load_config(args.config))
    paths = write_r2_feasibility(frame, Path(args.output))
    print(frame.to_string(index=False))
    max_floor = float(frame["large_c_xi_r2_floor"].max())
    candidates = [
        target
        for target in (0.60, 0.65, 0.70)
        if target >= max_floor + 0.04
        and frame[f"target_{target:.2f}_feasible".replace(".", "_")].all()
    ]
    print(f"\nmax_floor={max_floor:.10f}")
    print(f"recommended_common_target={min(candidates) if candidates else 'none'}")
    print("\n".join(str(path) for path in paths))


if __name__ == "__main__":
    main()
