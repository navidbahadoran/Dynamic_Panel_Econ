"""Run deterministic, resume-safe Monte Carlo chunks."""

from __future__ import annotations

import argparse

from dynamic_panel_econ.config import load_config
from dynamic_panel_econ.monte_carlo import run_monte_carlo


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--config", required=True)
    parser.add_argument("--n-jobs", type=int)
    modes = parser.add_mutually_exclusive_group()
    modes.add_argument("--resume", action="store_true")
    modes.add_argument("--overwrite", action="store_true")
    args = parser.parse_args()
    config = load_config(args.config)
    _, root = run_monte_carlo(
        config,
        resume=args.resume,
        overwrite=args.overwrite,
        n_jobs=args.n_jobs,
    )
    print(f"Completed Monte Carlo chunks under {root}")


if __name__ == "__main__":
    main()
