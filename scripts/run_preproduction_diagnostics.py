"""Run the long, resume-safe diagnostic pilots without launching production."""

from __future__ import annotations

import argparse
import csv
import subprocess
import sys
from pathlib import Path


def run(*arguments: str) -> None:
    command = [sys.executable, *arguments]
    print("RUN", subprocess.list2cmdline(command), flush=True)
    subprocess.run(command, check=True)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--execute-long-diagnostics",
        action="store_true",
        help="Explicitly permit long runs after all deterministic gates pass.",
    )
    args = parser.parse_args()
    # The feasibility study is quick and deterministic; refresh it before the
    # long pilots so the three top-level artifacts always describe this run.
    run("scripts/diagnose_r2_feasibility.py", "--config", "configs/mc/production.toml")
    run(
        "scripts/diagnose_coefficient_envelope.py",
        "--config",
        "configs/mc/rank_pilot_diagnostic.toml",
    )
    calibration_path = Path("results/mc/calibration/rank_stress_calibration.csv")
    with calibration_path.open(newline="", encoding="utf-8") as handle:
        infeasible = [
            row
            for row in csv.DictReader(handle)
            if row["target_r2_feasible"].lower() != "true"
        ]
    if infeasible:
        vectors = sorted({row["true_rank_vector"] for row in infeasible})
        raise RuntimeError(
            "long diagnostics remain blocked: common pooled-R2 calibration is infeasible for "
            f"true-rank vectors {vectors}"
        )
    if not args.execute_long_diagnostics:
        raise RuntimeError(
            "deterministic gates passed, but long diagnostics require "
            "--execute-long-diagnostics"
        )
    run(
        "scripts/run_mc.py",
        "--config",
        "configs/mc/rank_pilot_diagnostic.toml",
        "--n-jobs",
        "8",
        "--resume",
    )
    run("scripts/aggregate_rank_pilot.py", "--config", "configs/mc/rank_pilot_diagnostic.toml")
    run(
        "scripts/run_mc.py",
        "--config",
        "configs/mc/riesz_diagnostic.toml",
        "--n-jobs",
        "8",
        "--resume",
    )
    run(
        "scripts/aggregate_riesz_diagnostics.py",
        "--config",
        "configs/mc/riesz_diagnostic.toml",
    )
    marker = Path("results/mc/diagnostics/preproduction_diagnostics_complete.txt")
    marker.parent.mkdir(parents=True, exist_ok=True)
    marker.write_text("complete\n", encoding="utf-8")
    print(f"COMPLETE {marker}", flush=True)


if __name__ == "__main__":
    main()
