"""Run deterministic, resume-safe Monte Carlo chunks."""

from __future__ import annotations

import sys

from dynamic_panel_econ.cli import build_run_parser, resolve_run_args, resolved_config_text
from dynamic_panel_econ.mc_accounting import power_design_configs
from dynamic_panel_econ.monte_carlo import run_monte_carlo


def main() -> None:
    parser = build_run_parser()
    args = parser.parse_args()
    config = resolve_run_args(args)
    if args.print_resolved_config or args.dry_run:
        print(resolved_config_text(config), end="")
    if args.dry_run:
        print("Dry run only: no calibration, fitting, inference, or output was executed.")
        return
    if config["run"]["experiment"] == "power" and args.ic_multiplier is None:
        parser.error("--experiment power requires an explicitly approved --ic-multiplier")
    designs = power_design_configs(config) if config["run"]["experiment"] == "power" else [config]
    for design in designs:
        if design["run"]["experiment"] == "power":
            delta = design["run"]["nominal_delta"]
            mode = design["run"]["rank_mode"]
            design["run"]["name"] = f"{design['run']['name']}_power_{design['run']['power_block']}_{delta:g}_{mode}"
        _, root = run_monte_carlo(
            design,
            resume=args.resume,
            overwrite=args.overwrite,
            n_jobs=args.workers,
            cli_argv=sys.argv,
        )
        print(f"Completed Monte Carlo chunks under {root}")


if __name__ == "__main__":
    main()
