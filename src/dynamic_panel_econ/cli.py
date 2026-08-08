"""Comprehensive Monte Carlo CLI with CLI > TOML > defaults precedence."""

from __future__ import annotations

import argparse
import json
from copy import deepcopy
from pathlib import Path
from typing import Any

from .config import DEFAULTS, load_config, validate_config


def _csv_ints(value: str) -> list[int]:
    try:
        return [int(item.strip()) for item in value.split(",") if item.strip()]
    except ValueError as exc:
        raise argparse.ArgumentTypeError("expected comma-separated integers") from exc


def _csv_floats(value: str) -> list[float]:
    try:
        return [float(item.strip()) for item in value.split(",") if item.strip()]
    except ValueError as exc:
        raise argparse.ArgumentTypeError("expected comma-separated numbers") from exc


def _csv_strings(value: str) -> list[str]:
    return [item.strip() for item in value.split(",") if item.strip()]


def build_run_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description=(
            "Run the Revision-8 Monte Carlo. --pooled-r2-target calibrates the simulated DGP "
            "through c_xi; it is not an estimation or rank-selection tuning parameter."
        )
    )
    parser.add_argument("--config", help="Optional TOML configuration; CLI values override it.")
    simulation = parser.add_argument_group("simulation")
    simulation.add_argument("--dgp", type=int, choices=(1, 2, 3, 4))
    simulation.add_argument("--dgp-grid", type=_csv_ints)
    simulation.add_argument("--N", type=int, dest="panel_n")
    simulation.add_argument("--T", type=int, dest="panel_t")
    simulation.add_argument("--balanced-grid", type=_csv_ints)
    simulation.add_argument("--replications", type=int)
    simulation.add_argument("--burn-in", type=int)
    simulation.add_argument("--seed", type=int)

    dgp = parser.add_argument_group("DGP")
    for flag in (
        "rho-g", "mu-f-a", "kappa-f-a", "mu-f-b", "kappa-f-b", "rho-x",
        "delta-x", "rho-fx", "eta-x", "rho-s", "pi-h", "stability-bound",
        "mu-lambda-a-1", "mu-lambda-a-2", "sigma-lambda-a",
        "mu-lambda-b-1", "mu-lambda-b-2", "sigma-lambda-b",
    ):
        dgp.add_argument(f"--{flag}", type=float)
    calibration = parser.add_argument_group("calibration")
    calibration.add_argument(
        "--pooled-r2-target",
        type=float,
        help=(
            "Requested pooled DGP R-squared used to calibrate c_xi (baseline 0.65). "
            "For rank-stress r_B=0, c_xi=1 and induced R-squared is reported instead."
        ),
    )
    calibration.add_argument("--calibration-draws", type=int)
    calibration.add_argument("--calibration-seed", type=int)
    calibration.add_argument("--calibration-tolerance", type=float)

    rank = parser.add_argument_group("rank")
    rank.add_argument("--rank-mode", choices=("fixed", "selected"))
    rank.add_argument("--fixed-ranks", type=_csv_ints)
    rank.add_argument("--rank-caps", type=_csv_ints)
    rank.add_argument("--true-ranks", type=_csv_ints)
    rank.add_argument("--gamma", type=float)
    rank.add_argument("--epsilon", type=float)
    rank.add_argument("--threshold-multiplier", type=float)
    rank.add_argument("--eta-for-penalty", type=float)
    rank.add_argument("--spatial-dimension", type=int)
    rank.add_argument("--ic-multiplier", type=float)

    optimization = parser.add_argument_group("optimization")
    optimization.add_argument("--coefficient-bound", type=float)
    optimization.add_argument("--simulation-interior-margin", type=float)
    optimization.add_argument("--max-als-iterations", type=int)
    optimization.add_argument("--objective-tolerance", type=float)
    optimization.add_argument("--stationarity-tolerance", type=float)
    optimization.add_argument("--start-objective-stability-tol", type=float)
    optimization.add_argument("--max-cap-pilot-routes", type=int)
    optimization.add_argument("--start-envelope-fraction", type=float)

    riesz = parser.add_argument_group("Riesz and variance")
    riesz.add_argument("--riesz-solver", choices=("auto", "cg", "minres", "lsmr"))
    riesz.add_argument("--riesz-tolerance", type=float)
    riesz.add_argument("--riesz-max-iterations", type=int)
    riesz.add_argument("--target-support-tolerance", type=float)
    riesz.add_argument("--tangent-gram-min-eigenvalue-floor", type=float)
    riesz.add_argument("--variance-type", choices=("auto", "diagonal", "spatial"))
    riesz.add_argument("--c-sp", type=float)

    execution = parser.add_argument_group("execution and output")
    execution.add_argument("--workers", "--n-jobs", dest="workers", type=int)
    execution.add_argument("--parallel-level", choices=("replications", "none"))
    execution.add_argument("--chunk-size", type=int)
    modes = execution.add_mutually_exclusive_group()
    modes.add_argument("--resume", action="store_true")
    modes.add_argument("--overwrite", action="store_true")
    execution.add_argument("--output-root")
    execution.add_argument("--run-name")
    execution.add_argument("--targets", type=_csv_strings)
    execution.add_argument(
        "--save-candidate-details", action=argparse.BooleanOptionalAction, default=None
    )
    execution.add_argument("--print-resolved-config", action="store_true")
    execution.add_argument("--dry-run", action="store_true")
    execution.add_argument("--experiment", choices=("baseline", "power"))
    execution.add_argument("--alternative-grid", type=_csv_floats)
    execution.add_argument("--power-block", choices=("A", "B"))
    return parser


_DGP_NAMES = (
    "rho_g", "mu_f_a", "kappa_f_a", "mu_f_b", "kappa_f_b", "rho_x",
    "delta_x", "rho_fx", "eta_x", "rho_s", "pi_h", "stability_bound",
    "mu_lambda_a_1", "mu_lambda_a_2", "sigma_lambda_a",
    "mu_lambda_b_1", "mu_lambda_b_2", "sigma_lambda_b",
    "burn_in", "calibration_draws", "calibration_seed", "calibration_tolerance",
)


def _base_config(path: str | None) -> dict[str, Any]:
    return load_config(path, validate=False) if path else deepcopy(DEFAULTS)


def resolve_run_args(args: argparse.Namespace) -> dict[str, Any]:
    config = _base_config(args.config)
    run, dgp, estimation, inference = (
        config["run"], config["dgp"], config["estimation"], config["inference"]
    )
    values = vars(args)
    if args.dgp is not None and args.dgp_grid is not None:
        raise ValueError("use either --dgp or --dgp-grid, not both")
    if args.dgp is not None:
        run["dgps"] = [args.dgp]
    elif args.dgp_grid is not None:
        run["dgps"] = args.dgp_grid
    if args.balanced_grid is not None:
        run["cells"] = [[size, size] for size in args.balanced_grid]
    elif args.panel_n is not None or args.panel_t is not None:
        if args.panel_n is None or args.panel_t is None:
            raise ValueError("--N and --T must be supplied together")
        run["cells"] = [[args.panel_n, args.panel_t]]
    direct_run = {
        "replications": "replications", "seed": "master_seed", "workers": "n_jobs",
        "parallel_level": "parallel_level", "chunk_size": "chunk_size",
        "output_root": "output_root", "run_name": "name", "rank_mode": "rank_mode",
        "experiment": "experiment", "alternative_grid": "alternative_grid",
        "power_block": "power_block", "save_candidate_details": "save_candidate_details",
    }
    for source, destination in direct_run.items():
        if values.get(source) is not None:
            run[destination] = values[source]
    for name in _DGP_NAMES:
        if values.get(name) is not None:
            dgp[name] = values[name]
    if args.pooled_r2_target is not None:
        dgp["target_r2"] = args.pooled_r2_target
    estimation_map = {
        "fixed_ranks": "fixed_ranks", "rank_caps": "rank_caps",
        "gamma": "nuclear_gamma", "epsilon": "nuclear_epsilon",
        "threshold_multiplier": "threshold_multiplier",
        "eta_for_penalty": "eta_for_penalty", "spatial_dimension": "spatial_dimension",
        "ic_multiplier": "ic_multiplier", "coefficient_bound": "coefficient_bound",
        "simulation_interior_margin": "simulation_interior_margin",
        "max_als_iterations": "max_sweeps", "objective_tolerance": "objective_rtol",
        "stationarity_tolerance": "stationarity_tol",
        "start_objective_stability_tol": "start_objective_stability_tol",
        "max_cap_pilot_routes": "rank_adaptive_max_routes",
        "start_envelope_fraction": "cap_pilot_start_envelope_fraction",
    }
    for source, destination in estimation_map.items():
        if values.get(source) is not None:
            estimation[destination] = values[source]
    if args.true_ranks is not None:
        config.setdefault("rank_stress", {})["true_rank_vectors"] = [args.true_ranks]
    inference_map = {
        "riesz_solver": "riesz_solver", "riesz_tolerance": "riesz_tol",
        "riesz_max_iterations": "riesz_max_iter",
        "target_support_tolerance": "target_support_tolerance",
        "tangent_gram_min_eigenvalue_floor": "tangent_gram_min_eigenvalue_floor",
        "variance_type": "variance_type", "c_sp": "spatial_c", "targets": "targets",
    }
    for source, destination in inference_map.items():
        if values.get(source) is not None:
            inference[destination] = values[source]
    if dgp.get("calibration_seed") is None:
        dgp["calibration_seed"] = int(run["master_seed"])
    validate_config(config)
    return config


def resolved_config_text(config: dict[str, Any]) -> str:
    """Stable, valid TOML serialization for researcher inspection."""

    lines = ["# Fully resolved configuration; generated by scripts/run_mc.py."]

    def emit(prefix: str, values: dict[str, Any]) -> None:
        lines.append(f"\n[{prefix}]")
        nested = []
        for key, value in values.items():
            if isinstance(value, dict):
                nested.append((key, value))
            elif value is not None:
                lines.append(f"{key} = {json.dumps(value, ensure_ascii=True)}")
        for key, value in nested:
            emit(f"{prefix}.{key}", value)

    for section, values in config.items():
        emit(section, values)
    return "\n".join(lines) + "\n"


def write_cli_help(path: str | Path) -> None:
    Path(path).write_text(build_run_parser().format_help(), encoding="utf-8")
