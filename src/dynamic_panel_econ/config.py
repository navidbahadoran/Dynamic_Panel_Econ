"""TOML configuration loading, validation, resolution, and hashing."""

from __future__ import annotations

import hashlib
import json
import tomllib
from copy import deepcopy
from pathlib import Path
from typing import Any

DEFAULTS: dict[str, Any] = {
    "run": {
        "name": "monte_carlo",
        "experiment": "baseline",
        "rank_mode": "selected",
        "master_seed": 20260807,
        "dgps": [1, 2, 3, 4],
        "cells": [[50, 50]],
        "replications": 1000,
        "replication_start": 0,
        "chunk_size": 20,
        "parallel_level": "replications",
        "n_jobs": 1,
        "blas_threads": 1,
        "output_root": "results/mc",
        "save_candidate_details": True,
        "alternative_grid": [0.0],
        "power_block": "A",
    },
    "dgp": {
        "burn_in": 50,
        "rho_g": 0.5,
        "rho_s": 0.5,
        "rho_x": 0.5,
        "delta_x": 0.5,
        "eta_x": 0.3,
        "mu_f_a": 0.5,
        "kappa_f_a": 0.1,
        "mu_f_b": 0.6,
        "kappa_f_b": 0.2,
        "pi_h": 0.30,
        "target_r2": 0.65,
        "mu_lambda_a_1": 0.9,
        "mu_lambda_a_2": 1.1,
        "sigma_lambda_a": 0.08,
        "mu_lambda_b_1": 0.8,
        "mu_lambda_b_2": 1.2,
        "sigma_lambda_b": 0.25,
        "stability_bound": 0.85,
        "auto_adjust_group_gap": False,
        "min_abs_postscale_group_diff_a": 0.05,
        "prespecified_gap_grid": [0.2, 0.3, 0.4, 0.5],
        "calibration_draws": 3,
        "calibration_seed": None,
        "frozen_calibration_path": None,
        "calibration_tolerance": 1e-10,
        "rho_fx": 0.0,
    },
    "estimation": {
        "fixed_ranks": [1, 1, 1],
        "rank_caps": [3, 3, 3],
        "coefficient_bound": 9.0,
        "simulation_interior_margin": 1.0,
        "max_sweeps": 200,
        "objective_rtol": 1e-8,
        "stationarity_tol": 1e-6,
        "lstsq_rcond": 1e-10,
        "nuclear_gamma": 0.8,
        "nuclear_epsilon": 0.01,
        "nuclear_max_iter": 500,
        "nuclear_tol": 1e-7,
        "dykstra_max_iter": 100,
        "dykstra_tol": 1e-9,
        "eta_for_penalty": 4.0,
        "spatial_dimension": 1,
        "ic_multiplier": 1.0,
        "threshold_multiplier": 1.0,
        "start_objective_stability_tol": 1e-6,
        "rank_adaptive_improvement_tol": 1e-7,
        "rank_adaptive_removal_tol": 1e-7,
        "rank_adaptive_max_steps": 12,
        "rank_adaptive_max_routes": 6,
        "cap_pilot_start_envelope_fraction": 0.8,
        "dense_nuclear_gamma": 0.8944271909999159,
        "threshold_sensitivity_multipliers": [0.5, 1.0, 2.0],
        "ic_sensitivity_multipliers": [0.5, 1.0, 2.0],
        "larger_rank_caps": [4, 4, 4],
        "compute_rank_sensitivities": True,
        "compute_dense_grid_sensitivity": True,
        "compute_larger_cap_sensitivity": True,
    },
    "inference": {
        "riesz_solver": "auto",
        "riesz_tol": 1e-8,
        "riesz_max_iter": 1000,
        "riesz_target_rayleigh_floor": 1e-12,
        "split_relative_rank_floor": 1e-10,
        "compute_tangent_gram": False,
        "tangent_gram_tol": 1e-5,
        "tangent_gram_max_iter": 500,
        "tangent_gram_min_eigenvalue_floor": 1e-10,
        "target_support_tolerance": 1e-12,
        "spatial_c": 1.0,
        "variance_type": "auto",
        "targets": ["A_entry", "B_entry", "A_fixed_time_mean", "A_full_mean"],
    },
}


def _merge(base: dict[str, Any], update: dict[str, Any]) -> dict[str, Any]:
    out = deepcopy(base)
    for key, value in update.items():
        if isinstance(value, dict) and isinstance(out.get(key), dict):
            out[key] = _merge(out[key], value)
        else:
            out[key] = value
    return out


def load_config(path: str | Path, *, validate: bool = True) -> dict[str, Any]:
    with Path(path).open("rb") as handle:
        supplied = tomllib.load(handle)
    resolved = _merge(DEFAULTS, supplied)
    if validate:
        validate_config(resolved)
    return resolved


def validate_config(config: dict[str, Any]) -> None:
    run = config["run"]
    if run["parallel_level"] not in {"replications", "none"}:
        raise ValueError("parallel_level must be replications or none")
    if run["rank_mode"] not in {"fixed", "selected"}:
        raise ValueError("rank_mode must be fixed or selected")
    if run["experiment"] not in {"baseline", "power"}:
        raise ValueError("experiment must be baseline or power")
    if int(run["replications"]) < 1 or int(run["chunk_size"]) < 1:
        raise ValueError("replications and chunk_size must be positive")
    if int(run["replication_start"]) < 0:
        raise ValueError("replication_start must be nonnegative")
    for n, t in run["cells"]:
        if int(n) < 2 or int(t) < 2:
            raise ValueError("panel dimensions must be at least two")
    if any(int(dgp) not in {1, 2, 3, 4} for dgp in run["dgps"]):
        raise ValueError("DGP identifiers must be 1, 2, 3, or 4")
    if len(config["estimation"]["rank_caps"]) != 3:
        raise ValueError("baseline P=K=1 configuration requires three rank caps")
    if len(config["estimation"]["fixed_ranks"]) != 3:
        raise ValueError("baseline P=K=1 configuration requires three fixed ranks")
    estimation = config["estimation"]
    if run["rank_mode"] == "selected":
        multiplier = estimation["ic_multiplier"]
        if isinstance(multiplier, str) or float(multiplier) <= 0.0:
            raise ValueError("selected rank mode requires a fixed positive ic_multiplier")
    if float(estimation["coefficient_bound"]) <= float(
        estimation["simulation_interior_margin"]
    ):
        raise ValueError("coefficient_bound must exceed simulation_interior_margin")
    if float(estimation["start_objective_stability_tol"]) <= 0.0:
        raise ValueError("start_objective_stability_tol must be positive")
    if int(estimation["rank_adaptive_max_routes"]) < 2:
        raise ValueError("rank_adaptive_max_routes must be at least two")
    start_fraction = float(estimation["cap_pilot_start_envelope_fraction"])
    if not 0.0 < start_fraction < 1.0:
        raise ValueError("cap_pilot_start_envelope_fraction must be strictly between zero and one")
    if float(config["inference"]["split_relative_rank_floor"]) <= 0.0:
        raise ValueError("split_relative_rank_floor must be positive")
    if float(config["inference"]["tangent_gram_min_eigenvalue_floor"]) <= 0.0:
        raise ValueError("tangent_gram_min_eigenvalue_floor must be positive")
    if float(config["inference"]["target_support_tolerance"]) <= 0.0:
        raise ValueError("target_support_tolerance must be positive")
    if config["inference"]["riesz_solver"] not in {"auto", "cg", "minres", "lsmr"}:
        raise ValueError("riesz_solver must be auto, cg, minres, or lsmr")
    if config["inference"]["variance_type"] not in {"auto", "diagonal", "spatial"}:
        raise ValueError("variance_type must be auto, diagonal, or spatial")


def canonical_json(config: dict[str, Any]) -> str:
    return json.dumps(config, sort_keys=True, separators=(",", ":"), ensure_ascii=True)


def config_hash(config: dict[str, Any]) -> str:
    return hashlib.sha256(canonical_json(config).encode("utf-8")).hexdigest()[:16]


def write_resolved_config(config: dict[str, Any], path: str | Path) -> None:
    """Write an immutable, machine-readable resolved configuration."""

    destination = Path(path)
    destination.parent.mkdir(parents=True, exist_ok=True)
    payload = {"config_hash": config_hash(config), "config": config}
    text = json.dumps(payload, indent=2, sort_keys=True) + "\n"
    if destination.exists() and destination.read_text(encoding="utf-8") != text:
        raise FileExistsError(f"resolved config already exists with different contents: {destination}")
    destination.write_text(text, encoding="utf-8")
