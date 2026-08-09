import hashlib
import inspect
import json
from pathlib import Path

import numpy as np

import dynamic_panel_econ.estimation as estimation_module
import dynamic_panel_econ.inference as inference_module
import dynamic_panel_econ.rank_selection as rank_module
from dynamic_panel_econ.config import load_config
from dynamic_panel_econ.core import Coefficients, Design, max_abs
from dynamic_panel_econ.estimation import fit_fixed_rank
from dynamic_panel_econ.inference import infer_target, prepare_split_fits
from dynamic_panel_econ.lowrank import numerical_rank
from dynamic_panel_econ.monte_carlo import run_monte_carlo


def _intercept_fixture(level: float, n: int = 8, t: int = 6):
    zero = np.zeros((n, t), dtype=float)
    design = Design([zero.copy()], [zero.copy()])
    y = np.full((n, t), level, dtype=float)
    initial = Coefficients([zero.copy()], [zero.copy()], y.copy())
    return y, design, initial


def test_interior_fast_path_is_the_literal_inactive_constraint_solution():
    y, design, initial = _intercept_fixture(0.2)
    fit = fit_fixed_rank(y, design, (0, 0, 1), initial=initial, coefficient_bound=1.0)
    assert fit.converged
    assert fit.diagnostics["unconstrained_inside_box"] is True
    assert fit.diagnostics["constrained_fallback_used"] is False
    assert fit.objective < 1e-20
    assert max_abs(fit.theta) < 1.0


def test_outside_solution_uses_exact_box_fallback_and_retains_boundary_fit():
    y, design, initial = _intercept_fixture(2.0)
    fit = fit_fixed_rank(
        y,
        design,
        (0, 0, 1),
        initial=initial,
        coefficient_bound=0.5,
        max_sweeps=50,
        constrained_kkt_tolerance=1e-6,
    )
    assert fit.converged
    assert fit.diagnostics["unconstrained_outside_box"] is True
    assert fit.diagnostics["constrained_fallback_used"] is True
    assert fit.diagnostics["boundary_active"] is True
    assert fit.diagnostics["constrained_solver_status"] == "success"
    assert fit.diagnostics["max_constraint_violation"] <= 1e-8
    assert fit.diagnostics["constrained_KKT_residual"] <= 1e-6
    assert max_abs(fit.theta) <= 0.5 + 1e-8
    assert numerical_rank(fit.theta.H) == 1
    assert fit.objective >= fit.diagnostics["unconstrained_objective"]


def test_constrained_algorithm_contains_no_entrywise_clipping():
    source = inspect.getsource(estimation_module._fit_fixed_rank_constrained)
    source += inspect.getsource(estimation_module._solve_linear_box_subproblem)
    assert "clip" not in source
    assert "LinearConstraint" in source


def test_constrained_solver_failure_has_explicit_status(monkeypatch):
    y, design, initial = _intercept_fixture(2.0)

    def fail(design, outcome, constraints, initial, bound, **kwargs):
        return initial, False, "forced deterministic failure", 1

    monkeypatch.setattr(estimation_module, "_solve_linear_box_subproblem", fail)
    fit = fit_fixed_rank(
        y,
        design,
        (0, 0, 1),
        initial=initial,
        coefficient_bound=0.5,
    )
    assert not fit.converged
    assert fit.diagnostics["constrained_solver_status"] == "constrained_solver_failure"
    assert fit.diagnostics["constrained_fallback_used"] is True


def test_all_fixed_rank_paths_share_theorem_aligned_solver():
    assert rank_module.fit_fixed_rank is estimation_module.fit_fixed_rank
    assert inference_module.fit_fixed_rank is estimation_module.fit_fixed_rank

    y, design, initial = _intercept_fixture(2.0)
    full = fit_fixed_rank(
        y,
        design,
        (0, 0, 1),
        initial=initial,
        coefficient_bound=0.5,
        max_sweeps=50,
        constrained_kkt_tolerance=1e-6,
    )
    bundle = prepare_split_fits(
        full,
        y,
        design,
        np.repeat([0, 1], 4),
        time_seed=11,
        unit_seed=12,
        fit_options={
            "coefficient_bound": 0.5,
            "max_sweeps": 50,
            "constrained_kkt_tolerance": 1e-6,
        },
    )
    assert len(bundle.records) == 4
    assert all(record.fit.diagnostics["constrained_fallback_used"] for record in bundle.records)
    assert all(record.fit.diagnostics["constrained_solver_status"] == "success" for record in bundle.records)


def test_boundary_fit_uses_stated_inference_computation_without_zero_identity():
    y, design, initial = _intercept_fixture(2.0)
    fit = fit_fixed_rank(
        y,
        design,
        (0, 0, 1),
        initial=initial,
        coefficient_bound=0.5,
        max_sweeps=50,
        constrained_kkt_tolerance=1e-6,
    )
    direction = Coefficients(
        [np.zeros_like(y)],
        [np.zeros_like(y)],
        np.full_like(y, 1.0 / y.size),
    )
    result = infer_target(direction, fit, y, design, spatial=False)
    assert np.isfinite(result.estimate)
    assert result.diagnostics["boundary_active"] is True
    assert result.diagnostics["normal_equation_used_as_identity"] is False
    assert result.riesz.weighted_residual_identity is not None


def test_manifest_records_frozen_calibration_identity(monkeypatch, tmp_path):
    config = load_config("configs/mc/smoke.toml")
    config["run"].update(
        {
            "name": "frozen_manifest_test",
            "dgps": [1],
            "cells": [[50, 50]],
            "replications": 1,
            "chunk_size": 1,
            "parallel_level": "none",
            "n_jobs": 1,
            "output_root": str(tmp_path),
            "rank_mode": "fixed",
        }
    )
    config["inference"]["targets"] = []
    monkeypatch.setattr("dynamic_panel_econ.monte_carlo.run_replication", lambda *args: [])
    _, root = run_monte_carlo(config, overwrite=True)
    manifest = json.loads((root / "run_manifest.json").read_text(encoding="utf-8"))
    frozen = Path(config["dgp"]["frozen_calibration_path"])
    assert Path(manifest["frozen_calibration_file"]) == frozen
    assert manifest["frozen_calibration_hash"] == hashlib.sha256(frozen.read_bytes()).hexdigest()
    assert manifest["B"] == 10.0
    assert manifest["c_B"] == 1.0
    assert len(manifest["calibration_cells"]) == 1
    cell = manifest["calibration_cells"][0]
    assert cell["calibration_cell"] == {
        "dgp": 1,
        "N": 50,
        "T": 50,
        "true_rank_vector": None,
    }
    assert cell["c_h"] == 0.6546536707079772
    assert cell["C_Theta"] <= cell["B"] - cell["c_B"]
