from pathlib import Path

import numpy as np
import pytest

from dynamic_panel_econ.calibration import (
    calibrate_cell,
    deterministic_c_h,
    load_frozen_calibrations,
    population_h_raw_variance,
    population_u_tilde_variance,
)
from dynamic_panel_econ.config import load_config, validate_config
from dynamic_panel_econ.dgp import (
    DGPParameters,
    generate_panel,
    generate_rank_stress_panel,
    rank_one_raw_envelopes,
    stress_rescale_factor,
)
from dynamic_panel_econ.monte_carlo import calibrate_design

ROOT = Path(__file__).resolve().parents[1]
FROZEN = ROOT / "configs" / "mc" / "frozen_dgp_calibration.toml"


def test_population_h_scale_is_finite_and_independent_of_calibration_shocks():
    params = DGPParameters()
    for dgp in range(1, 5):
        assert population_u_tilde_variance(dgp, params) == 1.0
        expected = deterministic_c_h(dgp, 50, params, pi_h=0.30)
        assert np.isfinite(expected)
    first = calibrate_cell(1, 20, 20, 13, params=params, draws=2)
    second = calibrate_cell(1, 20, 20, 29, params=params, draws=2)
    assert first.c_h == second.c_h
    assert first.c_h == deterministic_c_h(1, 20, params, pi_h=0.30)
    assert first.diagnostics["c_h_source"] == "analytical_population_moments"


def test_rank_two_population_h_variance_includes_exact_support_rescaling():
    params = DGPParameters()
    base_variance = population_h_raw_variance(50, params, rank=1)
    raw_envelope = rank_one_raw_envelopes(1, params)["H_raw"]
    scale = stress_rescale_factor(raw_envelope, 2, (1.0, 1.0))
    expected = scale * scale * (base_variance + 1.0)
    assert population_h_raw_variance(
        50, params, rank=2, component_strengths=(1.0, 1.0)
    ) == pytest.approx(expected, abs=1e-15)


def test_frozen_calibration_is_selected_without_using_run_seed():
    base = {
        "dgp": {"frozen_calibration_path": str(FROZEN)},
        "run": {"dgps": [1], "cells": [[50, 50]], "master_seed": 1},
    }
    first = calibrate_design(base)
    base["run"]["master_seed"] = 999999
    second = calibrate_design(base)
    key = (1, 50, 50, None)
    assert first[key] == second[key]
    assert first[key].diagnostics["calibration_source"] == "frozen_ex_ante_table"


def test_frozen_envelopes_cover_generated_truths_and_proposed_common_box():
    frozen = load_frozen_calibrations(FROZEN)
    maximum = max(float(result.diagnostics["C_Theta"]) for result in frozen.values())
    assert maximum == pytest.approx(8.288745227963506, abs=1e-12)
    assert maximum < 10.0 - 1.0
    for key, result in frozen.items():
        ranks = key[3] or (1, 1, 1)
        assert float(result.diagnostics["C_A"]) == 0.85
        expected_beta = 0.0 if ranks[1] == 0 else (
            1.71466333698683 if key[0] == 4 else 1.7774613391789282
        )
        assert float(result.diagnostics["C_beta"]) == pytest.approx(
            expected_beta, abs=1e-12
        )
        assert float(result.diagnostics["C_H"]) == pytest.approx(
            3.0 * np.sqrt(3.0) * result.c_h * result.c_xi,
            abs=1e-12,
        )
    selected = {}
    for key, calibration in frozen.items():
        dgp, _n, _t, ranks = key
        design_key = (dgp, ranks)
        selected.setdefault(design_key, calibration)
    for (dgp, ranks), calibration in selected.items():
        diagnostics = calibration.diagnostics
        if ranks is None:
            panel = generate_panel(
                dgp,
                20,
                18,
                700 + dgp,
                c_h=calibration.c_h,
                c_xi=calibration.c_xi,
            )
        else:
            panel = generate_rank_stress_panel(
                dgp,
                20,
                18,
                ranks,
                800 + 10 * dgp + sum(ranks),
                component_strengths=(1.0, 1.0),
                c_h=calibration.c_h,
                c_xi=calibration.c_xi,
            )
        realized = max(float(np.max(np.abs(matrix))) for matrix in panel.theta0.matrices())
        assert realized <= float(diagnostics["C_Theta"]) + 1e-12
        assert float(diagnostics["C_Theta"]) <= 10.0 - 1.0


def test_zero_slope_frozen_cells_keep_normalization_and_induced_r2():
    frozen = load_frozen_calibrations(FROZEN)
    zero_slope = [
        result
        for key, result in frozen.items()
        if key[3] == (1, 0, 2)
    ]
    assert zero_slope
    assert all(result.c_xi == 1.0 for result in zero_slope)
    assert all(result.target_r2 is None for result in zero_slope)
    assert all(0.0 < result.achieved_r2 < 1.0 for result in zero_slope)
    assert all(result.diagnostics["r2_scale_identified"] is False for result in zero_slope)


def test_official_production_config_activates_frozen_table_and_common_bound():
    config = load_config("configs/mc/production.toml")
    assert config["dgp"]["frozen_calibration_path"] == str(FROZEN.relative_to(ROOT)).replace(
        "\\", "/"
    )
    assert config["estimation"]["coefficient_bound"] == 10.0
    assert config["estimation"]["simulation_interior_margin"] == 1.0
    config["dgp"]["frozen_calibration_path"] = None
    with pytest.raises(ValueError, match="production run requires frozen_calibration_path"):
        validate_config(config)
