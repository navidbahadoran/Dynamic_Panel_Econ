from __future__ import annotations

from pathlib import Path
from types import SimpleNamespace

import numpy as np
import pandas as pd
import pytest

from dynamic_panel_econ.calibration import (
    deterministic_c_h,
    load_frozen_calibrations,
)
from dynamic_panel_econ.cli import build_run_parser, resolve_run_args
from dynamic_panel_econ.config import DEFAULTS, load_config
from dynamic_panel_econ.core import Coefficients
from dynamic_panel_econ.dgp import (
    INITIAL_CONDITIONS,
    DGPParameters,
    _bounded,
    _draw_rank_stress_raw,
    _draw_raw,
    bounded_ar_envelope,
    rank_one_raw_envelopes,
    stress_rescale_factor,
)
from dynamic_panel_econ.mc_accounting import apply_retention_flags
from dynamic_panel_econ.monte_carlo import calibrate_design, classify_inference_status
from dynamic_panel_econ.targets import target_direction

ROOT = Path(__file__).resolve().parents[1]
FROZEN = ROOT / "configs" / "mc" / "frozen_dgp_calibration.toml"
MAINTAINED_CONFIGS = (
    "production.toml",
    "medium_preproduction.toml",
    "pilot.toml",
    "rank_pilot_diagnostic.toml",
    "rank_stress.toml",
    "rank_stress_medium.toml",
    "riesz_diagnostic.toml",
    "smoke.toml",
    "rank_stress_smoke.toml",
)


def test_locked_revision9_dgp_constants_and_b_support() -> None:
    params = DGPParameters()
    assert params.kappa_f_b == DEFAULTS["dgp"]["kappa_f_b"] == 0.15
    assert params.rho_fx == DEFAULTS["dgp"]["rho_fx"] == 0.5
    assert bounded_ar_envelope(params.rho_g) == pytest.approx(3.0, abs=1e-15)
    lower = params.mu_f_b - params.kappa_f_b * bounded_ar_envelope(params.rho_g)
    upper = params.mu_f_b + params.kappa_f_b * bounded_ar_envelope(params.rho_g)
    assert (lower, upper) == pytest.approx((0.15, 1.05), abs=1e-15)


def test_locked_initial_conditions_are_all_zero() -> None:
    assert INITIAL_CONDITIONS == {
        "y_i_minus_50": 0.0,
        "x_i_minus_50": 0.0,
        "g_a_minus_50": 0.0,
        "g_b_minus_50": 0.0,
        "g_h_minus_50": 0.0,
        "f_x_minus_50": 0.0,
    }


def test_locked_rank_selection_defaults_and_nuclear_path_length() -> None:
    estimation = DEFAULTS["estimation"]
    assert estimation["rank_caps"] == [3, 3, 3]
    assert estimation["coefficient_bound"] == 10.0
    assert estimation["simulation_interior_margin"] == 1.0
    assert estimation["ic_multiplier"] == 1.0
    assert estimation["threshold_multiplier"] == 1.0
    assert estimation["ic_sensitivity_multipliers"] == [0.5, 2.0]
    assert estimation["threshold_sensitivity_multipliers"] == [0.5, 2.0]
    assert estimation["nuclear_gamma"] == 0.8
    assert estimation["nuclear_epsilon"] == 0.01
    length = 1 + int(
        np.ceil(np.log(estimation["nuclear_epsilon"]) / np.log(estimation["nuclear_gamma"]))
    )
    assert length == 22
    assert estimation["dense_nuclear_gamma"] == pytest.approx(np.sqrt(0.8))


@pytest.mark.parametrize("filename", MAINTAINED_CONFIGS)
def test_maintained_configs_resolve_to_locked_revision9(filename: str) -> None:
    config = load_config(ROOT / "configs" / "mc" / filename)
    assert config["dgp"]["target_r2"] == 0.65
    assert config["dgp"]["kappa_f_b"] == 0.15
    assert config["dgp"]["rho_fx"] == 0.5
    assert config["estimation"]["rank_caps"] == [3, 3, 3]
    assert config["estimation"]["coefficient_bound"] == 10.0
    assert config["estimation"]["ic_multiplier"] == 1.0
    assert config["estimation"]["threshold_multiplier"] == 1.0
    assert config["estimation"]["ic_sensitivity_multipliers"] == [0.5, 2.0]
    assert config["estimation"]["threshold_sensitivity_multipliers"] == [0.5, 2.0]
    assert config["estimation"]["nuclear_gamma"] == 0.8
    assert config["estimation"]["nuclear_epsilon"] == 0.01


def test_canonical_cli_resolves_locked_values() -> None:
    args = build_run_parser().parse_args(
        [
            "--config",
            "configs/mc/production.toml",
            "--pooled-r2-target",
            "0.65",
            "--kappa-f-b",
            "0.15",
            "--coefficient-bound",
            "10",
            "--rank-caps",
            "3,3,3",
            "--ic-multiplier",
            "1.0",
            "--print-resolved-config",
            "--dry-run",
        ]
    )
    config = resolve_run_args(args)
    assert config["dgp"]["target_r2"] == 0.65
    assert config["dgp"]["kappa_f_b"] == 0.15
    assert config["estimation"]["coefficient_bound"] == 10.0
    assert config["estimation"]["rank_caps"] == [3, 3, 3]
    assert config["estimation"]["ic_multiplier"] == 1.0
    assert config["dgp"]["frozen_calibration_path"] == str(
        FROZEN.relative_to(ROOT)
    ).replace("\\", "/")


def test_frozen_revision9_identity_and_envelope() -> None:
    frozen = load_frozen_calibrations(FROZEN)
    assert len(frozen) == 52
    maximum = max(float(result.diagnostics["C_Theta"]) for result in frozen.values())
    assert maximum == pytest.approx(8.288745227963506, abs=1e-12)
    assert maximum <= 10.0 - 1.0
    assert 9.0 - maximum == pytest.approx(0.711254772036494, abs=1e-12)
    assert 10.0 - maximum == pytest.approx(1.711254772036494, abs=1e-12)
    for key, result in frozen.items():
        ranks = key[3] or (1, 1, 1)
        expected_beta = 0.0 if ranks[1] == 0 else (
            1.71466333698683 if key[0] == 4 else 1.7774613391789282
        )
        assert float(result.diagnostics["C_beta"]) == pytest.approx(expected_beta, abs=1e-12)
        expected_c_h = 0.7301712917987002 if ranks[2] == 2 else 0.6546536707079772
        assert result.c_h == pytest.approx(expected_c_h, abs=1e-15)
        if ranks[1] == 0:
            assert result.c_xi == 1.0
            assert result.target_r2 is None
            assert result.diagnostics["r2_scale_identified"] is False
        else:
            assert result.target_r2 == 0.65
            assert result.achieved_r2 == pytest.approx(0.65, abs=1e-10)


def test_production_uses_frozen_calibration_without_recalibration(monkeypatch) -> None:
    config = load_config(ROOT / "configs" / "mc" / "production.toml")
    config["run"]["dgps"] = [1]
    config["run"]["cells"] = [[50, 50]]
    monkeypatch.setattr(
        "dynamic_panel_econ.monte_carlo.calibrate_cell",
        lambda *args, **kwargs: pytest.fail("production attempted baseline recalibration"),
    )
    monkeypatch.setattr(
        "dynamic_panel_econ.monte_carlo.calibrate_rank_stress_cell",
        lambda *args, **kwargs: pytest.fail("production attempted stress recalibration"),
    )
    result = calibrate_design(config)
    assert result[(1, 50, 50, None)].diagnostics["calibration_source"] == (
        "frozen_ex_ante_table"
    )


def test_exact_rank_two_a_stress_formula() -> None:
    params = DGPParameters()
    seed, n, t = 731, 12, 10
    replay_rng = np.random.default_rng(seed)
    baseline = _draw_raw(4, n, t, replay_rng, params)
    added_loading = _bounded(replay_rng, n)
    added_factor = _bounded(replay_rng, t + params.burn_in)
    stress = _draw_rank_stress_raw(
        4, n, t, (2, 1, 1), np.random.default_rng(seed), params, (1.0, 1.0)
    )
    envelope = rank_one_raw_envelopes(4, params)["A_raw"]
    scale = stress_rescale_factor(envelope, 2, (1.0, 1.0))
    expected_raw = scale * (
        baseline.a_raw + added_loading[:, None] * added_factor[None, :]
    )
    expected_c_a = min(1.0, 0.85 / float(np.max(np.abs(expected_raw))))
    np.testing.assert_allclose(stress.a_raw, expected_raw, rtol=0.0, atol=0.0)
    np.testing.assert_allclose(stress.a, expected_c_a * expected_raw, rtol=0.0, atol=0.0)
    np.testing.assert_allclose(stress.beta, baseline.beta, rtol=0.0, atol=0.0)
    np.testing.assert_allclose(stress.h_raw, baseline.h_raw, rtol=0.0, atol=0.0)


def test_exact_zero_b_rank_two_h_stress_formula() -> None:
    params = DGPParameters()
    seed, n, t = 811, 12, 10
    replay_rng = np.random.default_rng(seed)
    baseline = _draw_raw(2, n, t, replay_rng, params)
    added_loading = _bounded(replay_rng, n)
    added_factor = _bounded(replay_rng, t + params.burn_in)
    stress = _draw_rank_stress_raw(
        2, n, t, (1, 0, 2), np.random.default_rng(seed), params, (1.0, 1.0)
    )
    envelope = rank_one_raw_envelopes(2, params)["H_raw"]
    scale = stress_rescale_factor(envelope, 2, (1.0, 1.0))
    expected_h = scale * (
        baseline.h_raw + added_loading[:, None] * added_factor[None, :]
    )
    np.testing.assert_allclose(stress.a, baseline.a, rtol=0.0, atol=0.0)
    np.testing.assert_array_equal(stress.beta, np.zeros_like(stress.beta))
    np.testing.assert_allclose(stress.h_raw, expected_h, rtol=0.0, atol=0.0)
    assert deterministic_c_h(
        2, t, params, pi_h=0.30, rank=2, component_strengths=(1.0, 1.0)
    ) == pytest.approx(0.7301712917987002, abs=1e-15)


def test_revision9_target_applicability() -> None:
    template = Coefficients(
        [np.ones((20, 20))], [np.ones((20, 20))], np.ones((20, 20))
    )
    groups = np.repeat([0, 1], 10)
    b_entry = target_direction("B_entry", template, groups, dgp=1)
    weak = target_direction("B_G2_minus_G1_fixed_time", template, groups, dgp=3)
    dgp4 = target_direction("B_G2_minus_G1_fixed_time", template, groups, dgp=4)
    assert b_entry.theorem_validation and b_entry.applicability == "theorem_covered"
    assert not weak.theorem_validation
    assert weak.applicability == "weak_target_stress_outside_assumption9"
    assert dgp4.theorem_validation and dgp4.applicability == "theorem_covered"


def _inference_result(*, full_boundary: bool, split_boundary: bool = False):
    split = {
        "converged": True,
        "stationarity_residual": 0.0,
        "max_envelope_ratio": 1.0,
        "rank_supported": True,
        "target_supported": True,
        "riesz_converged": True,
        "riesz_target_stable": True,
        "boundary_active": split_boundary,
    }
    return SimpleNamespace(
        failure_code=None,
        riesz=SimpleNamespace(converged=True, target_rayleigh_quotient=1.0),
        variance=1.0,
        corrected=True,
        diagnostics={"boundary_active": full_boundary, "split_fits": [dict(split) for _ in range(4)]},
    )


@pytest.mark.parametrize(
    ("full_boundary", "split_boundary"), [(True, False), (False, True)]
)
def test_boundary_activity_suppresses_inference_but_retains_point(
    full_boundary: bool, split_boundary: bool
) -> None:
    config = load_config(ROOT / "configs" / "mc" / "production.toml")
    status = classify_inference_status(
        _inference_result(full_boundary=full_boundary, split_boundary=split_boundary), config
    )
    assert status == "boundary_interiority_failure"
    row = pd.DataFrame(
        [
            {
                "dgp": 1,
                "N": 50,
                "T": 50,
                "method": "fixed_rank",
                "target": "B_entry",
                "primary_status": status,
                "estimate": 2.0,
                "standard_error": 1.0,
                "variance": 1.0,
            }
        ]
    )
    retained = apply_retention_flags(row).iloc[0]
    assert retained.point_estimate_valid
    assert retained.retained_for_bias_rmse
    assert not retained.inference_valid
    assert not retained.retained_for_coverage
    assert not retained.retained_for_rejection
