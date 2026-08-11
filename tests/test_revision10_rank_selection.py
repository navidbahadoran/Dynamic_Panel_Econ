from __future__ import annotations

import hashlib
from concurrent.futures import ProcessPoolExecutor
from multiprocessing import get_context

import numpy as np
import pytest

import dynamic_panel_econ.rank_selection as rank_module
from dynamic_panel_econ.config import DEFAULTS, validate_config
from dynamic_panel_econ.core import Coefficients, Design
from dynamic_panel_econ.estimation import FitResult
from dynamic_panel_econ.rank_selection import (
    RankPilotFailure,
    fit_revision10_spectral_pilot,
    revision10_assemble_rank_vector,
    revision10_normalized_spectrum,
    revision10_ridge,
    revision10_ridge_ratios,
    revision10_scale_weights,
    revision10_select_block_rank,
    select_ranks,
)


def _fit(
    theta: Coefficients,
    ranks: tuple[int, ...],
    objective: float = 1.0,
    *,
    converged: bool = True,
    stationarity: float = 0.0,
    envelope: float = 0.1,
) -> FitResult:
    return FitResult(
        theta=theta,
        ranks=ranks,
        objective=objective,
        converged=converged,
        iterations=2,
        objective_history=[objective + 0.1, objective],
        stationarity_residual=stationarity,
        max_envelope_ratio=envelope,
        factors=[],
        diagnostics={
            "stationarity_pass": stationarity <= 1e-6,
            "constrained_fallback_used": False,
            "boundary_active": False,
            "constrained_solver_status": "not_needed_interior_fast_path",
        },
    )


def _diagonal_theta(
    a: tuple[float, ...], b: tuple[float, ...], h: tuple[float, ...]
) -> Coefficients:
    return Coefficients([np.diag(a)], [np.diag(b)], np.diag(h))


def test_frozen_spectral_formula_retains_all_cap_plus_one_values() -> None:
    matrix = np.diag([6.0, 3.0, 1.0, 0.5])
    singular, normalized = revision10_normalized_spectrum(
        matrix,
        block_weight=2.0,
        reference_weight=4.0,
        n=4,
        t=4,
        count=3,
    )
    assert singular == pytest.approx((6.0, 3.0, 1.0))
    assert normalized == pytest.approx(
        tuple((2.0 / 4.0) ** 2 * value**2 / 16.0 for value in singular)
    )


def test_frozen_ridge_has_no_multiplier() -> None:
    assert revision10_ridge(7, 11) == pytest.approx(1.0 / np.log(77.0))


def test_rank_zero_anchor_competes_and_can_win() -> None:
    ratios = revision10_ridge_ratios((0.0, 0.0, 0.0), reporting_cap=2, n=10, t=10)
    ridge = 1.0 / np.log(100.0)
    assert ratios[0] == pytest.approx(ridge / (1.0 + ridge))
    assert revision10_select_block_rank(ratios) == 0


def test_positive_rank_unique_minimum() -> None:
    ratios = revision10_ridge_ratios((5.0, 0.0, 0.0), reporting_cap=2, n=10, t=10)
    assert revision10_select_block_rank(ratios) == 1


def test_rank_at_cap_uses_genuine_cap_plus_one_value() -> None:
    at_cap = revision10_ridge_ratios((4.0, 1.0, 0.0), reporting_cap=2, n=20, t=20)
    genuine_nonzero = revision10_ridge_ratios(
        (4.0, 1.0, 1.0), reporting_cap=2, n=20, t=20
    )
    assert revision10_select_block_rank(at_cap) == 2
    assert revision10_select_block_rank(genuine_nonzero) == 1
    assert at_cap[2] != genuine_nonzero[2]


def test_exact_tie_selects_smaller_rank_without_tolerance() -> None:
    assert revision10_select_block_rank((0.5, 0.25, 0.25, 0.8)) == 1


def test_blockwise_selection_and_vector_assembly_are_separate() -> None:
    ratios = (
        revision10_ridge_ratios((5.0, 0.0, 0.0), reporting_cap=2, n=10, t=10),
        revision10_ridge_ratios((0.0, 0.0, 0.0), reporting_cap=2, n=10, t=10),
        revision10_ridge_ratios((5.0, 2.0, 0.0), reporting_cap=2, n=10, t=10),
    )
    assert revision10_assemble_rank_vector(ratios) == (1, 0, 2)


def test_scale_equivariance_algebra_for_covariate_and_outcome_units() -> None:
    matrix = np.diag([4.0, 2.0, 0.5])
    _, baseline = revision10_normalized_spectrum(
        matrix, block_weight=3.0, reference_weight=2.0, n=3, t=3, count=3
    )
    covariate_scale = -5.0
    _, covariate_rescaled = revision10_normalized_spectrum(
        matrix / covariate_scale,
        block_weight=abs(covariate_scale) * 3.0,
        reference_weight=2.0,
        n=3,
        t=3,
        count=3,
    )
    outcome_scale = 7.0
    _, outcome_rescaled = revision10_normalized_spectrum(
        outcome_scale * matrix,
        block_weight=3.0,
        reference_weight=outcome_scale * 2.0,
        n=3,
        t=3,
        count=3,
    )
    assert covariate_rescaled == pytest.approx(baseline)
    assert outcome_rescaled == pytest.approx(baseline)
    baseline_ratios = revision10_ridge_ratios(baseline, reporting_cap=2, n=3, t=3)
    assert revision10_ridge_ratios(
        covariate_rescaled, reporting_cap=2, n=3, t=3
    ) == pytest.approx(baseline_ratios)
    assert revision10_select_block_rank(baseline_ratios) == revision10_select_block_rank(
        revision10_ridge_ratios(outcome_rescaled, reporting_cap=2, n=3, t=3)
    )


def test_scale_weights_use_uncentered_full_sample_rms_and_h_is_one() -> None:
    design = Design(
        [np.array([[1.0, 3.0], [5.0, 7.0]])],
        [np.array([[2.0, 4.0], [6.0, 8.0]])],
    )
    weights = revision10_scale_weights(design)
    assert weights == pytest.approx(
        (
            np.sqrt(np.mean(np.square(design.y_lags[0]))),
            np.sqrt(np.mean(np.square(design.x[0]))),
            1.0,
        )
    )


def test_unusable_reference_weight_is_unresolved_without_floor() -> None:
    design = Design([np.zeros((2, 2))], [np.ones((2, 2))])
    with pytest.raises(RankPilotFailure, match="unusable reference weight"):
        revision10_scale_weights(design)


def test_cap_plus_one_pilot_accepts_natural_lower_numerical_rank(monkeypatch) -> None:
    design = Design([np.ones((4, 4))], [np.ones((4, 4))])
    lower_rank = _diagonal_theta((1.0, 0.0, 0.0, 0.0), (0.0,) * 4, (1.0, 0.0, 0.0, 0.0))
    calls = 0

    def fake_fit(*args, **kwargs):
        nonlocal calls
        calls += 1
        return _fit(lower_rank, (3, 3, 3), 1.0 + calls * 1e-10)

    monkeypatch.setattr(rank_module, "fit_cap_plus_one", fake_fit)
    pilot, diagnostics = fit_revision10_spectral_pilot(
        np.ones((4, 4)),
        design,
        (2, 2, 2),
        seed=1,
        fit_options={},
        stationarity_tolerance=1e-6,
        start_objective_stability_tol=1e-6,
    )
    assert pilot.ranks == (3, 3, 3)
    assert diagnostics["numerical_rank_vector"] == (1, 0, 1)
    assert diagnostics["objective_stability_pass"] is True
    assert calls == 3


@pytest.mark.parametrize(
    ("objective", "converged", "stationarity", "envelope"),
    [
        (float("nan"), True, 0.0, 0.1),
        (1.0, True, 1.0, 0.1),
        (1.0, True, 0.0, 1.1),
    ],
)
def test_unresolved_pilot_failure_fixtures(
    monkeypatch, objective: float, converged: bool, stationarity: float, envelope: float
) -> None:
    design = Design([np.ones((4, 4))], [np.ones((4, 4))])
    theta = _diagonal_theta((1.0,) * 4, (1.0,) * 4, (1.0,) * 4)
    monkeypatch.setattr(
        rank_module,
        "fit_cap_plus_one",
        lambda *args, **kwargs: _fit(
            theta,
            (2, 2, 2),
            objective,
            converged=converged,
            stationarity=stationarity,
            envelope=envelope,
        ),
    )
    with pytest.raises(RankPilotFailure, match="rank_selection_numerically_unresolved"):
        fit_revision10_spectral_pilot(
            np.ones((4, 4)),
            design,
            (1, 1, 1),
            seed=1,
            fit_options={},
            stationarity_tolerance=1e-6,
            start_objective_stability_tol=1e-6,
        )


def test_pipeline_uses_pilot_then_one_final_post_refit(monkeypatch) -> None:
    design = Design([np.ones((4, 4))], [np.ones((4, 4))])
    pilot_theta = _diagonal_theta((8.0, 0.0, 0.0, 0.0), (0.0,) * 4, (8.0, 5.0, 0.0, 0.0))
    pilot = _fit(pilot_theta, (3, 3, 3), 0.8)
    sequence: list[str] = []

    def fake_pilot(*args, **kwargs):
        sequence.append("pilot")
        return pilot, {
            "reporting_rank_caps": (2, 2, 2),
            "pilot_rank_caps": (3, 3, 3),
            "feasibility": True,
            "finite_objective": True,
            "stationarity_residual": 0.0,
            "coefficient_box_activity": False,
            "all_start_objectives": [0.8, 0.8, 0.8],
            "all_start_stationarity_residuals": [0.0, 0.0, 0.0],
            "best_objective": 0.8,
            "second_best_valid_objective": 0.8,
            "best_two_objective_gap": 0.0,
            "objective_stability_pass": True,
        }

    final_theta = _diagonal_theta((2.0, 0.0, 0.0, 0.0), (0.0,) * 4, (2.0, 1.0, 0.0, 0.0))
    final = _fit(final_theta, (1, 0, 2), 0.7)

    def fake_final(*args, **kwargs):
        sequence.append("final_post_refit")
        assert args[2] == (1, 0, 2)
        return final, {"objective_stability_pass": True}

    monkeypatch.setattr(rank_module, "fit_revision10_spectral_pilot", fake_pilot)
    monkeypatch.setattr(rank_module, "fit_fixed_rank_multistart", fake_final)
    result = select_ranks(np.ones((4, 4)), design, (2, 2, 2), seed=4)
    assert sequence == ["pilot", "final_post_refit"]
    assert result.selected_ranks == (1, 0, 2)
    assert result.final_fit is final
    assert result.pilot_fit is pilot
    assert result.diagnostics["final_selected_rank_post_refit_status"] == "success"


def test_unresolved_pilot_never_calls_final_or_fallback(monkeypatch) -> None:
    design = Design([np.ones((4, 4))], [np.ones((4, 4))])
    final_called = False

    def fail_pilot(*args, **kwargs):
        raise RankPilotFailure("rank_selection_numerically_unresolved", {"feasibility": False})

    def forbidden_final(*args, **kwargs):
        nonlocal final_called
        final_called = True
        raise AssertionError("final/fallback selector must not be called")

    monkeypatch.setattr(rank_module, "fit_revision10_spectral_pilot", fail_pilot)
    monkeypatch.setattr(rank_module, "fit_fixed_rank_multistart", forbidden_final)
    with pytest.raises(RankPilotFailure):
        select_ranks(np.ones((4, 4)), design, (2, 2, 2), seed=4)
    assert final_called is False


def _parallel_signature(_: int) -> tuple[object, ...]:
    lag = np.arange(1, 17, dtype=float).reshape(4, 4) / 16.0
    covariate = np.flipud(lag) + 0.2
    outcome = 0.3 * lag + 0.2 * covariate + 0.1
    realization_hash = hashlib.sha256(
        np.ascontiguousarray(np.stack((outcome, lag, covariate))).view(np.uint8)
    ).hexdigest()
    result = select_ranks(
        outcome,
        Design([lag], [covariate]),
        (0, 0, 0),
        seed=7,
        max_sweeps=30,
        stationarity_tol=1.0,
        start_objective_stability_tol=1e-4,
    )
    blocks = result.diagnostics["blocks"]
    spectra = tuple(
        blocks[name]["pilot_singular_values_through_cap_plus_one"]
        for name in result.diagnostics["block_order"]
    )
    normalized = tuple(
        blocks[name]["normalized_lambda_hat_through_cap_plus_one"]
        for name in result.diagnostics["block_order"]
    )
    ratios = tuple(
        blocks[name]["ratios_R_M_0_through_cap"]
        for name in result.diagnostics["block_order"]
    )
    final_post_refit = tuple(
        tuple(matrix.ravel()) for matrix in result.final_fit.theta.matrices()
    )
    return (
        realization_hash,
        spectra,
        normalized,
        ratios,
        result.selected_ranks,
        result.final_fit.objective,
        "success",
        final_post_refit,
    )


def test_serial_parallel_revision10_semantics_are_identical() -> None:
    serial = [_parallel_signature(index) for index in range(3)]
    with ProcessPoolExecutor(max_workers=2, mp_context=get_context("spawn")) as executor:
        parallel = list(executor.map(_parallel_signature, range(3)))
    for serial_item, parallel_item in zip(serial, parallel, strict=True):
        assert parallel_item[0] == serial_item[0]
        assert parallel_item[4] == serial_item[4]
        assert parallel_item[6] == serial_item[6]
        for field in (1, 2, 3, 5, 7):
            np.testing.assert_allclose(
                np.asarray(parallel_item[field], dtype=float),
                np.asarray(serial_item[field], dtype=float),
                rtol=1e-12,
                atol=1e-12,
            )


def test_revision10_is_primary_and_legacy_ic_multiplier_is_ignored() -> None:
    from copy import deepcopy

    from dynamic_panel_econ.monte_carlo import _selection_options

    config = deepcopy(DEFAULTS)
    config["estimation"]["ic_multiplier"] = 0.0
    validate_config(config)
    active = _selection_options(config)
    assert "ic_multiplier" not in active
    assert "threshold_multiplier" not in active
    assert "nuclear_gamma" not in active
    assert "rank_adaptive_max_steps" not in active
    config["estimation"]["rank_selector_method"] = "revision9_ic"
    with pytest.raises(ValueError, match="legacy revision9_ic"):
        validate_config(config)
