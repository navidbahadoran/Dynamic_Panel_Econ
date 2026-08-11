from __future__ import annotations

from concurrent.futures import ProcessPoolExecutor
from multiprocessing import get_context

import numpy as np
import pytest

import dynamic_panel_econ.cap_plus_one as cap_module
from dynamic_panel_econ.cap_plus_one import balanced_factor_block, fit_cap_plus_one
from dynamic_panel_econ.core import Coefficients, Design, fitted_values, max_abs
from dynamic_panel_econ.estimation import fit_fixed_rank
from dynamic_panel_econ.lowrank import numerical_rank
from dynamic_panel_econ.rank_selection import fit_revision10_spectral_pilot


def _zero_design(n: int = 8, t: int = 8) -> Design:
    zero = np.zeros((n, t), dtype=np.float64)
    return Design([zero.copy()], [zero.copy()])


def _rank_product(rank: int, n: int = 8, t: int = 8) -> np.ndarray:
    if rank == 0:
        return np.zeros((n, t), dtype=np.float64)
    first = np.arange(1, n + 1, dtype=np.float64)[:, None] @ np.arange(
        1, t + 1, dtype=np.float64
    )[None, :]
    result = first / (n * t)
    if rank == 2:
        result += 0.2 * (
            np.linspace(-1.0, 1.0, n)[:, None]
            @ np.linspace(1.0, -1.0, t)[None, :]
        )
    return result


@pytest.mark.parametrize("rank", [0, 1, 2])
def test_width_four_naturally_represents_lower_rank(rank: int) -> None:
    y = _rank_product(rank)
    fit = fit_cap_plus_one(
        y,
        _zero_design(),
        (0, 0, 4),
        seed=41,
        max_sweeps=50,
    )
    assert fit.converged
    assert fit.factors[-1].loading.shape == (8, 4)
    assert fit.factors[-1].factor.shape == (8, 4)
    assert numerical_rank(fit.theta.H) == rank
    np.testing.assert_allclose(fit.theta.H, y, atol=1e-10, rtol=1e-10)
    assert fit.stationarity_residual <= 1e-6
    assert all(
        following <= previous + 1e-10
        for previous, following in zip(
            fit.objective_history, fit.objective_history[1:], strict=False
        )
    )


def test_balanced_gauge_preserves_product_objective_envelope_and_width() -> None:
    rng = np.random.default_rng(7)
    loading = rng.normal(size=(9, 4))
    factor = rng.normal(size=(7, 4))
    loading[:, 1] *= 1e8
    factor[:, 1] /= 1e8
    product = loading @ factor.T
    block = balanced_factor_block(product, 4)
    np.testing.assert_allclose(block.matrix(), product, atol=1e-12, rtol=1e-12)
    assert block.loading.shape == loading.shape
    assert block.factor.shape == factor.shape
    zero = np.zeros_like(product)
    before = Coefficients([zero.copy()], [], product)
    after = Coefficients([zero.copy()], [], block.matrix())
    design = Design([zero.copy()], [])
    y = rng.normal(size=product.shape)
    assert max_abs(before) == pytest.approx(max_abs(after), rel=1e-12, abs=1e-12)
    assert np.linalg.matrix_rank(after.H) <= 4
    np.testing.assert_allclose(
        fitted_values(before, design), fitted_values(after, design), atol=1e-12, rtol=1e-12
    )
    before_objective = np.mean(np.square(y - before.H)) / 2.0
    after_objective = np.mean(np.square(y - after.H)) / 2.0
    assert before_objective == pytest.approx(after_objective, rel=1e-12, abs=1e-12)


def test_nearly_collinear_factors_remain_finite_and_monotone() -> None:
    n = t = 10
    base = np.linspace(-1.0, 1.0, n)
    y = np.outer(base, base) + 1e-10 * np.outer(base**2, base[::-1])
    fit = fit_cap_plus_one(y, _zero_design(n, t), (0, 0, 4), seed=3, max_sweeps=50)
    assert fit.converged
    assert np.isfinite(fit.objective)
    assert all(np.isfinite(value) or np.isinf(value) for value in fit.diagnostics["condition_number_history"])
    assert fit.stationarity_residual <= 1e-6


def test_objective_safeguard_rolls_back_product_and_stored_value(monkeypatch) -> None:
    y = _rank_product(1)

    def destabilize(design, outcome, rcond):
        del outcome, rcond
        return np.full(design.shape[1], 1e6), 1.0

    monkeypatch.setattr(cap_module, "_stable_lstsq", destabilize)
    monkeypatch.setattr(
        cap_module,
        "_gauss_newton_refine",
        lambda y, design, blocks, coefficient_bound, interior_tolerance, current_objective: (
            blocks,
            current_objective,
            False,
        ),
    )
    fit = fit_cap_plus_one(y, _zero_design(), (0, 0, 4), seed=12, max_sweeps=1)
    recomputed = np.mean(np.square(y - fit.theta.H)) / 2.0
    assert fit.diagnostics["objective_safeguard_count"] == 1
    assert fit.objective_history[1] == fit.objective_history[0]
    assert fit.objective == pytest.approx(recomputed, rel=1e-14, abs=1e-14)


def test_safeguard_is_inactive_for_already_stable_solution() -> None:
    y = _rank_product(1)
    fit = fit_cap_plus_one(y, _zero_design(), (0, 0, 4), seed=12, max_sweeps=50)
    assert fit.converged
    assert fit.diagnostics["objective_safeguard_count"] == 0


def test_active_box_uses_literal_constraint_and_frozen_kkt() -> None:
    y = np.full((8, 6), 2.0)
    fit = fit_cap_plus_one(
        y,
        _zero_design(8, 6),
        (0, 0, 1),
        seed=5,
        coefficient_bound=0.5,
        max_sweeps=50,
        constrained_kkt_tolerance=1e-4,
    )
    assert fit.converged
    assert fit.diagnostics["constrained_fallback_used"] is True
    assert fit.diagnostics["constrained_algorithm"] == "deterministic_active_set_linear_box_QP"
    assert fit.diagnostics["max_constraint_violation"] <= 1e-8
    assert fit.diagnostics["constrained_KKT_residual"] <= 1e-4
    assert max_abs(fit.theta) <= 0.5 + 1e-8
    assert max_abs(fit.theta) >= 0.5 - 1e-6


def test_three_maintained_starts_pass_the_frozen_acceptance_gate() -> None:
    y = _rank_product(2, 10, 10)
    fit, diagnostics = fit_revision10_spectral_pilot(
        y,
        _zero_design(10, 10),
        (3, 3, 3),
        seed=123,
        fit_options={
            "max_sweeps": 50,
            "coefficient_bound": 10.0,
            "stationarity_tol": 1e-6,
            "constrained_kkt_tolerance": 1e-4,
        },
        stationarity_tolerance=1e-6,
        start_objective_stability_tol=1e-6,
    )
    assert fit.ranks == (4, 4, 4)
    assert diagnostics["maintained_start_count"] == 3
    assert diagnostics["valid_start_count"] == 3
    assert diagnostics["objective_stability_pass"] is True
    assert diagnostics["best_two_objective_gap"] <= 1e-6


def test_supplied_rank_legacy_solver_is_not_routed_through_cap_plus_one() -> None:
    y = _rank_product(1)
    design = _zero_design()
    initial = Coefficients(
        [np.zeros_like(y)], [np.zeros_like(y)], y.copy()
    )
    legacy = fit_fixed_rank(y, design, (0, 0, 1), initial=initial, max_sweeps=50)
    engineered = fit_cap_plus_one(y, design, (0, 0, 1), seed=9, max_sweeps=50)
    assert "solver_architecture" not in legacy.diagnostics
    np.testing.assert_allclose(legacy.theta.H, engineered.theta.H, atol=1e-10, rtol=1e-10)
    assert legacy.objective == pytest.approx(engineered.objective, abs=1e-12)


def _spawn_signature(_: int) -> tuple[float, float, tuple[float, ...]]:
    y = _rank_product(2, 7, 6)
    fit = fit_cap_plus_one(y, _zero_design(7, 6), (0, 0, 4), seed=19, max_sweeps=30)
    return fit.objective, fit.stationarity_residual, tuple(fit.theta.H.ravel())


def test_windows_spawn_determinism_for_engineered_solver() -> None:
    serial = _spawn_signature(0)
    with ProcessPoolExecutor(max_workers=2, mp_context=get_context("spawn")) as executor:
        parallel = list(executor.map(_spawn_signature, range(2)))
    for result in parallel:
        assert result[0] == pytest.approx(serial[0], abs=1e-14)
        assert result[1] == pytest.approx(serial[1], abs=1e-14)
        np.testing.assert_allclose(result[2], serial[2], atol=1e-12, rtol=1e-12)
