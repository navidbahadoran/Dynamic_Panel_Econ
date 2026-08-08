from __future__ import annotations

import numpy as np
import pytest

import dynamic_panel_econ.rank_selection as rank_module
from dynamic_panel_econ.core import Coefficients, Design, max_abs
from dynamic_panel_econ.estimation import FactorBlock, FitResult

RANKS = (1, 1, 1)


def _theta(scale: float = 0.5) -> Coefficients:
    matrices = []
    for block in range(3):
        left = np.arange(1, 7, dtype=float) + block
        right = np.arange(6, 0, -1, dtype=float) + block
        matrices.append(scale * np.outer(left, right) / np.max(np.outer(left, right)))
    return Coefficients([matrices[0]], [matrices[1]], matrices[2])


def _fit(
    objective: float,
    *,
    envelope_ratio: float = 0.2,
    theta: Coefficients | None = None,
) -> FitResult:
    value = _theta() if theta is None else theta
    blocks = [FactorBlock(np.empty((6, 0)), np.empty((6, 0))) for _ in RANKS]
    return FitResult(
        value,
        RANKS,
        objective,
        True,
        2,
        [objective + 0.1, objective],
        1e-8,
        envelope_ratio,
        blocks,
        {"runtime_seconds": 0.01},
    )


def _design() -> Design:
    ones = np.ones((6, 6))
    return Design([ones], [ones])


def test_confirmation_perturbations_preserve_rank_and_interior_without_clipping() -> None:
    base = _theta()
    for magnitude in rank_module.BEST_BASIN_PERTURBATION_MAGNITUDES:
        perturbed = rank_module._rank_preserving_perturbation(base, RANKS, magnitude)
        assert tuple(np.linalg.matrix_rank(matrix) for matrix in perturbed.matrices()) == RANKS
        assert max_abs(perturbed) < 9.0
        assert any(
            not np.allclose(before, after)
            for before, after in zip(base.matrices(), perturbed.matrices(), strict=True)
        )
        assert not any(np.any(np.abs(matrix) == 9.0) for matrix in perturbed.matrices())


@pytest.mark.parametrize(
    ("objectives", "expected"),
    [([1.0 + 2e-6, 1.0 + 4e-6, 1.1], True), ([1.1, 1.2, 1.3], False)],
)
def test_confirmation_requires_two_valid_matches(monkeypatch, objectives, expected) -> None:
    returned = iter(_fit(value) for value in objectives)
    monkeypatch.setattr(rank_module, "fit_fixed_rank", lambda *args, **kwargs: next(returned))
    _, records, passed = rank_module._confirm_best_basin(
        np.zeros((6, 6)),
        _design(),
        _fit(1.0),
        RANKS,
        fit_options={"coefficient_bound": 9.0},
        seed=123,
        stationarity_tolerance=1e-4,
        objective_tolerance=1e-5,
        start_envelope_fraction=0.8,
        diagnostic_context="test",
    )
    assert passed is expected
    assert sum(record["confirmation_pass"] for record in records) == (2 if expected else 0)
    assert all(record["route_type"] == "basin_confirmation" for record in records)


def test_fixed_rank_uses_three_deterministic_starts_and_requires_stability(monkeypatch) -> None:
    objectives = iter([1.0, 1.0 + 5e-6, 1.2])
    calls = []

    def fitted(*args, **kwargs):
        calls.append(kwargs)
        return _fit(next(objectives))

    monkeypatch.setattr(rank_module, "fit_fixed_rank", fitted)
    _, diagnostics = rank_module.fit_fixed_rank_multistart(
        np.zeros((6, 6)),
        _design(),
        RANKS,
        seed=99,
        fit_options={"coefficient_bound": 9.0},
        stationarity_tolerance=1e-4,
        start_objective_stability_tol=1e-5,
    )
    assert len(calls) == 3
    assert len({call["seed"] for call in calls}) == 3
    assert diagnostics["objective_stability_tolerance"] == 1e-5
    assert diagnostics["objective_stability_pass"]
    assert diagnostics["final_acceptance_basis"] == "original_route_stability"


def test_fixed_rank_bound_active_best_solution_remains_invalid(monkeypatch) -> None:
    objectives = iter([1.0, 1.0, 1.1])
    monkeypatch.setattr(
        rank_module,
        "fit_fixed_rank",
        lambda *args, **kwargs: _fit(next(objectives), envelope_ratio=1.01),
    )
    chosen, diagnostics = rank_module.fit_fixed_rank_multistart(
        np.zeros((6, 6)),
        _design(),
        RANKS,
        seed=99,
        fit_options={"coefficient_bound": 9.0},
        stationarity_tolerance=1e-4,
        start_objective_stability_tol=1e-5,
    )
    assert chosen.max_envelope_ratio >= 1.0
    assert not diagnostics["objective_stability_pass"]
    assert diagnostics["final_acceptance_basis"] == "failure"

