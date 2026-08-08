from __future__ import annotations

import numpy as np
import pytest

import dynamic_panel_econ.rank_selection as rank_module
from dynamic_panel_econ.config import load_config
from dynamic_panel_econ.core import Coefficients, Design, max_abs
from dynamic_panel_econ.estimation import FactorBlock, FitResult, NuclearFit


def _theta(ranks: tuple[int, int, int], scale: float = 1.0) -> Coefficients:
    matrices = []
    for rank in ranks:
        matrix = np.zeros((6, 6))
        for index in range(rank):
            matrix[index, index] = scale * (index + 1)
        matrices.append(matrix)
    return Coefficients([matrices[0]], [matrices[1]], matrices[2])


def _nuclear(rank: tuple[int, int, int], index: int) -> NuclearFit:
    theta = _theta(rank, 0.5)
    return NuclearFit(theta, 4.0 - index, 10.0 - index, True, 1, [10.0 - index], [[], [], []])


def _fit(ranks: tuple[int, int, int], objective: float = 1.0) -> FitResult:
    theta = _theta(ranks, 0.5)
    blocks = [FactorBlock(np.empty((6, 0)), np.empty((6, 0))) for _ in ranks]
    return FitResult(theta, ranks, objective, True, 1, [objective], 0.0, 0.1, blocks)


def test_complete_path_routes_deduplicate_rank_vectors_and_keep_occurrences() -> None:
    preliminary = [
        _nuclear((0, 0, 0), 0),
        _nuclear((1, 0, 0), 1),
        _nuclear((1, 0, 0), 2),
        _nuclear((1, 1, 0), 3),
        _nuclear((1, 1, 1), 4),
    ]
    routes, catalog = rank_module._cap_pilot_routes(
        preliminary, 0.1, (3, 3, 3), max_routes=6, coefficient_bound=9.0
    )
    ranks = [route.rank for route in routes]
    assert len(routes) >= 4
    assert len(ranks) == len(set(ranks))
    repeated = next(item for item in catalog if item["rank_vector"] == (1, 0, 0))
    assert repeated["first_path_index"] == 1
    assert repeated["best_path_index"] == 2
    assert repeated["occurrence_count"] == 2


def test_cap_start_rescaling_is_one_common_factor_without_clipping() -> None:
    theta = _theta((2, 1, 1), 12.0)
    rescaled, diagnostics = rank_module._rescale_cap_start(theta, 9.0, 0.8)
    expected = min(1.0, 0.8 * 9.0 / max_abs(theta))
    assert diagnostics["applied_common_scale"] == pytest.approx(expected)
    assert max_abs(rescaled) == pytest.approx(7.2)
    for before, after in zip(theta.matrices(), rescaled.matrices(), strict=True):
        np.testing.assert_allclose(after, expected * before)
        assert np.linalg.matrix_rank(after) == np.linalg.matrix_rank(before)


def test_strict_cap_pilot_acceptance_requires_two_valid_outer_routes(monkeypatch) -> None:
    zero = _theta((0, 0, 0))
    preliminary = [NuclearFit(zero, 1.0, 0.0, True, 1, [0.0], [[], [], []])]

    def one_valid(y, design, ranks, starts, *args, **kwargs):
        fit = _fit(ranks)
        reasons = [] if ranks == (0, 0, 0) else ["stationarity_high"]
        return fit, False, reasons, {"objective_stability_pass": not reasons}

    monkeypatch.setattr(rank_module, "_fit_candidate", one_valid)
    with pytest.raises(rank_module.RankPilotFailure) as caught:
        rank_module.fit_rank_adaptive_cap_pilot(
            np.zeros((6, 6)),
            Design([np.ones((6, 6))], [np.ones((6, 6))]),
            (2, 2, 2),
            preliminary,
            0.1,
            seed=1,
            fit_options={"coefficient_bound": 9.0},
            stationarity_tolerance=1e-6,
            start_objective_stability_tol=1e-6,
            improvement_tolerance=1e-7,
            removal_tolerance=1e-7,
            max_steps=0,
        )
    assert caught.value.diagnostics["attempted_route_count"] == 4
    assert caught.value.diagnostics["valid_route_count"] == 1
    assert caught.value.diagnostics["objective_stability_pass"] is False


def test_preflight_controls_preserve_revision8_bound_and_tau_formula() -> None:
    config = load_config("configs/mc/smoke.toml")
    assert config["estimation"]["coefficient_bound"] == 9.0
    assert config["estimation"]["simulation_interior_margin"] == 1.0
    assert config["estimation"]["rank_adaptive_max_routes"] == 6
    assert config["estimation"]["cap_pilot_start_envelope_fraction"] == 0.8
    n = t = 50
    assert np.sqrt(n * t) / np.log(n * t) == pytest.approx(50.0 / np.log(2500.0))
