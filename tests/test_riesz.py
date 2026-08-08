import numpy as np

from dynamic_panel_econ.core import Coefficients, Design, adjoint, fitted_values, inner
from dynamic_panel_econ.estimation import FactorBlock, FitResult
from dynamic_panel_econ.inference import (
    infer_corrected_target,
    infer_target,
    prepare_riesz_system,
    prepare_split_fits,
    solve_riesz,
    split_correct_target,
    tangent_gram_spectrum,
)
from dynamic_panel_econ.lowrank import tangent_project
from dynamic_panel_econ.targets import target_direction


def _fit(theta: Coefficients, ranks: tuple[int, ...]) -> FitResult:
    blocks = [FactorBlock(np.empty((theta.shape[0], 0)), np.empty((theta.shape[1], 0))) for _ in ranks]
    return FitResult(theta, ranks, 0.0, True, 1, [0.0], 0.0, 0.01, blocks)


def test_empirical_riesz_equation_for_random_tangent_direction():
    rng = np.random.default_rng(4)
    n, t = 8, 8
    design = Design([rng.normal(size=(n, t))], [rng.normal(size=(n, t))])
    matrices = [rng.normal(size=(n, 1)) @ rng.normal(size=(1, t)) for _ in range(3)]
    fitted = Coefficients([matrices[0]], [matrices[1]], matrices[2])
    d = Coefficients([rng.normal(size=(n, t))], [rng.normal(size=(n, t))], rng.normal(size=(n, t)))
    result = solve_riesz(d, fitted, design, (1, 1, 1), tolerance=1e-9, max_iter=1000)
    random_collection = Coefficients(
        [rng.normal(size=(n, t))], [rng.normal(size=(n, t))], rng.normal(size=(n, t))
    )
    delta = tangent_project(random_collection, fitted, (1, 1, 1))
    lhs = np.vdot(fitted_values(delta, design), result.weights)
    rhs = inner(d, delta)
    np.testing.assert_allclose(lhs, rhs, rtol=1e-6, atol=1e-7)
    assert result.equation_residual < 1e-6
    quotient = np.vdot(result.weights, result.weights) / sum(
        np.vdot(matrix, matrix) for matrix in result.q.matrices()
    )
    np.testing.assert_allclose(result.target_rayleigh_quotient, quotient)


def test_weighted_residual_identity_is_computed_from_fit_residual():
    rng = np.random.default_rng(8)
    n, t = 4, 4
    design = Design([rng.normal(size=(n, t))], [rng.normal(size=(n, t))])
    matrices = [rng.normal(size=(n, 1)) @ rng.normal(size=(1, t)) for _ in range(3)]
    fitted = Coefficients([matrices[0]], [matrices[1]], matrices[2])
    d = adjoint(np.ones((n, t)), design)
    residual = rng.normal(size=(n, t))
    result = solve_riesz(d, fitted, design, (1, 1, 1), residuals=residual)
    expected = abs(np.sum(result.weights * residual)) / (np.linalg.norm(result.weights) * np.linalg.norm(residual))
    assert result.weighted_residual_identity == expected


def test_tangent_gram_spectrum_uses_nonredundant_coordinates():
    rng = np.random.default_rng(81)
    n, t = 5, 4
    design = Design([rng.normal(size=(n, t))], [rng.normal(size=(n, t))])
    matrices = [rng.normal(size=(n, 1)) @ rng.normal(size=(1, t)) for _ in range(3)]
    fitted = Coefficients([matrices[0]], [matrices[1]], matrices[2])
    result = tangent_gram_spectrum(fitted, design, (1, 1, 1), tolerance=1e-8)
    assert result["tangent_gram_coordinate_dimension"] == 3 * (n + t - 1)
    assert result["tangent_gram_eigensolver_converged"] is True
    assert result["tangent_gram_largest_eigenvalue"] >= result["tangent_gram_smallest_eigenvalue"]


def test_rank_zero_entry_targets_are_unsupported_before_riesz(monkeypatch):
    rng = np.random.default_rng(91)
    n = t = 6
    design = Design([rng.normal(size=(n, t))], [rng.normal(size=(n, t))])
    rank_one = rng.normal(size=(n, 1)) @ rng.normal(size=(1, t))
    zero = np.zeros((n, t))
    y = rng.normal(size=(n, t))

    def forbidden(*args, **kwargs):
        raise AssertionError("unsupported target must not call solve_riesz")

    monkeypatch.setattr("dynamic_panel_econ.inference.solve_riesz", forbidden)
    for name, theta, ranks in (
        ("A_entry", Coefficients([zero], [rank_one], rank_one), (0, 1, 1)),
        ("B_entry", Coefficients([rank_one], [zero], rank_one), (1, 0, 1)),
    ):
        fit = _fit(theta, ranks)
        direction = target_direction(name, theta).direction
        result = infer_target(direction, fit, y, design, spatial=False)
        assert result.failure_code == "target_unsupported_selected_rank"
        assert result.riesz.iterations == 0
        assert result.riesz.target_tangent_norm == 0.0


def test_supported_target_has_positive_tangent_norm():
    rng = np.random.default_rng(92)
    n = t = 6
    design = Design([rng.normal(size=(n, t))], [rng.normal(size=(n, t))])
    matrix = rng.normal(size=(n, 1)) @ rng.normal(size=(1, t))
    theta = Coefficients([matrix], [np.zeros((n, t))], np.zeros((n, t)))
    fit = _fit(theta, (1, 0, 0))
    result = infer_target(
        target_direction("A_fixed_time_mean", theta).direction,
        fit,
        rng.normal(size=(n, t)),
        design,
        spatial=False,
    )
    assert result.riesz.target_tangent_norm > 0.0


def test_full_tangent_gram_floor_suppresses_interval():
    rng = np.random.default_rng(93)
    n = t = 5
    design = Design([rng.normal(size=(n, t))], [rng.normal(size=(n, t))])
    matrices = [rng.normal(size=(n, 1)) @ rng.normal(size=(1, t)) for _ in range(3)]
    theta = Coefficients([matrices[0]], [matrices[1]], matrices[2])
    fit = _fit(theta, (1, 1, 1))
    system = prepare_riesz_system(theta, design, fit.ranks)
    system.spectrum = {
        "tangent_gram_eigensolver_converged": True,
        "tangent_gram_smallest_eigenvalue": 1e-12,
    }
    result = infer_target(
        target_direction("A_fixed_time_mean", theta).direction,
        fit,
        rng.normal(size=(n, t)),
        design,
        spatial=False,
        tangent_gram_min_eigenvalue_floor=1e-8,
        riesz_system=system,
    )
    assert result.failure_code == "tangent_gram_nearly_singular"
    assert np.isnan(result.standard_error)


def test_split_fits_are_cached_and_split_gram_floor_suppresses(monkeypatch):
    rng = np.random.default_rng(94)
    n = t = 6
    design = Design([rng.normal(size=(n, t))], [rng.normal(size=(n, t))])
    matrix = rng.normal(size=(n, 1)) @ rng.normal(size=(1, t))
    theta = Coefficients([matrix], [np.zeros((n, t))], np.zeros((n, t)))
    fit = _fit(theta, (1, 0, 0))
    y = fitted_values(theta, design) + 0.1 * rng.normal(size=(n, t))
    groups = np.repeat([0, 1], n // 2)
    calls = 0

    def deterministic_split_fit(sub_y, sub_design, ranks, *, initial, **kwargs):
        nonlocal calls
        calls += 1
        return _fit(initial, ranks)

    monkeypatch.setattr("dynamic_panel_econ.inference.fit_fixed_rank", deterministic_split_fit)
    bundle = prepare_split_fits(
        fit,
        y,
        design,
        groups,
        time_seed=11,
        unit_seed=12,
    )
    assert calls == 4
    full_system = prepare_riesz_system(theta, design, fit.ranks)
    direction_a = target_direction("A_full_mean", theta).direction
    direction_b = target_direction("A_G1_time_average", theta, groups).direction
    first = infer_corrected_target(
        direction_a, fit, full_system, y, design, bundle, spatial=False
    )
    second = infer_corrected_target(
        direction_b, fit, full_system, y, design, bundle, spatial=False
    )
    assert calls == 4
    legacy = split_correct_target(
        direction_a,
        fit,
        y,
        design,
        groups,
        time_seed=11,
        unit_seed=12,
        spatial=False,
    )
    assert calls == 8
    assert first.estimate == legacy.estimate
    assert bundle.coefficient_fit_count == 4
    assert second.diagnostics["split_coefficient_fit_count"] == 4

    bundle.records[0].riesz_system.spectrum = {
        "tangent_gram_eigensolver_converged": True,
        "tangent_gram_smallest_eigenvalue": 0.0,
    }
    failed = infer_corrected_target(
        direction_a,
        fit,
        full_system,
        y,
        design,
        bundle,
        spatial=False,
        tangent_gram_min_eigenvalue_floor=1e-8,
    )
    assert failed.failure_code == "split_tangent_gram_nearly_singular"
    assert np.isnan(failed.standard_error)
