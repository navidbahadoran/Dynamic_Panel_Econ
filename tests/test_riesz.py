import numpy as np

from dynamic_panel_econ.core import Coefficients, Design, adjoint, fitted_values, inner
from dynamic_panel_econ.inference import solve_riesz, tangent_gram_spectrum
from dynamic_panel_econ.lowrank import tangent_project


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
