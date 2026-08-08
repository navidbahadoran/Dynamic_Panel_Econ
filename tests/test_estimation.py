from itertools import pairwise

import numpy as np

from dynamic_panel_econ.core import Coefficients, Design, fitted_values
from dynamic_panel_econ.estimation import fit_fixed_rank, fit_nuclear, lambda_maximum
from dynamic_panel_econ.lowrank import numerical_rank


def synthetic(seed=1, noise=0.0):
    rng = np.random.default_rng(seed)
    n, t = 10, 9
    design = Design([rng.normal(size=(n, t))], [rng.normal(size=(n, t))])
    matrices = [rng.normal(size=(n, 1)) @ rng.normal(size=(1, t)) for _ in range(3)]
    theta = Coefficients([matrices[0]], [matrices[1]], matrices[2])
    y = fitted_values(theta, design) + noise * rng.normal(size=(n, t))
    return y, design, theta


def test_joint_als_rank_zero_and_monotone_objective():
    y, design, theta = synthetic()
    initial = Coefficients([theta.A[0]], [np.zeros_like(theta.B[0])], theta.H)
    fit = fit_fixed_rank(y, design, (1, 0, 1), initial=initial, max_sweeps=100, objective_rtol=1e-9)
    np.testing.assert_array_equal(fit.theta.B[0], np.zeros_like(theta.B[0]))
    assert numerical_rank(fit.theta.A[0]) == 1
    assert all(b <= a + 1e-10 for a, b in pairwise(fit.objective_history))


def test_noiseless_joint_fit_recovers_fitted_values():
    y, design, theta = synthetic()
    fit = fit_fixed_rank(y, design, (1, 1, 1), initial=theta, max_sweeps=200, objective_rtol=1e-10)
    np.testing.assert_allclose(fitted_values(fit.theta, design), y, atol=2e-5)


def test_lambda_maximum_gives_zero_nuclear_solution():
    y, design, _ = synthetic(noise=0.1)
    fit = fit_nuclear(y, design, 1.001 * lambda_maximum(y, design), max_iter=30)
    assert max(np.linalg.norm(matrix) for matrix in fit.theta.matrices()) < 1e-8
    assert all(b <= a + 1e-10 for a, b in pairwise(fit.objective_history))
