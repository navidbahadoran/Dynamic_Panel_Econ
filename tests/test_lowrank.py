import numpy as np

from dynamic_panel_econ.core import Coefficients
from dynamic_panel_econ.lowrank import tangent_project, tangent_project_matrix


def test_tangent_projector_idempotent_and_self_adjoint():
    rng = np.random.default_rng(18)
    fitted_matrix = rng.normal(size=(7, 1)) @ rng.normal(size=(1, 6))
    x, y = rng.normal(size=(7, 6)), rng.normal(size=(7, 6))
    px = tangent_project_matrix(x, fitted_matrix, 1)
    ppx = tangent_project_matrix(px, fitted_matrix, 1)
    py = tangent_project_matrix(y, fitted_matrix, 1)
    np.testing.assert_allclose(ppx, px, atol=1e-11)
    np.testing.assert_allclose(np.vdot(px, y), np.vdot(x, py), atol=1e-11)

    fitted = Coefficients([fitted_matrix], [fitted_matrix], fitted_matrix)
    zero = Coefficients([x], [x], x)
    projected = tangent_project(zero, fitted, (0, 1, 1))
    np.testing.assert_array_equal(projected.A[0], np.zeros_like(x))
