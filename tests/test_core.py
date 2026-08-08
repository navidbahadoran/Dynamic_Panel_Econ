import numpy as np

from dynamic_panel_econ.core import Coefficients, Design, adjoint, fitted_values, inner


def test_adjoint_identity():
    rng = np.random.default_rng(4)
    design = Design([rng.normal(size=(5, 4))], [rng.normal(size=(5, 4))])
    delta = Coefficients(
        [rng.normal(size=(5, 4))], [rng.normal(size=(5, 4))], rng.normal(size=(5, 4))
    )
    residual = rng.normal(size=(5, 4))
    np.testing.assert_allclose(
        np.vdot(fitted_values(delta, design), residual),
        inner(delta, adjoint(residual, design)),
        rtol=1e-13,
        atol=1e-13,
    )
