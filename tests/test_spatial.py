import numpy as np

from dynamic_panel_econ.inference import corrected_scores
from dynamic_panel_econ.spatial import (
    bartlett_matrix,
    bartlett_quadratic,
    bartlett_quadratic_dense,
)


def test_cutoff_zero_is_diagonal():
    rng = np.random.default_rng(5)
    scores = rng.normal(size=(9, 7))
    assert bartlett_quadratic(scores, 0) == np.sum(scores**2)


def test_dense_and_lag_sum_agree_and_matrix_is_psd():
    rng = np.random.default_rng(5)
    scores = rng.normal(size=(9, 7))
    np.testing.assert_allclose(bartlett_quadratic(scores, 4), bartlett_quadratic_dense(scores, 4))
    eigenvalues = np.linalg.eigvalsh(bartlett_matrix(20, 5))
    assert eigenvalues.min() >= -1e-12


def test_corrected_score_uses_distinct_split_residuals():
    shape = (3, 2)
    full_w = np.full(shape, 2.0)
    time_w = np.full(shape, 3.0)
    unit_w = np.full(shape, 5.0)
    full_r = np.full(shape, 7.0)
    time_r = np.full(shape, 11.0)
    unit_r = np.full(shape, 13.0)
    actual = corrected_scores(full_w, full_r, time_w, time_r, unit_w, unit_r)
    np.testing.assert_array_equal(actual, 3 * full_w * full_r - time_w * time_r - unit_w * unit_r)
