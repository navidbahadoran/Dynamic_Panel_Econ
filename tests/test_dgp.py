import numpy as np
import pytest

from dynamic_panel_econ.dgp import (
    DGPParameters,
    coefficient_envelopes,
    generate_panel,
    generate_rank_stress_panel,
)
from dynamic_panel_econ.lowrank import numerical_rank


def test_dgp_reproducibility_shapes_ranks_and_stability():
    for dgp in range(1, 5):
        first = generate_panel(dgp, 12, 10, 1234)
        second = generate_panel(dgp, 12, 10, 1234)
        np.testing.assert_array_equal(first.y, second.y)
        assert first.y.shape == (12, 10)
        assert first.design.y_lags[0].shape == (12, 10)
        assert max(abs(first.theta0.A[0].ravel())) <= 0.85 + 1e-12
        assert [numerical_rank(matrix) for matrix in first.theta0.matrices()] == [1, 1, 1]


def test_dgp4_exact_truths_and_group_means():
    panel = generate_panel(4, 20, 16, 91)
    t0 = 16 // 2 - 1
    g1, g2 = panel.groups == 0, panel.groups == 1
    assert panel.truths["A_G1_fixed_time_true"] == pytest.approx(panel.theta0.A[0][g1, t0].mean())
    assert panel.truths["A_G2_fixed_time_true"] == pytest.approx(panel.theta0.A[0][g2, t0].mean())
    assert panel.truths["A_G2_minus_G1_fixed_time_true"] == pytest.approx(
        panel.theta0.A[0][g2, t0].mean() - panel.theta0.A[0][g1, t0].mean()
    )
    assert "A_G2_minus_G1_time_average_raw_true" in panel.truths


def test_dgp4_requires_even_n():
    with pytest.raises(ValueError):
        generate_panel(4, 9, 8, 1)


def test_predetermined_covariate_uses_lagged_not_current_primitive():
    panel = generate_panel(3, 300, 300, 2026)
    current = np.corrcoef(panel.design.x[0].ravel(), panel.u_tilde.ravel())[0, 1]
    lagged = np.corrcoef(panel.design.x[0].ravel(), panel.u_tilde_lag.ravel())[0, 1]
    assert abs(lagged) > abs(current) + 0.05
    assert abs(current) < 0.08


def test_spatial_recursion_only_from_dgp2_onward():
    iid = generate_panel(1, 120, 300, 87).u_tilde
    spatial = generate_panel(2, 120, 300, 87).u_tilde
    iid_corr = np.corrcoef(iid[:-1].ravel(), iid[1:].ravel())[0, 1]
    spatial_corr = np.corrcoef(spatial[:-1].ravel(), spatial[1:].ravel())[0, 1]
    assert abs(iid_corr) < 0.05
    assert spatial_corr == pytest.approx(0.5, abs=0.06)


@pytest.mark.parametrize("ranks", [(1, 1, 1), (2, 1, 1), (1, 0, 2)])
def test_rank_stress_generator_has_exact_ranks(ranks):
    panel = generate_rank_stress_panel(2, 20, 18, ranks, 55)
    assert tuple(numerical_rank(matrix) for matrix in panel.theta0.matrices()) == ranks
    assert np.max(np.abs(panel.theta0.A[0])) <= 0.85 + 1e-12


@pytest.mark.parametrize("dgp", [1, 2, 3, 4])
@pytest.mark.parametrize("ranks", [(1, 1, 1), (2, 1, 1), (1, 0, 2)])
def test_every_stress_design_satisfies_deterministic_common_envelope(dgp, ranks):
    c_h, c_xi = 0.7, 2.0
    theoretical = coefficient_envelopes(
        dgp, DGPParameters(), c_h=c_h, c_xi=c_xi, ranks=ranks
    )
    panel = generate_rank_stress_panel(
        dgp,
        20,
        18,
        ranks,
        551 + dgp,
        c_h=c_h,
        c_xi=c_xi,
        coefficient_bound=9.0,
        simulation_interior_margin=1.0,
    )
    assert theoretical["all"] <= 8.0
    assert panel.diagnostics["coefficient_envelope_condition_pass"] is True
    assert panel.diagnostics["realized_coefficient_envelope"] <= theoretical["all"] + 1e-12
    if ranks[0] > 1:
        assert panel.diagnostics["stress_rescale_factor_A"] < 1.0
    if ranks[2] > 1:
        assert panel.diagnostics["stress_rescale_factor_H"] < 1.0
