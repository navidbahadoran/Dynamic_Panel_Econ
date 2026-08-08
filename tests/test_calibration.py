import pytest

from dynamic_panel_econ.calibration import (
    CalibrationFeasibilityError,
    calibrate_cell,
    calibrate_rank_stress_cell,
)


def test_calibration_is_deterministic_and_hits_feasible_targets():
    # 0.65 is used here because the exact DGP's asymptotic R2 floor can exceed 0.50.
    first = calibrate_cell(1, 20, 20, 13, target_r2=0.65, draws=2)
    second = calibrate_cell(1, 20, 20, 13, target_r2=0.65, draws=2)
    assert first == second
    assert first.achieved_h_share == pytest.approx(0.30, abs=1e-12)
    assert first.achieved_r2 == pytest.approx(0.65, abs=1e-8)


def test_infeasible_r2_target_fails_loudly():
    with pytest.raises(RuntimeError, match=r"DGP 1, N=30, T=30.*target pooled R2.*lower bound"):
        calibrate_cell(1, 30, 30, 13, target_r2=0.05, draws=2)


@pytest.mark.parametrize("ranks", [(1, 1, 1), (2, 1, 1)])
def test_rank_stress_is_separately_calibrated_on_actual_matrices(ranks):
    result = calibrate_rank_stress_cell(
        1,
        50,
        50,
        ranks,
        20260807,
        component_strengths=(1.0, 1.0),
        target_r2=0.65,
        draws=3,
    )
    assert result.achieved_h_share == pytest.approx(0.30, abs=1e-12)
    assert result.achieved_r2 == pytest.approx(0.65, abs=1e-8)
    assert result.diagnostics["true_rank_vector"] == ranks


def test_zero_slope_rank_stress_reports_common_r2_infeasibility():
    with pytest.raises(CalibrationFeasibilityError, match="target pooled R2=0.650000 is infeasible"):
        calibrate_rank_stress_cell(
            1,
            50,
            50,
            (1, 0, 2),
            20260807,
            component_strengths=(1.0, 1.0),
            target_r2=0.65,
            draws=3,
        )
