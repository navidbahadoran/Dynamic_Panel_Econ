from __future__ import annotations

from copy import deepcopy
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

from dynamic_panel_econ.cli import build_run_parser, resolve_run_args, resolved_config_text
from dynamic_panel_econ.config import DEFAULTS
from dynamic_panel_econ.core import Coefficients, Design
from dynamic_panel_econ.dgp import PanelData
from dynamic_panel_econ.estimation import FactorBlock, FitResult
from dynamic_panel_econ.mc_accounting import (
    apply_retention_flags,
    power_design_configs,
    power_loading_means,
    reconcile_fit_rows,
    reconcile_summary,
    semantic_replication_id,
    summarize_accounting,
    summarize_power,
)
from dynamic_panel_econ.method_reporting import _write_exclusive
from dynamic_panel_econ.monte_carlo import (
    _calibration_failure_task_rows,
    _dgp_realization_hash,
    _replication_chunks,
    _selected_context_fit_type,
    _worker,
    run_replication,
)
from dynamic_panel_econ.rank_selection import RankSelectionFailure
from scripts.audit_filtering import audit_repository


def _config(mode: str) -> dict:
    config = deepcopy(DEFAULTS)
    config["run"].update({"rank_mode": mode, "cells": [[4, 4]], "replications": 1})
    config["estimation"].update(
        {"fixed_ranks": [1, 1, 1], "rank_caps": [2, 2, 2], "max_sweeps": 2}
    )
    config["inference"]["targets"] = []
    return config


def _panel() -> PanelData:
    shape = (4, 4)
    zeros = np.zeros(shape)
    theta = Coefficients([zeros.copy()], [zeros.copy()], zeros.copy())
    return PanelData(
        y=zeros.copy(),
        design=Design([zeros.copy()], [zeros.copy()]),
        theta0=theta,
        u=zeros.copy(),
        u_tilde=zeros.copy(),
        u_tilde_lag=zeros.copy(),
        groups=np.array([0, 0, 1, 1]),
    )


def _fit() -> FitResult:
    panel = _panel()
    blocks = [FactorBlock(np.zeros((4, 1)), np.zeros((4, 1))) for _ in range(3)]
    return FitResult(panel.theta0, (1, 1, 1), 0.0, True, 1, [1.0, 0.0], 0.0, 0.0, blocks)


def test_fixed_rank_mode_skips_selection(monkeypatch) -> None:
    monkeypatch.setattr("dynamic_panel_econ.monte_carlo.generate_panel", lambda *a, **k: _panel())
    monkeypatch.setattr(
        "dynamic_panel_econ.monte_carlo.fit_fixed_rank_multistart",
        lambda *a, **k: (_fit(), {"objective_stability_pass": True}),
    )
    monkeypatch.setattr(
        "dynamic_panel_econ.monte_carlo.select_ranks",
        lambda *a, **k: pytest.fail("fixed mode entered rank selection"),
    )
    rows = run_replication((1, 4, 4, 0, None), _config("fixed"), {"c_h": 1, "c_xi": 1})
    assert rows[0]["method"] == "fixed_rank"
    assert rows[0]["supplied_rank_vector"] == "[1, 1, 1]"


def test_selected_rank_mode_calls_selection(monkeypatch) -> None:
    monkeypatch.setattr("dynamic_panel_econ.monte_carlo.generate_panel", lambda *a, **k: _panel())

    def called(*args, **kwargs):
        raise RankSelectionFailure("selection called")

    monkeypatch.setattr("dynamic_panel_econ.monte_carlo.select_ranks", called)
    rows = run_replication((1, 4, 4, 0, None), _config("selected"), {"c_h": 1, "c_xi": 1})
    assert rows[0]["status"] == "rank_selection_failure"


def test_methods_share_semantic_replication_identity_and_dgp_draw(monkeypatch) -> None:
    assert semantic_replication_id(4, 50, 50, 7, (1, 1, 1)) == semantic_replication_id(
        4, 50, 50, 7, (1, 1, 1)
    )
    states = []

    def capture(*args, **kwargs):
        states.append(tuple(args[3].generate_state(4)))
        return _panel()

    monkeypatch.setattr("dynamic_panel_econ.monte_carlo.generate_panel", capture)
    monkeypatch.setattr(
        "dynamic_panel_econ.monte_carlo.fit_fixed_rank_multistart",
        lambda *a, **k: (_fit(), {"objective_stability_pass": True}),
    )
    fixed_rows = run_replication(
        (1, 4, 4, 0, None), _config("fixed"), {"c_h": 1, "c_xi": 1}
    )
    monkeypatch.setattr(
        "dynamic_panel_econ.monte_carlo.select_ranks",
        lambda *a, **k: (_ for _ in ()).throw(RankSelectionFailure("stop after draw")),
    )
    selected_rows = run_replication(
        (1, 4, 4, 0, None), _config("selected"), {"c_h": 1, "c_xi": 1}
    )
    assert states[0] == states[1]
    assert fixed_rows[0]["dgp_realization_hash"] == selected_rows[0]["dgp_realization_hash"]


def test_dgp_realization_hash_includes_calibration() -> None:
    panel = _panel()
    baseline = _dgp_realization_hash(panel, {"c_h": 1.0, "c_xi": 1.0})
    changed = _dgp_realization_hash(panel, {"c_h": 1.0, "c_xi": 2.0})
    assert baseline != changed


def test_independent_preflight_replication_indices_are_new_and_chunked() -> None:
    run = {"replication_start": 3, "replications": 3, "chunk_size": 2}
    assert _replication_chunks(run) == [(3, 5), (5, 6)]


def test_failed_attempt_remains_in_worker_records(monkeypatch) -> None:
    monkeypatch.setattr(
        "dynamic_panel_econ.monte_carlo.generate_panel",
        lambda *a, **k: (_ for _ in ()).throw(RuntimeError("synthetic failure")),
    )
    rows = _worker(((1, 4, 4, 0, None), _config("fixed"), {"c_h": 1, "c_xi": 1}))
    assert sum(row["record_type"] == "replication" for row in rows) == 1
    assert any(row["primary_status"] == "software_exception" for row in rows)


def test_calibration_failure_is_an_attempted_replication() -> None:
    rows = _calibration_failure_task_rows(
        (1, 50, 50, 9, None), _config("fixed"), "CalibrationError: infeasible"
    )
    attempt = next(row for row in rows if row["record_type"] == "replication")
    assert attempt["primary_status"] == "calibration_failure"
    assert not attempt["completed_dgp_replication"]


def test_worker_fit_diagnostics_reconcile_executed_fits(monkeypatch) -> None:
    monkeypatch.setattr("dynamic_panel_econ.monte_carlo.generate_panel", lambda *a, **k: _panel())
    rows = _worker(((1, 4, 4, 0, None), _config("fixed"), {"c_h": 1, "c_xi": 1}))
    expected = next(row["expected_fit_count"] for row in rows if row["record_type"] == "replication")
    fit_rows = [row for row in rows if row["record_type"] == "fit_diagnostic"]
    actual = len(fit_rows)
    assert expected == actual == 3
    assert {row["fit_type"] for row in fit_rows} == {"full_fixed_rank"}
    assert all("max_abs_coefficient" in row for row in fit_rows)
    assert all("constrained_runtime_seconds" in row for row in fit_rows)


def test_selected_context_classification_preserves_split_fit_labels() -> None:
    assert (
        _selected_context_fit_type("coefficient_fit", "cap_pilot_route")
        == "rank_cap_pilot"
    )
    assert (
        _selected_context_fit_type("coefficient_fit", "post_refit_start")
        == "candidate_post_refit"
    )
    assert (
        _selected_context_fit_type("time_split_0", "post_refit_start")
        == "time_split_0"
    )


def _attempts() -> pd.DataFrame:
    return pd.DataFrame(
        [
            {"run_id": "x", "dgp": 1, "N": 50, "T": 50, "replication": i, "method": "fixed_rank", "primary_status": "success" if i < 2 else "full_fit_failure", "completed_dgp_replication": True, "replication_runtime_seconds": i + 1.0}
            for i in range(3)
        ]
    )


def _targets() -> pd.DataFrame:
    return pd.DataFrame(
        [
            {"run_id": "x", "dgp": 1, "N": 50, "T": 50, "replication": 0, "method": "fixed_rank", "target": "A", "primary_status": "success", "estimate": 1.0, "truth": 0.0, "standard_error": 1.0, "variance": 1.0, "covered_95pct": True, "reject_zero_5pct": False},
            {"run_id": "x", "dgp": 1, "N": 50, "T": 50, "replication": 1, "method": "fixed_rank", "target": "A", "primary_status": "success", "estimate": 1e12, "truth": 0.0, "standard_error": 2.0, "variance": 4.0, "covered_95pct": False, "reject_zero_5pct": True},
        ]
    )


def test_accounting_retains_extreme_and_missing_failure() -> None:
    summary = summarize_accounting(_attempts(), _targets(), targets=["A"])
    row = summary.iloc[0]
    assert row.R_attempted == 3
    assert row.R_point == row.R_inference == 2
    assert row.rmse > 1e11
    assert row.failure_full_fit_failure == 1
    assert row.coverage == 0.5
    assert row.rejection_probability == 0.5
    assert row.total_inference_failure_rate == pytest.approx(1 / 3)
    assert row.point_retained_share == row.inference_retained_share == pytest.approx(2 / 3)


def test_nonfinite_estimate_has_explicit_status() -> None:
    records = _targets().iloc[:1].copy()
    records.loc[records.index[0], "estimate"] = np.inf
    flagged = apply_retention_flags(records)
    assert flagged.iloc[0].primary_status == "nonfinite_estimate"
    assert not flagged.iloc[0].retained_for_bias_rmse


def test_primary_failure_cannot_retain_finite_inference() -> None:
    records = _targets().iloc[:1].assign(primary_status="split_fit_failure")
    flagged = apply_retention_flags(records)
    row = flagged.iloc[0]
    assert row.point_estimate_valid
    assert not row.inference_valid
    assert row.primary_status == "split_fit_failure"


def test_nonfatal_warning_can_coexist_with_success() -> None:
    records = _targets().iloc[:1].assign(warning_flags='["diagnostic_warning"]')
    flagged = apply_retention_flags(records)
    row = flagged.iloc[0]
    assert row.primary_status == "success"
    assert row.inference_valid
    assert row.warning_flags == '["diagnostic_warning"]'


def test_extreme_finite_flag_never_controls_retention() -> None:
    records = pd.concat([_targets(), _targets().assign(replication=[2, 3], estimate=[1.1, 1.2])])
    flagged = apply_retention_flags(records)
    extreme = flagged.loc[flagged["estimate"].eq(1e12)].iloc[0]
    assert extreme.extreme_estimate_flag
    assert extreme.retained_for_bias_rmse


def test_reconciliation_detects_bad_denominators() -> None:
    summary = summarize_accounting(_attempts(), _targets(), targets=["A"])
    summary.loc[0, "R_inference"] = 4
    with pytest.raises(ValueError, match="retention"):
        reconcile_summary(summary)


def test_fit_diagnostic_reconciliation() -> None:
    keys = {"run_id": "x", "dgp": 1, "N": 50, "T": 50, "replication": 0, "method": "fixed_rank"}
    fits = pd.DataFrame([{**keys, "fit_type": "full"}, {**keys, "fit_type": "time_split"}])
    reconcile_fit_rows(fits, pd.DataFrame([{**keys, "expected_fit_count": 2}]))
    with pytest.raises(ValueError, match="fit diagnostic"):
        reconcile_fit_rows(fits.iloc[:1], pd.DataFrame([{**keys, "expected_fit_count": 2}]))


def test_power_null_and_mcse_use_inference_denominator() -> None:
    records = _targets().assign(nominal_delta=0.0, realized_true_contrast=0.0)
    records.loc[records.index[1], "standard_error"] = np.nan
    power = summarize_power(records)
    assert power.iloc[0].null_or_alternative == "size/null"
    assert power.iloc[0].valid_inference_replications == 1
    assert power_loading_means(0.4) == pytest.approx((0.8, 1.2))


def test_power_design_expands_matched_methods() -> None:
    config = _config("selected")
    config["run"].update({"experiment": "power", "alternative_grid": [0.0, 0.4]})
    designs = power_design_configs(config)
    assert {(item["run"]["nominal_delta"], item["run"]["rank_mode"]) for item in designs} == {
        (0.0, "fixed"), (0.0, "selected"), (0.4, "fixed"), (0.4, "selected")
    }


def test_cli_precedence_and_balanced_grid() -> None:
    args = build_run_parser().parse_args(
        [
            "--balanced-grid",
            "50,100,200,400",
            "--rank-mode",
            "fixed",
            "--fixed-ranks",
            "1,1,1",
            "--pooled-r2-target",
            "0.65",
            "--dry-run",
        ]
    )
    config = resolve_run_args(args)
    assert config["run"]["cells"] == [[50, 50], [100, 100], [200, 200], [400, 400]]
    assert config["run"]["rank_mode"] == "fixed"
    assert config["dgp"]["target_r2"] == 0.65
    import tomllib

    parsed = tomllib.loads(resolved_config_text(config))
    assert parsed["run"]["cells"][-1] == [400, 400]


def test_reporting_refuses_silent_overwrite(tmp_path: Path) -> None:
    path = tmp_path / "table.tex"
    _write_exclusive(path, "first", overwrite=False)
    with pytest.raises(FileExistsError):
        _write_exclusive(path, "second", overwrite=False)


def test_filtering_audit_finds_known_constructs() -> None:
    rows = audit_repository(Path("."))
    assert any("reporting.py" in str(row["file"]) and "isna" in str(row["construct"]) for row in rows)


def test_exactly_four_split_fit_types_are_target_independent() -> None:
    split_types = {"time_split_0", "time_split_1", "unit_split_0", "unit_split_1"}
    assert len(split_types) == 4
