# RR5e.0 pytest error map

Date: 2026-08-11

The blocked preflight collected `174 passed, 8 errors`. All eight errors occurred during pytest setup, before the affected test body ran. Pytest's configured base temporary directory is `results/pytest-parallel-audit` (`pyproject.toml`: `--basetemp=results/pytest-parallel-audit`). At session/setup preparation, `_pytest.pathlib` called `shutil.rmtree`; the terminal attempted operation was `os.rmdir` on the base directory and raised `PermissionError: [WinError 5] Access is denied`.

| Test | Filesystem path | Exception | Attempted operation | Phase | Target type | Git state | Classification |
|---|---|---|---|---|---|---|---|
| `tests/test_audit.py::test_serial_two_worker_end_to_end_equality_and_broad_schema` | `results/pytest-parallel-audit` | `PermissionError: [WinError 5]` | pytest `shutil.rmtree` / terminal `os.rmdir` | setup | directory | untracked and ignored | pytest-only basetemp |
| `tests/test_audit.py::test_rank_stress_end_to_end_runs_feasible_true_vectors` | `results/pytest-parallel-audit` | `PermissionError: [WinError 5]` | pytest `shutil.rmtree` / terminal `os.rmdir` | setup | directory | untracked and ignored | pytest-only basetemp |
| `tests/test_constrained_estimator.py::test_manifest_records_frozen_calibration_identity` | `results/pytest-parallel-audit` | `PermissionError: [WinError 5]` | pytest `shutil.rmtree` / terminal `os.rmdir` | setup | directory | untracked and ignored | pytest-only basetemp |
| `tests/test_mc_interface_accounting.py::test_reporting_refuses_silent_overwrite` | `results/pytest-parallel-audit` | `PermissionError: [WinError 5]` | pytest `shutil.rmtree` / terminal `os.rmdir` | setup | directory | untracked and ignored | pytest-only basetemp |
| `tests/test_resumability.py::test_uninterrupted_equals_interrupted_and_resumed` | `results/pytest-parallel-audit` | `PermissionError: [WinError 5]` | pytest `shutil.rmtree` / terminal `os.rmdir` | setup | directory | untracked and ignored | pytest-only basetemp |
| `tests/test_resumability.py::test_resume_refuses_fingerprint_or_task_set_mismatch` | `results/pytest-parallel-audit` | `PermissionError: [WinError 5]` | pytest `shutil.rmtree` / terminal `os.rmdir` | setup | directory | untracked and ignored | pytest-only basetemp |
| `tests/test_resumability.py::test_corrupt_and_incomplete_bundles_are_quarantined` | `results/pytest-parallel-audit` | `PermissionError: [WinError 5]` | pytest `shutil.rmtree` / terminal `os.rmdir` | setup | directory | untracked and ignored | pytest-only basetemp |
| `tests/test_resumability.py::test_manifest_exposes_pending_running_and_terminal_sets` | `results/pytest-parallel-audit` | `PermissionError: [WinError 5]` | pytest `shutil.rmtree` / terminal `os.rmdir` | setup | directory | untracked and ignored | pytest-only basetemp |

The inaccessible `.pytest_cache` produced cache-provider warnings, not any of the eight test errors. It was also untracked, ignored, and pytest-only.
