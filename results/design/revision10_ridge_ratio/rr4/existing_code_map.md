# RR4 existing-code map

Audit date: 2026-08-11. This map was written before RR4 source modification.

## Reused unchanged

- `src/dynamic_panel_econ/estimation.py::fit_fixed_rank` is the maintained literal
  coefficient-box constrained joint least-squares solver. Its positive-rank blocks are
  represented by factor columns, so the represented matrix has rank at most the supplied
  number of columns and columns may collapse naturally. The interior ALS fast path and the
  alternating exact linear-box-QP fallback both solve the same constrained objective.
- `src/dynamic_panel_econ/estimation.py::_fit_fixed_rank_constrained` enforces the box inside
  every alternating subproblem; it does not clip completed coefficient matrices.
- `src/dynamic_panel_econ/rank_selection.py::fit_fixed_rank_multistart` supplies the maintained
  final fixed-rank multistart objective/stationarity checks.
- `src/dynamic_panel_econ/inference.py` and its callers construct tangent spaces, Riesz objects,
  split fits, corrections, and variances from the fit passed downstream. RR4 will continue to
  pass the final fixed-rank post-refit, never the spectral pilot.
- `src/dynamic_panel_econ/monte_carlo.py::iter_outer_task_results`, `_outer_worker`, and
  `_worker` retain deterministic outer replication parallelism, Windows spawn safety, one
  native thread per worker, and no nested pools.
- Existing supplied-rank execution through
  `src/dynamic_panel_econ/monte_carlo.py::run_replication` and
  `src/dynamic_panel_econ/rank_selection.py::fit_fixed_rank_multistart` remains separate.

## Modified for Revision 10

- `src/dynamic_panel_econ/rank_selection.py`: add explicit fitted-value RMS weights, normalized
  cap+1 spectra, the fixed ridge, zero anchor, positive-rank ratios, smallest-index argmin,
  block-vector assembly, cap+1 multistart pilot acceptance, and one final post-refit.
- `src/dynamic_panel_econ/monte_carlo.py`: route primary selected-rank mode to Revision 10,
  serialize its diagnostics, emit the frozen unresolved status, and keep a cap selection valid.
- `src/dynamic_panel_econ/config.py` and `src/dynamic_panel_econ/cli.py`: expose an explicit
  selector method. Revision 10 is primary; Revision 9 is legacy-only. Legacy IC/path fields do
  not affect Revision-10 ranks.
- `src/dynamic_panel_econ/mc_accounting.py`: recognize the frozen unresolved status.
- `README.md`: document the primary selector and legacy disposition.

## Bypassed by primary Revision 10

- `src/dynamic_panel_econ/rank_selection.py::nuclear_path`, `build_candidates`,
  `fit_rank_adaptive_cap_pilot`, `information_criterion`, `revision8_kappa`, neighboring-rank
  completion, threshold ranks, and sensitivity calculations.
- Candidate-coverage, IC-gap, threshold-sensitivity, and penalty-sensitivity fields in the
  Revision-9 record path.

## Retained for Revision-9 reproducibility

- The pre-RR4 `select_ranks` implementation will remain callable under an explicit
  `revision9_ic` method name.
- Nuclear screening, rank-adaptive cap pilot, candidate enumeration/local completion,
  Revision-8/9 IC, and historical result fields remain in place and are not deleted or
  reinterpreted.

## Current result/config locations

- Primary simulation record construction: `src/dynamic_panel_econ/monte_carlo.py::_rank_record`,
  `_failure_record`, and `run_replication`.
- Accounting: `src/dynamic_panel_econ/mc_accounting.py`.
- Defaults and TOML validation: `src/dynamic_panel_econ/config.py::DEFAULTS` and
  `validate_config`.
- CLI precedence: `src/dynamic_panel_econ/cli.py::build_run_parser` and `resolve_run_args`.
