# RR4 implementation map

## Selector formulas

- Block order: `rank_selection.py::revision10_block_names`.
- RMS weights: `revision10_scale_weights`.
- Pilot singular values and normalized spectra: `revision10_normalized_spectrum`.
- Ridge: `revision10_ridge`.
- Zero and positive ratios: `revision10_ridge_ratios`.
- Smallest-index argmin: `revision10_select_block_rank`.
- Rank-vector assembly: `revision10_assemble_rank_vector`.

## Numerical estimators

- Cap+1 pilot: `rank_selection.py::fit_revision10_spectral_pilot` calls the existing literal
  joint box solver `estimation.py::fit_fixed_rank` with every requested factor-column bound equal
  to the reporting cap plus one. Pilot acceptance deliberately does not require exact numerical
  rank, so columns may collapse.
- Complete selector: `rank_selection.py::select_ranks`.
- Final estimator: `rank_selection.py::fit_fixed_rank_multistart`, called once by `select_ranks`
  at the selected vector. The returned `final_fit`, not `pilot_fit`, is passed downstream.

## Orchestration and records

- Active method routing: `config.py::DEFAULTS`, `validate_config`,
  `monte_carlo.py::_selection_options`, and `run_replication`.
- Revision-10 rank record: `monte_carlo.py::_revision10_rank_record`.
- Pilot/final fit diagnostics: `monte_carlo.py::_selected_fit_diagnostic_records`.
- Failure accounting: `mc_accounting.py::PRIMARY_STATUSES` and
  `monte_carlo.py::run_replication`.

No inference, DGP, target, split, Riesz, or variance formula was changed.
