# Legacy Revision-9 disposition

- The former `select_ranks` implementation is preserved as
  `rank_selection.py::select_ranks_revision9`.
- It is reachable only when `estimation.rank_selector_method = "revision9_ic"`.
- Nuclear screening, the rank-adaptive cap pilot, threshold ranks, candidate enumeration,
  neighboring completion, IC calculation, and sensitivity calculations are unchanged inside that
  legacy function.
- Primary selected mode defaults to `revision10_ridge_ratio`.
- `monte_carlo.py::_selection_options` does not pass any nuclear, threshold, IC, candidate, or
  sensitivity field to the Revision-10 function. Legacy fields therefore cannot affect its ranks.
- Existing Revision-9 tests explicitly invoke the legacy function and cap-hit behavior.

No historical output or source implementation was deleted.
