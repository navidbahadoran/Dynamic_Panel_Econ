# RR4 numerical failure rules

## Cap+1 pilot

Each maintained start must return a finite objective, a box-feasible coefficient collection, and
the maintained stationarity/KKT acceptance. Lower numerical rank than cap+1 is permitted. Valid
start objectives are sorted and the best two must have normalized gap no greater than the existing
`start_objective_stability_tol`.

If fewer than two valid stable starts remain, or a later weight/spectrum/ratio input is nonfinite,
the exact status is `rank_selection_numerically_unresolved`. The attempted diagnostics are retained.
No rank, Revision-9 IC result, truth, threshold rank, path rank, primary point estimate, or primary
inference is substituted.

These checks are numerical evidence only and do not certify the theoretical global objective-gap
condition.

## Final post-refit

The selected-rank fit must preserve its literal supplied rank, feasibility, maintained
stationarity/KKT rule, and multistart objective stability. Failure has the distinct status
`selected_rank_post_refit_numerically_unresolved`; the pilot is never reported in its place.

Boundary-active point/inference handling after a successful fit remains unchanged.
