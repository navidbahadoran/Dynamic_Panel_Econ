# Final numerical requirements and unresolved rule

## Theoretical condition

The computed joint cap+1 pilot must have normalized objective gap from the unknown global enlarged-cap
infimum equal to `o_p(zeta_NT)`. This condition is not replaced by an observable heuristic.

## Observable acceptance diagnostics

Implementation must record and apply prespecified acceptance rules for:

- feasibility of every rank and coefficient constraint;
- every start objective and stationarity/KKT residual;
- best stable objective pair and normalized objective gap across starts;
- objective-stability pass/fail and use of any additional start;
- coefficient-box activity;
- numerical ranks and singular values through every pilot cap;
- deterministic starting-value identities and ordering;
- runtime and termination status.

Numerical tolerances must be fixed before RR5 results are examined. These diagnostics make the theoretical
global-gap condition credible but do not prove it.

## Failure and point-estimate handling

If the cap+1 pilot is missing, infeasible, nonfinite, violates the acceptance tolerances, or remains unstable
across the required starts, record exactly

`rank_selection_numerically_unresolved`.

No rank vector is selected and no primary selected-rank point estimate, standard error, interval, or test is
reported for that fit. Numerical artifacts may be retained solely in diagnostic records and must not be
labeled as selected-rank estimates.

There is no silent fallback to Revision-9 IC, cap-pilot thresholding, a nuclear-path rank, supplied truth, or
another candidate. A separately prespecified supplied-rank analysis, if one exists for a different purpose,
must remain explicitly labeled and cannot replace the failed selected-rank result.

After a valid rank selection, failure of the final literal fixed-rank post-refit similarly suppresses primary
selected-rank output under a distinct post-refit numerical-failure status; the spectral pilot is never used as
the coefficient estimate.
