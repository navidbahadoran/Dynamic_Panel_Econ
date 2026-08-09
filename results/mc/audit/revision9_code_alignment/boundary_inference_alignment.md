# Boundary point-estimation and inference alignment

Locked Revision 9 distinguishes a valid finite constrained point estimate from eligibility for theorem-based inference.

The estimator is unchanged: it first tries unconstrained supplied-rank ALS; an interior solution is accepted as the constrained fast path; otherwise reconstructed entries are constrained directly to `[-10,10]` in the alternating convex subproblems. There is no entrywise clipping or surrogate norm constraint.

After alignment:

- a successfully solved finite boundary-active fit has `point_estimate_valid=true`;
- it remains retained for point reporting, bias, RMSE, and Monte Carlo SD;
- if the full coefficient fit is boundary-active, the target receives `boundary_interiority_failure`;
- if any of the four coefficient fits required by corrected inference is boundary-active, the same status is used;
- `inference_valid=false`, so the record is excluded from coverage and rejection denominators;
- numerical constrained-solver, feasibility, KKT, and nonfinite failures remain separate statuses.

The status is part of the canonical accounting vocabulary and legacy `coefficient_bound_active` records map to it. Deterministic tests cover both full-fit and split-fit boundary activity and verify point retention with inference suppression.
