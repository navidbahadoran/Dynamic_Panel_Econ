# Constrained-solver validation status

Status: **NOT IMPLEMENTED / NOT VALIDATED**.

The prerequisite deterministic support audit failed: the maintained empirical
calibration scale `c_h` has no finite deterministic upper bound over admissible
primitive draws.  The task explicitly requires stopping at that point and
forbids inventing a calibration cap.  Accordingly:

- no constrained optimization code was added;
- no estimator/failure vocabulary was changed;
- no boundary-event replay was performed;
- no new DGP draw was generated;
- no constrained-solver tests or claimed KKT results were fabricated;
- the cancelled medium, N=400, rank-stress, power, and production runs were not
  started or resumed.

The present code's ALS stationarity residual is an unconstrained tangent-space
first-order diagnostic.  It is not a KKT residual for the entrywise-box
constrained low-rank problem and is not reported as one.

Before implementation can resume, the author must decide how the maintained
DGP/calibration specification is to provide a uniform finite H envelope without
silently changing the DGP.  The exact manuscript rule for inference on a valid
finite-sample boundary solution also needs confirmation because the manuscript
is not in this working tree.

