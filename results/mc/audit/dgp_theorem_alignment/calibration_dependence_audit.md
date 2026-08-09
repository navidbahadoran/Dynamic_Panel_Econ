# Calibration dependence audit

## Current code before this change

The old calibration drew raw panels with seed keys `"calibration"` for baseline
cells and `"rank_stress_calibration"` for stress cells.  It computed

`c_h=sqrt([pi_H/(1-pi_H)] var_sample(u_tilde)/var_sample(H_raw))`

and then solved the pooled-R2 equation for `c_xi` on those same calibration
draws.  Both constants therefore depended on realized calibration-sample
`u_tilde`, Gaussian covariate innovations, bounded coefficient primitives, and
the resulting calibration outcomes.

Production panels use a distinct semantic seed family beginning with
`"production"`.  `run_replication` receives an already computed calibration
record and passes its `c_h` and `c_xi` into the DGP generator.  Consequently:

- the old production `H_0=c_xi c_h H_raw` did **not** depend on that same
  replication's `u_it`, `epsilon_it`, or `e_it`;
- it did depend on independent realized calibration-experiment shocks through
  both `c_h` and `c_xi`;
- A and B never depended on either calibration scale;
- H depended on production bounded H primitives and the independently
  calibrated scales, while production u was `c_xi u_tilde`.

Thus there was no same-shock mechanical endogeneity to declare.  The problem
was that a realized, unbounded calibration statistic became part of a
structural coefficient scale, preventing a uniform ex-ante H envelope.

## Revised candidate conditioning structure

`c_h` is now a deterministic function of analytical population moments only.
It is independent of all calibration and production shocks.

`c_xi` remains chosen by a separate independent calibration experiment because
an exact population pooled-R2 calculation is impractical.  Its resulting value
is serialized in the candidate frozen table.  When a run explicitly supplies
`dgp.frozen_calibration_path`, `calibrate_design` reads that table and does not
solve a new R2 calibration.  Changing the run's production master seed does not
change either frozen scale.

No active medium or production configuration points to the candidate table in
this task.  The table is deliberately marked `candidate_not_activated`.

## Maintained conditional-moment interpretation

Gaussian idiosyncratic disturbances are retained.  Nothing in this calibration
change asserts that Gaussian u violates the theorem.  The change only removes
realized disturbance variance from the scale of H.  Confirmation of the
paper's exact conditioning sigma-field remains a manuscript matter, but the
candidate construction is stronger operationally: coefficient scaling is
fixed before, and independent of, every production replication.

