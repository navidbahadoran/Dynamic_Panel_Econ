# Medium pre-production diagnostic preflight: NO-GO

Date: 2026-08-07

Approved implementation commit: `75e36bc757fb399d1fff4fb748af270b605c13bb`

The 1,000-replication production run was not started. The 1,200-replication baseline medium
diagnostic and 1,080-replication rank-stress medium diagnostic were also not dispatched because
the mandatory rank-cap-pilot gate failed before a selected rank was available.

## Exact-lattice gate

- DGP: 1
- N=T: 50
- replication: 0
- pooled-R2 target: 0.65
- calibration draws: 50
- coefficient bound: 9
- simulation interior margin: 1
- baseline IC multiplier: 1
- baseline threshold multiplier: 1
- result: `rank_pilot_failure`
- elapsed time: 213.5 seconds for calibration plus the replication gate

The low-lambda nuclear-path route began at rank `(1,3,3)` and all starts were invalid, including
coefficient-envelope activity. The second route began at `(1,1,3)` but produced only one valid
start; consequently no stable pair existed under the approved objective-stability rule. The cap
pilot therefore had fewer than two valid outer routes.

Because rank selection failed, this gate produced no target rows and hence no target-support,
Riesz, tangent-Gram, or interval diagnostics. Launching the requested medium run unchanged would
not meet the stated purpose of comparing candidate coverage, IC choice, target support, and Riesz
conditioning.

## Decision required

Proceeding requires explicit author approval for a change outside the currently approved run,
such as revising the cap-pilot nuclear-path starting routes or changing the diagnostic numerical
screening specification. The coefficient bound, objective-stability rule, estimator, and
Revision-8 rank penalty were not changed during this preflight.
