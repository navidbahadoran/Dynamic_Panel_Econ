# Final small independent preflight audit

Preliminary diagnostics only; three replications per DGP are not substantive Monte Carlo evidence.

- Matched DGP hashes: PASS (12/12 semantic draws, all three methods equal).
- Attempt accounting: PASS (36 = 36 method-replication evaluations).
- Fixed rank: 11/12 successful; failures {'coefficient_bound_hit': 1}.
- Fixed target retention: point 198/216; inference 198/216.
- Selected 3e-6 pilot success: 9/12.
- Selected 3e-6 ranks: {"[1, 1, 1]": 2, "[1, 1, 2]": 4, "[2, 1, 0]": 1, "[2, 1, 1]": 2}.
- Selected 4e-6 pilot success: 9/12.
- Selected 4e-6 ranks: {"[1, 0, 0]": 1, "[1, 0, 1]": 1, "[1, 1, 1]": 6, "[1, 1, 2]": 1}.
- Broad-target split reuse: PASS (exactly four coefficient fits per successful replication).
- Fit-level runtime: UNAVAILABLE in these executed outputs; instrumentation was corrected only for future runs, and no extra evaluations were launched.
- Decision: NO-GO because both pilots completed only 9/12 and fit-runtime reporting was incomplete; no medium or production run was launched.

See the adjacent CSV and JSON files for lossless detail.
