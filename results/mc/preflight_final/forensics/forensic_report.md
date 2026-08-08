# Final-preflight failure forensics

No new DGP draws, medium diagnostics, power, rank stress, or production runs were launched.

## Central diagnosis

All pre-IC artifacts match for 12/12 semantic replications. The same pilot failures occur under both IC multipliers: dgp1_N50_T50_r00002_truth1-1-1, dgp3_N50_T50_r00000_truth1-1-1, dgp3_N50_T50_r00002_truth1-1-1.
The nuclear path, six pilot routes, route objectives/validity, cap-pilot rank, initial candidate ranks, Q-hat values, and dimensions are identical. DGP 2 replication 2 has an expected post-IC local-completion difference: the 3e-6 path adds two neighbors after its initial IC choice, while 4e-6 does not. Thus there is no stochastic reproducibility bug, but the literal final candidate pools are not identical in that cell.

## Pilot failures

- `dgp1_N50_T50_r00002_truth1-1-1`: 6/6 routes fully valid, 10 stable route pairs, but rank-adaptive cap pilot objective stability failed: best-two normalized gap=0.109596. Primary classification: large disagreement between otherwise valid local optima.
- `dgp3_N50_T50_r00000_truth1-1-1`: 5/6 routes fully valid, 6 stable route pairs, but rank-adaptive cap pilot objective stability failed: best-two normalized gap=0.283875. Primary classification: large disagreement between otherwise valid local optima.
- `dgp3_N50_T50_r00002_truth1-1-1`: 5/6 routes fully valid, 6 stable route pairs, but rank-adaptive cap pilot objective stability failed: best-two normalized gap=0.0294083. Primary classification: large disagreement between otherwise valid local optima.

The failures are not caused by iteration caps, stationarity, or rank collapse. They are caused by materially different local objective basins among otherwise credible routes. DGP 3 replication 0 also has one bound-invalid route, but five valid routes remain, so that is not the pilot-level cause.

## Retention

- `4e-6`: 216 attempted -> 144 point retained -> 126 inference retained. Pilot failures remove 54, the cap hit removes 18, and unsupported selected-rank targets remove 18 at inference. Twenty split-fit diagnostic statuses affect three replications but remove zero records because the stored finite estimates/SEs were retained; this status/retention inconsistency is reported rather than hidden.
- `3e-6`: 216 attempted -> 36 point and inference retained. This is entirely three pilot failures (54 targets) plus seven cap hits (126 targets).

## Fixed-rank bound event

Only one fixed-rank start is configured. Its exact-seed replay converges and passes stationarity but moves outside B. No authorized alternative start establishes an interior comparable solution; the event remains a coefficient-bound failure.
The replay also confirms that supplied-rank ALS treats `coefficient_bound` as an acceptance diagnostic rather than projecting iterates into the box: the envelope rises from about 0.089 to 14.439 while B=9. This is an implementation gap relative to a literal box-constrained optimization, but changing it would change the numerical estimator and was not authorized here.

## Runtime instrumentation

All 302 replayed coefficient-fit rows have genuine positive runtimes: True. Each has fit identity, requested/realized rank, objective, iterations, convergence, stationarity, envelope, and start/route context: True. No historical runtimes were fabricated. Split-fit runtime plumbing is verified by tests; none of the four authorized failure replays reaches split inference.

## Recommendation

The smallest plausible next numerical action is deterministic polishing/continuation of the already-converged route endpoints before comparing outer-route objectives, using the unchanged loss, box, ranks, starts, and tolerances. Do not relax objective stability. This is a recommendation only and was not implemented.

Saved rank outcomes remain: `3e-6` exact 2/9, under 1/9, over 7/9, cap hits 7/9; `4e-6` exact 6/9, under 2/9, over 1/9, cap hits 1/9. The initial IC choice is over the same candidate values; only the subsequent local-completion pool differs in DGP 2 replication 2, without changing either reported winner.

`3e-6` is not recommended for further preflight. `4e-6` remains a candidate for a later independent validation, but is not selected or frozen as a production default.
