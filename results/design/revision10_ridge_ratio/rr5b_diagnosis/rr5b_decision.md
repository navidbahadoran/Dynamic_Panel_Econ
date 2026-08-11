# RR5b decision

## Classification: RR5b-NUMERICAL-SOLVER

All 180 cap+1 pilots failed the frozen multistart acceptance gate. Of 540 starts, 53 were individually valid. In 177 replications fewer than two starts were valid because termination/stationarity failed; in the remaining three, the two valid objectives disagreed by `3.66e-4`, `7.41e-4`, and `4.78e-3`, all above `1e-6`. No objective, feasibility, box, rank-support, boundary, or nonfinite condition independently rejected a replication.

This is numerical pilot failure, not statistical ridge-ratio failure: normalized spectra and ratios were never computed, and no rank was selected. Code inspection found no acceptance-scale mismatch or implementation bug. The current cap+1 solver architecture does not make two starts satisfy the frozen numerical credibility gate and is computationally impractical for larger runs in its present form.

Paper consequence: the rank-selection theory is not contradicted by RR5. Category A applies—the numerical implementation requires engineering before empirical evidence can exercise the theory—and category D applies operationally to the current three-start cap+1 estimator architecture. No manuscript or scientific rule is changed in RR5b.
