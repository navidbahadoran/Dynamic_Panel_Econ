# Revision-10 empirical rank diagnostics

The empirical application should report:

- reporting caps and cap+1 spectral-pilot caps;
- normalized pilot singular values through each block's cap+1 value;
- `R_{M,0},...,R_{M,bar r_M}`, selected block ranks, and smallest-versus-second-smallest ratio margins;
- pilot feasibility, multistart objectives, stationarity residuals, objective stability, box activity, and termination status;
- final fixed-rank post-refit feasibility, stationarity, objective stability, and boundary/interiority diagnostics.

If the pilot fails the maintained numerical acceptance criteria, the status is `rank_selection_numerically_unresolved`; no fallback rank, primary selected-rank point estimate, or selected-rank inference is reported. Neighboring supplied-rank estimates may remain separate sensitivity diagnostics but are not part of the selector and do not imply rank-robust inference.

Candidate counts, candidate coverage, IC values or gaps, penalty/threshold multipliers, and local completion are not Revision-10 empirical diagnostics.
