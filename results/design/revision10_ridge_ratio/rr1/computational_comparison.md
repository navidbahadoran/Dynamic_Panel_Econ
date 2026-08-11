# Computational and numerical implications

## Workload

Revision 9 statistically selects ranks using a nuclear path, one cap pilot, thresholded candidate
construction, many unpenalized fixed-rank candidate post-refits, an IC comparison, and local
one-coordinate completion.

The provisional ridge-ratio route would require:

1. one joint rank-at-most-cap+1 pilot;
2. singular values through `bar_r_M+1` for every block;
3. deterministic blockwise ratios and an argmin;
4. one final unpenalized exact-selected-rank joint post-refit.

The nuclear path could remain as initialization and a diagnostic source of singular spaces. The old cap
pilot could be replaced by the cap+1 spectral pilot. Thresholded path ranks, candidate construction,
candidate post-refit enumeration, the final IC, penalty sensitivities, and IC local completion would no
longer be statistically necessary.

## Numerical condition

Let `L_NT` be the manuscript's normalized half-squared loss and let `L_NT^*(cap+1)` be the global
infimum over the enlarged fixed cap and literal box constraint. The computed spectral pilot must satisfy

`L_NT(Theta_tilde^(cap+1))-L_NT^*(cap+1)=o_p(zeta_NT)`.

This is the same form used by the existing cap-pilot proof; no fixed-rank numerical tolerance changes.
The condition contains no unknown true coefficient or rank, but the exact global infimum is not directly
observable for a nonconvex fit. Implementation can report all start objectives, best-start gaps,
stationarity, constraint activity, numerical rank, and objective stability without knowing truth. Those are
credible diagnostics, not a proof of the global objective gap. A certified optimization lower bound would
make the condition directly checkable; RR1 does not introduce one.
