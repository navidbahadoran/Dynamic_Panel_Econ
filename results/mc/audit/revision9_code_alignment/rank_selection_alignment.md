# Locked Revision-9 rank-selection alignment

## Baseline rule

The implementation retains exactly

\[
b_{NT}=(NT)^{1/(8+\eta)}\log(NT),\qquad
\kappa_{NT}=b_{NT}^2\log(NT)^{d_s+3},
\]

with baseline multiplier `1.0`. The information criterion is `log(Q_hat)+kappa_NT*d(r)/(NT)` with `Q_hat=SSE/(NT)` from the unpenalized constrained post-refit. The only sensitivity multipliers are now `0.5` and `2.0`.

The singular-value threshold is

\[
\tau_{NT}=\sqrt{NT}/\log(NT),
\]

with only `0.5` and `2.0` sensitivities. The baseline threshold multiplier is `1.0`.

## Screening and candidates

The nuclear path uses `gamma=0.8`, `epsilon=0.01`, hence `L=1+ceil(log(.01)/log(.8))=22`. Dense sensitivity uses `sqrt(0.8)`. Nuclear estimates remain screening proposals/warm starts only.

The cap pilot solves the same rank-at-most `(3,3,3)` and max-entry-10 objective as the paper. Candidate construction remains: threshold every path estimate and cap pilot; union their rank vectors; add valid one-coordinate neighbors within the cap lattice; remove duplicates; unpenalized constrained fixed-rank post-refit; apply IC only to valid post-refits; locally complete omitted neighbors until none improves the IC. No Cartesian exhaustive search is used by the baseline selector.

## Caps and cap inference

Production, medium, pilot, and rank-stress maintained configurations resolve to `(3,3,3)`. All true vectors `(1,1,1)`, `(2,1,1)`, and `(1,0,2)` lie inside. The selected-to-cap diagnostic uses the resolved caps. A selected cap hit returns `rank_at_cap` before target inference; the selected rank is not altered.

Historical `2,2,2`, `3e-6`, and `4e-6` files remain only in explicitly named preflight/replay artifacts and are not referenced by the canonical paper-validation command.
