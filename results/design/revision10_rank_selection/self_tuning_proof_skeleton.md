# Self-tuning proof skeleton

## Target theorem

For the eventual, fully specified selector `r_hat_ST`, prove

`P(r_hat_ST=r_0 | C) ->_p 1`.

Then on the rank-recovery event invoke the existing supplied-rank recovery, target expansion, two-way correction, Riesz, spatial-variance, and selected-rank inference results without altering their proofs.

## Proposed state machine

Initialize `r^(0)=0`. At step t evaluate all feasible `r^(t)+e_M`, using unpenalized box-constrained joint post-refits. Form matrix-specific G_M statistics. If none passes the ultimately proved envelope boundary, stop. Otherwise add the coordinate with largest standardized excess; break exact ties by the fixed manuscript block order. Ranks only rise and caps bound the number of steps. This rule is deterministic conditional on post-refits and handles zero rank.

## Lemmas needed

1. **Dimension algebra.** `Delta_d_M=N+T-2r_M-1>0` below a valid cap.
2. **Numerical fidelity.** Uniform post-refit objective error is negligible relative to local pass/fail margins; a numerically unresolved neighbor cannot trigger an addition.
3. **Overfit envelope.** At every state containing truth coordinatewise, every extra-rank optimized RSS gain divided by Delta_d is no larger than the block self-normalizer, simultaneously with conditional probability tending to one. Existing finite caps, matrix count, score bounds, and `O_p(zeta_NT)` improvement localize the task, but do not alone provide a tuning-free critical envelope.
4. **Underfit step.** At every reachable underfit state, prediction identification and a missing order-sqrt(NT) singular direction yield at least one positive population neighbor gain of order NT; the stochastic and numerical remainders are lower order, so its G statistic crosses the boundary.
5. **No false step before truth.** The simultaneous envelope prevents adding a coordinate beyond its true block rank even while another block is underfit. This uniform mixed-state version must be included in the self-normalizer lemma.
6. **Finite induction.** Lemmas 4-5 imply each step adds a missing true coordinate and never a false one. After `sum_M r_M,0` steps the state is truth; Lemma 3 makes it stop.

## Existing theory reused

The manuscript already gives order-one normalized loss separation for underfit models, `O_p(zeta_NT)` overfit improvements, strong singular values, prediction identification, uniform concentration, finite caps, and joint fixed-rank numerical requirements. These establish the rates but not the exact observable self-normalizer or hereditary neighbor property.

## Two new primitive conditions (maximum)

1. Uniform block self-normalizer validity, including mixed under/over states and simultaneous block/cap control.
2. Hereditary one-coordinate detectability along the upward lattice.

Neither condition may impose iid/Gaussian errors, cross-sectional or serial independence, homoskedasticity, exogenous future covariates, or stronger-than-order-sqrt(NT) signals. If deriving them requires any such restriction, the design becomes **THEORY COST TOO HIGH**.

## Open proof gap

An ordinary spatial-HAC variance estimates a variance, not the upper tail of an optimized operator-norm score. A valid observable, tuning-free envelope under conditional mixing/NED has not yet been derived. The boundary also must vanish in false-addition probability without inserting an outcome-chosen constant. Until that gap and the mixed-state induction are resolved, theorem drafting is not authorized.
