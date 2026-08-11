# Solver semantic-equivalence audit

## Mathematical problem

Both numerical paths minimize the unchanged loss

\[
 (2NT)^{-1}\left\|Y-\sum_m Z_m\odot U_mV_m'\right\|_F^2
\]

over width-four products with block rank at most four and `max_m ||U_m V_m'||_max <= 10`. No penalty, clipping, barrier surrogate, singular-value cutoff, factor floor or approximate statistical objective was added.

## Gauge and lower-rank behavior

Balanced factorization reconstructs each product through its thin SVD. Positive directions use `P sqrt(S)` and `Q sqrt(S)`; exact zero directions use a zero loading and an orthonormal time direction. Thus the product is unchanged while rank 0–4 remains available. Deterministic tests found maximum product errors of numerical roundoff and unchanged objective/envelope. Rank-0/1/2 products were represented at full width four and recovered their intended product ranks.

## Linear algebra and QPs

Interior column equilibration is the bijection `x=Dz`; the solution maps back to the original coordinate. Matrix-free Gauss-Newton uses the exact derivative of the same fitted-value map, zero damping, and a monotone line search on the same objective.

The constrained solver uses the same quadratic objective and inequalities `-B <= Cx <= B`. Internal equilibration maps back before commit and before the frozen factor-KKT calculation. Active-set equality systems use SVD least squares and deterministic constraint ordering. The active-box fixture achieved maximum coefficient 0.49999995 under `B=0.5`, violation below `1e-8`, and KKT below `1e-4`.

## Frozen acceptance gate

- interior product projected-gradient: `1e-6`;
- constrained factor-space KKT: `1e-4`;
- coefficient bound: `B=10` in the Revision-10 selector;
- best-two normalized objective agreement: `1e-6`;
- three maintained deterministic starts;
- finite objective, literal feasibility, terminal status and unresolved behavior unchanged.

The real deterministic width-four three-start fixture produced three valid starts and passed the unchanged objective gate. No rank result is created when the gate fails.

## Legacy regression

The new kernel is not called by ordinary supplied-rank estimation, split fits, final selected-rank post-refits or `revision9_ic`. On a deterministic fixture where the old solver converges, old and new product fits agree to numerical precision and solve the same objective. All existing supplied-rank and Revision-9 tests pass; historical evidence is untouched.

## Latent rollback repair

The isolated cap+1 path snapshots factors and objective before a sweep and restores both if the original objective increases. A forced deterministic test verifies the returned product and stored objective remain consistent. The legacy solver was not modified.
