# Semantic equivalence of the preferred solver

## Product and objective invariance

For any block factors (M=UV') and invertible (G), define

\[
 \widetilde U=UG,\qquad \widetilde V=VG^{-T}.
\]

Then

\[
 \widetilde U\widetilde V'=UG(VG^{-T})'=UGG^{-1}V'=UV'=M.
\]

QR renormalization is this transformation with (G=R^{-1}) when nonsingular. Balanced SVD refactorization directly reconstructs the same product: if (M=P\Sigma Q'), then `(P sqrt(Sigma))(Q sqrt(Sigma))'=M`. Zero padding adds columns whose outer products are zero. A dormant direction with one factor column exactly zero also adds the zero matrix. Thus the fitted coefficient products, rank and every entrywise box inequality are unchanged at the transformation point.

Because

\[
 L(\Theta)=\frac1{2NT}\left\|Y-\sum_m Z_m\odot M_m\right\|_F^2
\]

depends only on the products, every such refactorization leaves the normalized objective exactly unchanged in exact arithmetic. Implementations must verify equality to the established floating-point tolerance and reject a transformation that exceeds it.

## Feasible-set invariance

Every width-(d_m) product has rank at most (d_m), and every matrix of rank at most (d_m) has a width-(d_m) factorization. Reparameterization and zero padding neither increase the product rank nor remove a product direction. No singular-value threshold is used to declare or select rank. Hence the represented set remains

\[
 \{M_m:\operatorname{rank}(M_m)\le d_m,\ \|M_m\|_{\max}\le B\}.
\]

Internal subproblem equilibration `x=Dz` with nonsingular diagonal (D) is bijective. The transformed constraints `-B <= C D z <= B` map exactly to `-B <= C x <= B`, and the transformed objective `||X D z-y||^2` maps exactly to `||Xx-y||^2`. Mapping back before committing therefore changes neither the subproblem nor global feasible product.

An exact active-set solve and SLSQP target the same convex row/time QP. Changing the numerical algorithm, warm start, or cached factorization does not change its minimizer set. A monotone line search evaluates the same objective along directions within the same factorized feasible set; it adds no term to the criterion.

## Stationarity and KKT

The interior projected-gradient residual is a function of the coefficient product, the supplied cap tangent representation, and the data. Product-preserving refactorization therefore leaves it unchanged, up to deterministic floating-point evaluation error. Its frozen threshold remains `1e-6`.

For constrained coordinate subproblems, exact KKT zero is invariant under a nonsingular coordinate change: if `x=Dz`, gradients transform as `grad_z=D' grad_x`, constraints as `C_z=C D`, and the same multipliers make the transformed Lagrangian gradient zero if and only if the original one is zero. Feasibility and complementarity are also identical.

However, the *finite normalized Euclidean factor-KKT residual* is not invariant under an arbitrary non-orthogonal global gauge because Euclidean norms change under `G`. This matters at the frozen nonzero threshold `1e-4`. The preferred constrained solver therefore uses equilibration only internally, maps back to the maintained factor coordinates, and evaluates the existing residual there. It does not apply balanced non-orthogonal gauge transformations between constrained sweeps. Orthogonal sign/rotation changes are permitted because they preserve Euclidean norms and the product.

This distinction is mandatory: balanced global scaling is semantically exact for the mathematical estimator and exact KKT zero, but it is not guaranteed to preserve a finite coordinate-dependent diagnostic. It is consequently not part of the constrained baseline without a new paper-level clarification.

## Objective-stability gate

The three start identities and RNG streams are unchanged. Each start must still satisfy the same finite objective, convergence, path-specific residual, feasibility and envelope checks. The best two valid objectives are still compared as

\[
 |L_1-L_2|/\max(1,|L_1|)\le10^{-6}.
\]

Better solution of the same problem may change the numerical iterates and allow agreement; it does not relax or replace the gate.

## Implementation invariants

A future implementation is semantically acceptable only if tests establish:

1. identical initial start factors/products and seeds;
2. coefficient-product preservation for every gauge operation;
3. unchanged normalized objective and literal box;
4. width `cap+1` retained with no statistical threshold;
5. identical row/time QP matrices before and after reversible equilibration;
6. residuals evaluated by the frozen functions and thresholds;
7. no approximate singular values enter the ridge ratios;
8. unresolved pilots still produce no rank or primary inference.
