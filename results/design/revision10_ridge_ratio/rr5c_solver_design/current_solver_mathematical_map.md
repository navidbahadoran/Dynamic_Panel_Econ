# Current cap+1 solver: mathematical map

## Frozen coefficient problem

Let the coefficient blocks be (M_m\in\mathbb R^{N\times T}), with regressors (Z_m), and let (q) be the number of blocks. The pilot solves

\[
 \min_{\{M_m\}} L(M)=\frac1{2NT}\left\|Y-\sum_{m=1}^q Z_m\odot M_m\right\|_F^2,
 \quad \operatorname{rank}(M_m)\le d_m=\bar r_m+1,
 \quad \max_m\|M_m\|_{\max}\le B.
\]

RR5 has three blocks, reporting caps `(3,3,3)`, pilot widths `d=(4,4,4)`, and `B=10`. Each block is represented as

\[
 M_m=U_mV_m',\qquad U_m\in\mathbb R^{N\times d_m},\quad
 V_m\in\mathbb R^{T\times d_m}.
\]

This fixed-width factorization represents rank *at most* (d_m); the pilot deliberately disables exact numerical-rank rejection.

## Starts and product reconstruction

There are exactly three maintained deterministic seeded starts. The first uses the supplied pilot seed; two integer seeds are drawn deterministically from its NumPy generator. With no initial coefficient product, every entry of every (U_m,V_m) is independent normal with scale `0.1`. If an initial product is supplied elsewhere, balanced square-root SVD factors are used and deficient widths receive seeded `1e-4` perturbations. The pilot itself supplies no initial coefficient product.

After every update, coefficient products are reconstructed only as `U_m @ V_m.T`, in block order `A1, B1, H` for the RR5 design. The normalized objective is always SSE divided by `2NT`.

## Interior fast path

For unit (i), stack (u_i=(u_{1i}',\ldots,u_{qi}')'). With factors fixed, solve the unconstrained least-squares problem

\[
 \min_{u_i}\frac12\|y_i-X_i(V)u_i\|_2^2,
 \quad X_i(V)=\left[\operatorname{diag}(z_{1i})V_1\;\cdots\;
 \operatorname{diag}(z_{qi})V_q\right].
\]

All units are updated in increasing index order using `numpy.linalg.lstsq(..., rcond=1e-10)`. Each block is then re-gauged by a reduced QR of its loading,

\[
 U_m=Q_mR_m,\qquad (U_m,V_m)\leftarrow(Q_m,V_mR_m'),
\]

which preserves (U_mV_m'). For time (t), stack (v_t) and solve

\[
 \min_{v_t}\frac12\|y_{:t}-X_t(U)v_t\|_2^2,
 \quad X_t(U)=\left[\operatorname{diag}(z_{1:t})U_1\;\cdots\;
 \operatorname{diag}(z_{q:t})U_q\right],
\]

again in increasing time order. In RR5 each joint row/time problem has 12 variables.

If the global objective increases beyond `1e-11*max(1, previous objective)`, the loop stops before appending the new value. The implementation does **not** restore the prior factors: it therefore returns the post-sweep product while `FitResult.objective` remains the preceding history value. RR5 did not serialize enough history to establish whether this branch fired. This is a latent product/objective bookkeeping and rollback defect, not an established RR5 failure cause. Otherwise the loop stops when relative objective change, with denominator `max(abs(previous),1e-14)`, is at most `1e-8`, or after 200 sweeps. The final interior residual is

\[
 \frac{\|P_{T(\Theta;d)}\nabla L(\Theta)\|_F}
 {\max(\|Y\|_F/\sqrt{NT},1)},
\]

where the product-space gradient is projected onto the rank-(d) tangent representation. A result is marked converged if objective stopping occurred or this residual is at most `1e-6`; pilot validity separately requires the residual itself to be at most `1e-6`. This path is retained only if the final unconstrained product envelope is strictly below `B-1e-8`.

## Constrained fallback

If the unconstrained envelope is not interior, the unconstrained coefficient products are uniformly scaled just inside the box and refactorized by square-root SVD. The same seed supplies any rank-deficiency perturbation.

With (V) fixed, each unit subproblem is the convex linearly constrained least-squares problem

\[
 \min_{u_i}\frac1{2T}\|X_i(V)u_i-y_i\|_2^2
 \quad\text{s.t.}\quad -B\mathbf1\le C(V)u_i\le B\mathbf1,
\]

where (C(V)=\operatorname{blockdiag}(V_1,\ldots,V_q)). With (U) fixed, each time subproblem is the analogous problem with divisor (N) and (C(U)=\operatorname{blockdiag}(U_1,\ldots,U_q)). Thus the constraints are exactly the entrywise product box for the row or time slice being updated.

SLSQP starts from the current row/time vector, uses analytic gradients, `ftol=1e-10`, and at most 200 iterations per subproblem. A candidate subproblem is accepted only if finite, feasible within `1e-8`, nonincreasing within numerical slack, and either SLSQP reports success or its normalized active-set KKT residual is at most `max(1e-6,sqrt(1e-10))=1e-5`. A rejected subproblem stops the fit. Loading rows are updated first, followed by the same blockwise QR gauge, then time rows. There is no QR after the time half-sweep.

The outer constrained loop stops only when relative global objective change, with denominator `max(1,abs(previous))`, is at most `1e-8` and the maximum row/time factor-space box KKT residual is at most `1e-4`, or after 200 sweeps. Final convergence requires no subproblem failure, product feasibility within `1e-8`, and KKT at most `1e-4`.

## Pilot multistart gate

Each start must have a finite objective, `fit.converged=True`, its path-appropriate frozen stationarity/KKT pass, and product envelope ratio at most `1+1e-8`. Exact factor/product rank is not required. Valid starts are sorted by objective; at least two are required and their normalized objective difference must be at most `1e-6`. Failure raises `rank_selection_numerically_unresolved` before spectra or ratios are computed.

## Nonidentifiability and conditioning

- **Rotation:** for any invertible (G_m), (U_mV_m'=(U_mG_m)(V_mG_m^{-T})').
- **Scaling:** diagonal (G_m) can make one factor arbitrarily large and the other arbitrarily small without changing the product.
- **Near-collinearity:** nearly dependent factor columns make the joint row/time designs poorly conditioned.
- **Collapsed columns:** a zero or tiny direction lies on a lower-rank stratum. If both sides are zero, ordinary alternating updates cannot reactivate it.
- **Rank-deficient factor Grams:** (U_m'U_m) or (V_m'V_m) can be singular at lower rank; redundant columns also create flat factor-space directions.
- **Cross-block confounding:** the joint 12-column row/time design may be nearly dependent across `A1`, `B1`, and `H`, even when each block is internally normalized.
- **Product box versus factor gauge:** the box constrains (U_mV_m'), not either factor. Equivalent gauges can therefore have very different subproblem conditioning while representing the same feasible product.

The existing QR removes loading-side scale after each loading half-sweep, but puts all triangular scale into the time factor. It does not balance singular scale between the two factors or remove cross-block/redundant-direction conditioning.
