# Scale and unit equivariance

## Audit of the alternatives

### A. Raw coefficient spectrum — rejected

The RR1 provisional object `sigma_j(M_hat)^2/(NT)` is not a unit-equivariant reported-rank input. If
`x_k^*=c x_k`, the equivalent coefficient is `B_k^*=B_k/c`, so its raw squared singular values are divided
by `c^2`, while the fixed ridge is unchanged. Consistency still holds for each fixed nonzero `c`, but the
finite ratio and its margin depend arbitrarily on the choice of covariate units.

### B. Existing fitted-value weights — selected, with a common reference

The manuscript already defines `w_A,l`, `w_B,k`, and `w_H=1` to put blocks on fitted-value scales. The
product `w_M M` has outcome units for every block, but `w_M^2 sigma_j(M)^2/(NT)` is not dimensionless.
Dividing every block by the common existing outcome-unit reference `s_NT=w_A,1` resolves that issue.

### C. Coefficient envelope — rejected

The common numerical envelope `B` is a coordinate constraint, not a block unit. Dividing by `B` would not
cancel `B_k -> B_k/c` unless the envelope were made block- and unit-specific, which would alter the
maintained coefficient-bound formulation.

### D. Final normalization — frozen

For `q_M=w_M/w_A,1`, freeze

`lambda_hat_M,j=q_M^2 sigma_j(M_hat^pil)^2/(NT)`.

It has no tuning constant, uses only existing preprocessing quantities, and leaves rank unchanged whenever
`q_M>0`.

## Algebraic equivariance

Let the outcome be expressed in units `y^*=d y` and covariate `x_k` in units `x_k^*=c_k x_k`, where all
constants are fixed and nonzero. The equivalent coefficients satisfy

`A_l^*=A_l`, `B_k^*=(d/c_k)B_k`, and `H^*=dH`.

The weights satisfy

`w_A,l^*=|d|w_A,l`, `w_B,k^*=|c_k|w_B,k`, `w_H^*=1`, and
`s_NT^*=|d|s_NT`.

Consequently, each transformed scaled block equals its original scaled block up to a sign:

`q_A,l^* A_l^*=q_A,l A_l`,
`q_B,k^* B_k^*=sign(d/c_k) q_B,k B_k`, and
`q_H^* H^*=sign(d)q_H H`.

Singular values remove the sign. Every `lambda_hat`, every ratio, and the selected rank are therefore
exactly invariant whenever the corresponding pilot solution is reparameterized coherently.

The literal box constraint deserves a precise qualification. A fixed numerical condition
`||Theta||_max<=B` is itself unit-dependent. Rescaling a regressor but holding that numerical feasible set
unchanged defines a different constrained statistical problem; no spectral normalization can make those two
optimization problems identical. Exact finite optimization equivariance therefore means either (i) map the
box under the same unit reparameterization, or (ii) convert inputs back to the stored canonical scientific
units before applying the unchanged `B`. The Revision-10 protocol freezes the latter interpretation. Under
the maintained strict interior condition, fixed coherent unit representations also have the same asymptotic
rank conclusion. Box activity is always reported as a numerical/theorem diagnostic.

## Stochastic orders under existing assumptions

The maintained moment/concentration bounds give upper stochastic bounds for the empirical second moments.
For the reference weight, sequential exogeneity and the conditional innovation-variance lower bound give
`E_0(y_it^2|x_it,F_t-1,C)>=Var_0(u_it|x_it,F_t-1,C)>=c`; stability and the moment assumptions give the upper
bound. Removing finitely many presample boundary observations does not change the order, and the maintained
concentration argument transfers these bounds to `w_A,1`.

Prediction identification gives the corresponding lower bound for each selectable regressor block. If its
truth is positive rank, take an admissible difference that scales only that true block by a small fixed amount.
Then identification and `||M_0||_F^2>=cNT`, together with `|M_0,it|<=B`, imply
`E_0 sum_it z_M,it^2>=cNT`. If its truth is zero and its reporting cap is positive, use a small constant rank-one
alternative. (A block with reporting cap zero has known rank and is omitted from selection.) Empirical
concentration again yields the sample lower bound. Thus every nontrivial existing weight is nondegenerate.
Because the block collection is fixed,

`0<c_q<=min_M q_M<=max_M q_M<=C_q<infinity`

with probability tending to one, conditionally on a regular common-shock realization. These bounds are a
consequence of the maintained assumptions, not a new DGP condition.

If the generic pilot rate is `max_M ||M_hat^pil-M_0||_op=O_p(sqrt(NT zeta_NT))`, Weyl and the strong-signal
condition yield, uniformly over fixed blocks and indices,

- `lambda_hat_M,j=Theta_p(1)`, bounded away from zero, for `j<=r_M,0`;
- `lambda_hat_M,j=O_p(zeta_NT)` for `j>r_M,0`.

Thus the new normalization preserves exactly the signal/noise orders used in RR1.

## Anchor audit

After normalization, `lambda_hat_M,j` is dimensionless, so the deterministic anchor
`lambda_hat_M,0=1` is coherent across blocks and units. It is a convention, not a cutoff. For rank zero its
ratio tends to zero; for positive rank its numerator is bounded away from zero. No alternate anchor or
selected constant is needed.
