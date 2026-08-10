# Formal statement of the Phase-1 ST2 gap

This analysis uses the external author-supplied manuscript
`E:\OneDrive\Desktop\ver7_revision9_Montecarlo_appendix_design.tex`, SHA-256
`2525d019ec5d7b28457585ff57c44670d7fafc6c199c60acc753a76e85070138`.
The prompt's `(1)` duplicate-download name is not present; no older repository manuscript was
substituted.

Let `X_M(D)` be the N by T fitted-value array with entries
`[X_M(D)]_it=z_it,M D_it`, where `z_it,A(ell)=y_i,t-ell`,
`z_it,B(k)=x_it,k`, and `z_it,H=1`. For an exact rank-r post-refit let
`e(r)=Y-Y(Theta_hat(r))` and define:

- raw score matrix: `S_M^raw(r)=X_M^* e(r)=Z_M elementwise e(r)`;
- normal score: `S_M^N(r)=P_U,M,perp S_M^raw(r) P_V,M,perp`;
- operator score: `T_M(r)=||S_M^N(r)||_op`;
- local rank-one gain:
  `Delta_M^loc(r)=sup_D {2<e(r),X_M(D)>-||X_M(D)||_F^2}`, where D is an
  admissible rank-one normal direction;
- post-refit gain:
  `Delta_M(r)=RSS(r)-RSS(r+e_M)`.

The raw score is a matrix, its operator norm is a scalar dual norm, the local gain is a
quadratic optimization, and the post-refit gain allows all existing factors/loadings to change.
They are not interchangeable. Because a local normal update is feasible in the enlarged model,
`Delta_M(r) >= Delta_M^loc(r)` for exact global post-refits. The reverse inequality requires
profiling, curvature, and localization; it does not follow from first-order conditions.

Let `T_r` be the joint tangent space of all currently fitted blocks and let `P_XT` project in
fitted-value norm onto `X(T_r)`. Define the residualized design
`X_M^eff(D)=(I-P_XT)X_M(D)`, the effective score by
`<S_M^eff,D>=<e(r),X_M^eff(D)>`, and

`mu_M(r)=inf_{D in C_M(r), ||D||_F=1} ||X_M^eff(D)||_F^2`.

If `mu_M(r)>0`, exact quadratic maximization gives

`Delta_M^loc,eff(r) <= ||P_U,perp S_M^eff(r) P_V,perp||_op^2/mu_M(r)`.

The constant is exactly one under the RSS convention: `sup_x {2 a x-mu x^2}=a^2/mu`.
With absolute RSS suboptimality at most `epsilon_RSS` for each endpoint, an observed gain differs
from its ideal endpoint gain by at most `2 epsilon_RSS`.

The exact missing probabilistic object is an observable envelope `Ehat_M(r)` satisfying,
conditionally on the common shocks,

`P_0(max_(r,M in N_0) {Delta_M^ideal(r)-Ehat_M(r)}>0 | C) ->_p 0`,

where `N_0` includes every unnecessary increment, including a block already at truth while
another block is underfit. It must also obey `max Ehat_M(r)=o_p(NT)` so a genuine order-NT gain
dominates it. The finite caps make the displayed maximum finite, but Revision 9 proves only
stochastic orders with unknown constants, not this observable simultaneous quantile.
