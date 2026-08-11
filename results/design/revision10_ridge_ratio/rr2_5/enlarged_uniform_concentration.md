# Uniform concentration on the enlarged pilot class

Let `n=NT`, `m=N+T`, `b=b_NT`, and let `s_M^+=bar_r_M+1+r_M,0` and
`S_+=sum_M s_M^+`. Assume the localized identification extension on `D_max^+`.

## Score

The manuscript's localization, truncation, blocking, and matrix Bernstein argument is independent of the
coefficient ranks. It continues to give for every fixed block

`max_M ||S_M||_op
 =O_p[sqrt((N+T)log n)+b^2(log n)^((d_s+3)/2)]
 =O_p(sqrt(n zeta_NT))`.

For `Delta in D_max^+`, nuclear/operator duality and
`||Delta_M||_*<=sqrt(s_M^+)||Delta_M||_F` yield

`|sum_it u_it Y_it(Delta)|
 <=max_M ||S_M||_op sum_M ||Delta_M||_*
 <=sqrt(S_+) max_M ||S_M||_op ||Delta||_F`.

Since `S_+` is fixed,

`sup_(Delta in D_max^+\{0})
 |sum_it u_it Y_it(Delta)|/||Delta||_F
 =O_p(sqrt(NT zeta_NT))`.

Relative to Revision 9, only the fixed multiplier changes from its reporting-class counterpart to
`sqrt(S_+)`.

## Quadratic-form net

On a dyadic Frobenius range, a standard singular-vector parameterization for rank at most `s_M^+` gives an
`epsilon`-net satisfying, for a universal numerical constant `C`,

`log |N_r^+|
 <=C [sum_M s_M^+] (N+T+1) log(C/epsilon)
 <=C_+ m log(Cb)`

when `epsilon=c_0/b^2`. Here `C_+` depends only on the fixed block count and the fixed enlarged ranks. It has
no `N,T` dependence. The finite number `|J|^2` of regressor-block pairs is unchanged.

The entrywise box still gives `max_itM |delta_M,it|<=2B/r` for a normalized net point in a range with lower
endpoint `r`. Hence the quadratic weights retain exactly the manuscript bounds

`sum_it |w_it^(MM')|<=1`,
`sum_it (w_it^(MM'))^2<=C/r^2`, and `max_it |w_it^(MM')|<=C/r^2`.

Set

`r_0^2=C_0 n zeta_NT=C_0 b^2(N+T)(log n)^(d_s+2)`.

Choose the fixed `C_0` large enough as a function of `C_+`, `|J|^2`, `B`, and the unchanged identification
constant. The Bernstein exponent remains at least `c C_0 m log n`, so it dominates the enlarged net entropy,
the finite block-pair union, and the `O(log n)` dyadic ranges. The regressor multiplication map retains its
`Cb` Lipschitz bound, which is rank-free. Thus the net argument extends with only fixed constants.

Using the enlarged population identification constant `c_+` gives, uniformly on `D_max^+`,

`sum_it Y_it(Delta)^2
 >=(c_+/2)||Delta||_F^2-CNT zeta_NT`.

The small-norm case is absorbed by the additive term exactly as in Revision 9; the large-norm case follows
from population curvature and the uniform centered quadratic bound.

## Rate conclusion

The enlarged ranks change `sqrt(S_+)`, `C_+`, and the sufficiently large fixed peeling constant `C_0` only.
They do not change `b_NT`, `zeta_NT`, any logarithmic exponent, or the growth condition. The tangent-space Gram
part of `prop:uniform_concentration` is not enlarged and remains the supplied-rank result on `T_0`.
