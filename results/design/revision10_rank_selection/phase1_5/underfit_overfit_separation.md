# Underfit, overfit, zero-rank, and cap separation

## Underfit

For any profile class with `s<r_M,0`, Eckart-Young and the manuscript's strong-signal assumption
give

`inf_(rank(M)<=s) ||M-M_0||_F^2 >= sigma_(s+1)(M_0)^2 >= c NT`.

Allowing all other blocks to their caps cannot remove this coefficient distance. Conditional
prediction identification therefore gives an order-NT population loss gap for the class. The
existing uniform score and empirical-curvature bounds make stochastic remainders
`O_p(NT zeta_NT)`.

What does not follow is an order-NT difference between *successive optimized profile classes*.
The precise second new condition is:

`min_M min_(0<=s<r_M,0) Delta_M^prof,pop(s)/(NT) >= c_gain>0`.

Under it, `Delta_M^prof(s)=c_M,s NT+o_p(NT)`. Any valid envelope satisfying
`max Ehat_M(s)=o_p(NT)` gives `Delta/Ehat -> infinity` for every missing rank.

## Overfit

For `s>=r_M,0`, truth is feasible in both profiled classes. The manuscript's basic inequality,
uniform score bound, and curvature imply

`Delta_M^prof(s)=O_p(NT zeta_NT)`

uniformly over the fixed L comparisons. A valid observable envelope must dominate these gains
simultaneously with conditional probability tending to one. An `O_p` rate alone cannot specify
the boundary constant.

## True rank zero

If `r_M,0=0`, the comparison at s=0 is already an unnecessary increment. Nuisance-at-cap
profiling makes the full truth feasible under the null profile. The same simultaneous envelope
therefore yields

`P(Delta_M^prof(0)>Ehat_M(0)|C)->_p 0`.

Hence no first rank is added, including for the middle zero block in truth `(1,0,2)`.

## Caps

The theorem may retain `r_0 in R_max`. If `r_M,0<cap_M`, the first null increment stops block M.
If `r_M,0=cap_M`, all lower increments are detected and the block ends at its cap without a null
comparison above it. If truth lies outside a cap, exact recovery is impossible and the procedure
should report cap truncation; no automatic cap expansion is defined.
