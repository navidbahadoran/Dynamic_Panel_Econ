# Rank zero and truth at the reporting cap

## True rank zero

If `r_M,0=0`, every pilot singular value through `bar_r_M+1` is post-truth. Hence

`lambda_hat_M,j=O_p(zeta_NT)=o_p(a_NT)`, `j=1,...,bar_r_M+1`.

For the fixed dimensionless anchor `lambda_hat_M,0=1`,

`R_M,0={lambda_hat_M,1+a_NT}/{1+a_NT}
       ={a_NT+o_p(a_NT)}/{1+a_NT}->_p0`.

For every `j=1,...,bar_r_M`,

`R_M,j={a_NT+o_p(a_NT)}/{a_NT+o_p(a_NT)}->_p1`.

Thus zero is the unique asymptotic minimum, with no hypothesis test, cutoff, or user-selected constant.
This directly covers the zero middle block in the maintained rank-stress truth `(1,0,2)`.

## True rank at the reporting cap

If `r_M,0=bar_r_M`, the denominator of `R_M,bar_r_M` contains the last genuine signal singular value.
The numerator uses `lambda_hat_M,bar_r_M+1` from a matrix allowed to have rank `bar_r_M+1`. Weyl and the
pilot rate make this extra value `O_p(zeta_NT)`, rather than setting it mechanically to zero. Consequently,

`R_M,bar_r_M=O_p(a_NT)->0`.

Every lower-index competitor remains bounded away from zero by the positive-rank argument. Selection at the
cap is therefore statistical. The extra pilot rank never expands the reported set
`{0,...,bar_r_M}`.

If a computed pilot numerically collapses below the cap+1 allowance, that fact is recorded as numerical rank;
it is not used as a statistical threshold. A failed validity or stability check invokes the frozen unresolved
status rather than a mechanical cap decision.
