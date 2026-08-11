# Blockwise ratio separation

Fix a block and write `r=r_M,0`, `bar_r=bar_r_M`. On events whose conditional probability tends to one,
the normalized positive spectrum satisfies

`underline_lambda <= lambda_hat_j <= overline_lambda`, for `1<=j<=r`,

for fixed `0<underline_lambda<=overline_lambda<infinity`, while

`max_{j>r} lambda_hat_j=O_p(zeta_NT)=o_p(a_NT)`.

The bounds hold jointly because the block set and all caps are fixed.

## Positive truth

Suppose `0<r<=bar_r`.

Before truth, for `1<=j<r` and eventually `a_NT<=1`,

`R_j=(lambda_hat_j+1+a_NT)/(lambda_hat_j+a_NT)
    >=underline_lambda/(overline_lambda+1)`.

For the anchored competitor,

`R_0=(lambda_hat_1+a_NT)/(1+a_NT)>=underline_lambda/2`.

Thus every pre-truth competitor is bounded away from zero by an explicit fixed constant.

At truth,

`R_r={O_p(zeta_NT)+a_NT}/{Theta_p(1)+a_NT}=O_p(a_NT)->0`.

After truth, for `r<j<=bar_r`, both adjacent spectral values are `o_p(a_NT)`, hence

`R_j={a_NT+o_p(a_NT)}/{a_NT+o_p(a_NT)}->_p1`.

Therefore one may take an asymptotic competitor separation constant smaller than

`c_* = min{underline_lambda/2, underline_lambda/(overline_lambda+1), 1/2}`.

With probability tending to one, `R_r<c_*/2` and every `R_j`, `j!=r`, exceeds `c_*/2`. The truth is the
unique minimizer; the tie rule is irrelevant asymptotically.

## Dependence and scaling

All statements are conditional on a regular common-shock realization. Dependence, predetermined regressors,
and rectangular dimensions enter only through the generic pilot rate and maintained scale-factor bounds.
There is no IC, penalty dominance, candidate coverage, or additional spectral-gap condition.
