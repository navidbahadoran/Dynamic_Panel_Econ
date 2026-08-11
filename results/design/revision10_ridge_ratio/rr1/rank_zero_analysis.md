# Rank-zero extension

Introduce the deterministic normalized anchor

`lambda_hat_M,0=1`.

Then the entire ratio family can be written

`R_M,j=(lambda_hat_M,j+1+a_NT)/(lambda_hat_M,j+a_NT),
        j=0,...,bar_r_M`,

where `a_NT=1/log(NT)`. For j=0 this is exactly

`R_M,0=(lambda_hat_M,1+a_NT)/(1+a_NT)`.

If `r_M,0=0`, every pilot singular value through `bar_r_M+1` obeys
`lambda_hat_M,j=O_p(zeta_NT)=o_p(a_NT)`. Therefore

`R_M,0=[a_NT+o_p(a_NT)]/[1+a_NT] ->_p 0`,

whereas, for every `j>=1`,

`R_M,j=[a_NT+o_p(a_NT)]/[a_NT+o_p(a_NT)] ->_p 1`.

Thus rank zero is the unique asymptotic ratio minimum. The anchor 1 is a normalization convention,
not a fitted threshold. Any fixed positive anchor would have the same limit, so fixing one removes
rather than introduces a tuning choice. This handles the zero component in truth `(1,0,2)`.

For positive truth, `lambda_hat_M,1` is bounded below with probability tending to one, so `R_M,0`
is bounded away from zero and cannot compete with the vanishing true-rank ratio.
