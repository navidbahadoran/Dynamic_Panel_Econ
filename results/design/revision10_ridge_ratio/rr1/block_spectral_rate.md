# Block spectral rate from the existing pilot theorem

Fix block M and let `E_M=M_hat^pil-M_0`. The enlarged-cap proof gives

`max_M ||E_M||_op=O_p(sqrt(NT zeta_NT))`.

Weyl's inequality implies, simultaneously for the fixed collection of blocks and fixed indices,

`|sigma_j(M_hat^pil)-sigma_j(M_0)| <= ||E_M||_op`.

Define `lambda_hat_M,j=sigma_j(M_hat^pil)^2/(NT)`.

## Positive singular values

For `j<=r_M,0`, the manuscript signal condition gives

`c <= sigma_j(M_0)/sqrt(NT) <= C`.

Since `zeta_NT->0`, Weyl yields

`sigma_j(M_hat^pil)/sqrt(NT)
 =sigma_j(M_0)/sqrt(NT)+O_p(sqrt(zeta_NT))`.

Consequently, with probability tending to one, `lambda_hat_M,j` is bounded above and below by
fixed positive constants: `lambda_hat_M,j=Theta_p(1)`.

## Post-rank singular values

For `j>r_M,0`, `sigma_j(M_0)=0`, so

`lambda_hat_M,j
 <= ||E_M||_op^2/(NT)
 =O_p(zeta_NT)`.

Because the block count and caps are fixed, both conclusions hold jointly by a finite union. No
new spectral-rate assumption is required, and the proof remains conditional on the common shocks.
