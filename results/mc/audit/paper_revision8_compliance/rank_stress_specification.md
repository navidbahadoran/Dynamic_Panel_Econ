# Exact current rank-stress specification

This file records current code behavior only. It does not recommend changing it before the paper is made reproducible.

## Common construction

Let `E_Araw`, `E_B`, and `E_Hraw` denote the deterministic rank-one raw envelopes. For any block `M` with requested rank `r>=1`, the code begins with the baseline rank-one matrix `M^(1)` and adds

\[
M^{raw}=M^{(1)}+\sum_{j=2}^{r}s_j\ell_j q_j',
\]

where, independently for each added component and of the baseline primitives,

- `ell_j,i iid U[-sqrt(3),sqrt(3)]` over units;
- `q_j,t iid U[-sqrt(3),sqrt(3)]` over every burn-in and observed column;
- the added factor is iid over time, not an AR(1);
- current `component_strengths=(1,1)`, so `s_j=1` for the rank-two designs.

The entire expanded block is multiplied by one common deterministic scalar

\[
s_M(r)=\frac{E_M}{E_M+\sum_{j=2}^{r}3|s_j|}.
\]

There is no entrywise clipping. This common rescaling makes the expanded raw support no larger than the corresponding rank-one support envelope. Requested rank zero produces an identically zero block and uses reported rescale factor one.

The rank-one raw envelopes are:

- `E_Araw=(loading-A envelope)*(0.5+0.1*3)`;
- `E_B=(loading-B envelope)*(0.6+0.2*3)`;
- `E_Hraw=sqrt(3)*3=3sqrt(3)`.

For DGPs 1-3 the loading envelopes are `1+0.1sqrt(3)` for A and `1+0.4sqrt(3)` for B. For DGP4 they are `1.1+0.08sqrt(3)` for A and `1.2+0.25sqrt(3)` for B.

After expansion and common rescaling, A alone receives the realized full-horizon stability treatment

\[
c_a=\min\{1,0.85/\max_{i,t\in burn+observed}|A^{raw}_{it}|\},
\qquad A=c_a A^{raw}.
\]

B receives no further scaling. H becomes `H_0=c_xi*c_H*H_raw`. The deterministic reported final envelopes remain `C_A=0.85`, the applicable rank-one `C_B`, and `C_H=|c_xi*c_H|*3sqrt(3)` for positive ranks.

## True rank `(2,1,1)`

- A receives the additional independent bounded outer product and is commonly rescaled by
  `E_Araw/(E_Araw+3)` before the realized `c_a` stability scaling.
- B is the unmodified rank-one baseline (`s_B=1`).
- H is the unmodified rank-one baseline (`s_H=1`).
- Calibration is performed separately for every `DGP x N x T x rank vector` using actual stress raw draws.
- `c_H` is obtained analytically from population moments at `pi_H=0.30` for the actual H rank (here rank one).
- `c_xi` is found from deterministic common-random-number calibration to the configured pooled-R2 target (currently 0.65 in the frozen implementation).

## True rank `(1,0,2)`

- A is the unmodified rank-one baseline followed by the standard realized `c_a` stability scaling.
- B is identically zero.
- H receives one added independent bounded outer product and the common rescaling
  `E_Hraw/(E_Hraw+3)`.
- `c_H` is recalculated from population moments for the actual rescaled rank-two H. Specifically, its population observed-entry variance is
  `s_H(2)^2 * (mean_t[1-0.5^(2k_t)] + 1)`, where the `+1` is the variance of the iid added factor-loading term.
- Because B is zero, the scale calibration is declared unidentified: `r2_scale_identified=false`.
- The established normalization is exactly `c_xi=1`.
- No target pooled R2 is claimed. The pooled R2 induced at this normalization is calculated and reported.

## Paper details required for exact reproduction

The paper currently must add all of the following: component strengths; primitive distribution of every added loading and factor; iid-over-time status of added factors; independence from baseline primitives and across components/blocks; use of the full burn-plus-observed horizon; the envelope formula; blockwise common rescaling formula; absence of clipping; A's post-expansion realized stability scaling; separate calibration by DGP/cell/rank vector; analytical population definition of `c_H` for actual H rank; deterministic CRN seed and number of calibration draws; root-selection rule for `c_xi`; and the zero-B normalization/induced-R2 rule.

Exact random-number reproduction also requires the master seed (`20260807` in current configs), RNG family, semantic seed labels, draw order, group assignment (first half G1, second half G2), and the initial conditions documented separately.
