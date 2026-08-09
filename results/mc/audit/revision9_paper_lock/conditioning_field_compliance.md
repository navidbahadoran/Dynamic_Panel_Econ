# Monte Carlo conditioning-field compliance

## Main DGPs

The current primitives admit the Revision-9 partition

\[
\mathcal C^{MC}=\sigma(g_a,g_b,g_h,f_x)\vee\mathcal D_N,
\]

with

\[
\mathcal D_N=\sigma(\lambda_a,\lambda_b,\lambda_h,\lambda_x,
\sigma_i,\sigma_{e,i},G_i),
\]

where the DGP4 bounded loading perturbations are included in `lambda_a` and `lambda_b` (or may be listed separately before those loadings are formed). The group assignment `G_i` is deterministic in the code but may be included.

Conditioning on these variables fixes `A_0`, `B_0`, `H_raw`, unit-specific scales, and the bounded design loadings. The realized stability factor `c_a` is measurable with respect to `lambda_a` and the `g_a` path; `c_H` and `c_xi` are frozen cell constants.

The time-varying Gaussian disturbance innovations `epsilon_it` and covariate innovations `e_it` remain random. In DGPs 2--4 the spatial disturbance recursion is a measurable transformation of the still-random `epsilon_it`; in DGPs 3--4 its lag enters x without changing the partition.

No main-DGP primitive in the user's CONDITION ON list is missing. `sigma_i` is generated in code through `sigma_i^2 ~ U(0.5,1.5)`, and similarly for `sigma_e,i`.

**Main DGP mapping: MATCH.**

## Rank-stress extension

For the rank-stress theorem/design mapping, the stated field must additionally include each added bounded time-factor path (`tilde f_A` or `tilde g_H`) among the common paths and each added bounded unit loading (`tilde lambda_A` or `tilde lambda_H`) in `D_N`. These primitives are absent from the supplied main-DGP list but are required to make the higher-rank true matrices fixed conditionally.

**Rank-stress mapping without this explicit extension: PAPER UNDERSPECIFIED.**
