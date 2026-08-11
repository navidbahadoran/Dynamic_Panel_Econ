# Contribution after removing headline automatic rank selection

The paper remains coherent and potentially strong without an automatic selector. Its unified contribution is inference for low-rank dynamic-panel coefficient matrices when different structural blocks have different ranks and latent spaces.

The contribution package is:

- joint estimation of separately ranked autoregressive, covariate-effect, and interactive-effect matrices;
- target-specific empirical tangent-space Riesz representation;
- full-panel inference for local and fixed-time targets;
- two-way split correction for broad averages;
- conditional spatial dependence with the paper's spatial variance construction;
- distinct ranks and latent spaces rather than a single common factor dimension;
- rectangular `N,T` asymptotics under the exact `1/N+1/T` growth condition.

These results do not logically require automatic rank estimation. The supplied-rank condition must be prominent, and the paper must explain specification sensitivity, but this is a narrower and more credible limitation than presenting an unusable selector.

Contribution language should avoid “adaptive,” “data-driven rank,” and “rank-robust.” The paper should say that it develops inference conditional on low-rank specifications and studies sensitivity across prespecified nearby ranks.

Automatic-rank literature should become context rather than a headline contribution. The core novelty remains the inferential architecture across several separately low-rank dynamic-panel matrices.
