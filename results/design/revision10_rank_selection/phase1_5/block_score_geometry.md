# Block score geometry and local quadratic gain

For a residual array e and perturbation D in block M, RSS expands exactly as

`RSS(Theta)-RSS(Theta+D)=2<S_M^raw(e),D>-Gamma_MM(D,D)`,

where `Gamma_MM(D,D)=sum_it z_it,M^2 D_it^2` and
`S_M^raw(e)=Z_M elementwise e`.

If the already fitted tangent coordinates are jointly reoptimized, replace the raw design by its
empirical fitted-value residualization against all existing tangent designs. The resulting Schur
complement curvature is `Gamma_M|T(D,D)=||X_M^eff(D)||_F^2`. For a rank-one normal direction
`D=a b'`, `||a||=||b||=1`, duality gives

`sup_(a,b) <S_M^eff,a b'> = ||P_U,perp S_M^eff P_V,perp||_op`.

Thus, whenever `Gamma_M|T(D,D)>=mu_M||D||_F^2`,

`Delta_M^loc,eff <= ||P_U,perp S_M^eff P_V,perp||_op^2/mu_M`.

This is a local upper bound. A global post-refit may rotate all singular spaces and move farther
than one tangent chart, so a uniform localization argument is required before applying it to the
post-refit RSS gain.

## A blocks

`S_A(ell),it=e_it y_i,t-ell`. Sequential exogeneity centers this product because lagged outcomes
belong to past panel information. Its scale depends on dynamic propagation, the conditional
second moments of lagged outcomes, and spatial covariance of the contemporaneous innovations.

## B blocks

`S_B(k),it=e_it x_it,k`. Current x is included in the sequential-exogeneity conditioning set, so
the true innovation product is centered. Its variance and curvature depend on the particular
covariate's conditional moments and dependence; different B blocks need not share a scale.

## H block

`S_H,it=e_it`. There is no regressor multiplier. Its variance geometry is the innovation's own
conditional spatial covariance and is generally different from both A and B.

A single residual MSE cannot standardize these three geometries under conditional
heteroskedasticity. A common scalar can only be a conservative maximum of valid block envelopes;
it is not rate-sharp block by block. Matrix-specific normalizers are necessary for the intended
design.
