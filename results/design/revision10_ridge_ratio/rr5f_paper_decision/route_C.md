# Route C: replace the cap+1 pilot

Route C is scientifically possible but is a new theory/design phase, not a repair authorized by RR5f. The project has studied the following candidate classes.

## Nuclear/full-spectrum pilot

This object is computationally available from the existing nuclear path and supplies more than `bar_r+1` singular values. Revision-9 diagnostics found that at least one nuclear-path proposal equaled truth in 24/24 preflight draws, so the path was useful for screening. That evidence does not establish a final selector, an operator-rate theorem for the penalized estimator, or unbiased recovery of its noise spectrum.

Required new result: a uniform block operator-norm rate for the chosen nuclear/full-spectrum estimator under the paper's dynamic-panel dependence, coefficient box, and rectangular `N,T` growth, followed by a perturbation/ridge-ratio consistency theorem for its biased spectrum.

Potential new assumptions and tuning: penalty-level control, approximation/bias conditions, and a rule fixing which path point defines the spectrum. Dependence enters the penalized score and residual replacement. Rectangular theory must preserve the paper's unrestricted `1/N+1/T` condition rather than import balance restrictions. Computation is likely tractable because the path already exists, but statistical definition and tuning are unresolved. Prior diagnostics encouraged it as screening only and explicitly did not authorize promotion.

## Residual-score spectrum

Required new result: define a genuine post-cap noise spectral object, prove its normalization and uniform perturbation rate, and show that its estimated extra eigenvalue separates signal and noise for all three blocks, including rank zero and truth at cap.

Potential new assumptions and tuning: residual replacement, operator-valued conditional variance control, and possibly a spectral regularization or truncation rule. Conditional spatial dependence, lagged outcomes, heteroskedastic block scores, and rectangular panels complicate concentration. The prior extra-spectrum audit classified this as a new statistical construction with its own proof; it was not implemented or validated.

## Other documented extra-spectrum constructions

The prior audit found no existing manuscript object that simultaneously provides a genuine `(bar_r_M+1)`st value and the required generic operator rate. Any alternative would need its own estimator definition, rate theorem, normalization, dependence analysis, rectangular proof, tuning freeze, implementation audit, and independent validation. RR5f does not invent one.

## Recommendation

Not recommended for the present paper timeline. Route C may be a separate research project, but it has materially higher theory, tuning, implementation, and validation cost than repositioning the current paper around supplied-rank inference.
