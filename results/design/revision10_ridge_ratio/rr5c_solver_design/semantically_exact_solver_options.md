# Semantically exact cap+1 solver options

All options below retain the same product problem, widths, box, normalized objective, three maintained starts, and frozen acceptance gate. They change only how the same finite-dimensional problem is represented or solved.

| Option | Exactness argument | Conditioning effect | Complexity | Expected runtime effect |
|---|---|---|---|---|
| Balanced SVD gauge after each accepted half-sweep | Replace `M=UV'` by thin SVD factors `U*=P sqrt(S)`, `V*=Q sqrt(S)` with zero padding to width `d`; product, rank and box are identical | Shares singular scale and orthogonalizes both sides | Moderate | Extra small/block SVD; likely fewer ill-conditioned sweeps |
| Product-preserving QR/polar gauge | Apply invertible `G` and inverse transpose to the paired factor | Product is algebraically unchanged | Controls column norms without reconstructing full products | Low | Small overhead; less robust at deficient rank than SVD |
| Deterministic dormant-direction handling | At a zero singular direction, keep one side zero and install a deterministic orthonormal direction on the other; product remains exactly zero and width remains available | Avoids two-sided-zero absorbing states | Moderate and delicate | May reduce stalls at lower-rank strata |
| Rank-revealing coordinate basis, without statistical thresholding | Use numerical pivoting only to choose a factor coordinate basis; never delete a feasible product direction or use it for rank selection | Reduces redundant coordinate columns | Moderate | Faster/stabler solves; must pass product-invariance tests |
| QR/SVD least-squares solves | Solve the same linear least-squares subproblem without forming normal equations | Better stability for nearly dependent columns | Low; NumPy `lstsq` already uses a stable backend | Interior benefit may be limited; explicit diagnostics useful |
| Exact convex-QP active-set solver for box-product subproblems | Solve the identical quadratic objective under identical linear inequalities to KKT accuracy | Exploits small dense QP structure and warm active sets | High | Expected major fallback speedup versus generic SLSQP |
| Exact variable equilibration/preconditioning | Apply invertible `x=Dz`, transform design and constraints, then map back | Reduces scale disparity; feasible `x` set is unchanged | Moderate | Fewer QP/LS iterations |
| Warm active-set and primal iterate | Starts each unchanged convex subproblem from the previous solution/active set | Improves sequential sweep continuity | Moderate | Likely substantial constrained speedup; current code warm-starts only the primal vector |
| Cache per-half-sweep invariant constraint blocks and symbolic layouts | Reuse identical `C(V)` for all loading rows and `C(U)` for all time rows; numeric outcomes unchanged | No estimator effect | Low | Avoids repeated assembly/factorization overhead |
| Reusable work arrays | Storage reuse only | No estimator effect | Low | Reduces allocation overhead |
| Monotone global safeguard | Commit an update/half-sweep only when the same objective does not increase beyond roundoff; otherwise backtrack on the same factor-coordinate line | Prevents loss-increasing numerical steps | Moderate | May avoid wasted sweeps; line search adds evaluations |
| Deterministic line search | Evaluate the frozen objective along a factor-coordinate segment and choose a nonincreasing step | Same feasible factor path and objective; no penalty is added | Moderate | Helps unstable steps; should be inactive for exact coordinate minimizers |
| Recompute frozen residuals at accepted product | Diagnostic computation only; definitions and thresholds unchanged | Detects false termination consistently | Low | Small overhead, improved reliability |
| Factor/Gram condition telemetry | Read-only diagnostics: singular values, norms, QP iterations and active-set size | Establishes mechanisms without changing results | Low | Small overhead |

## Important qualifications

The current solver already warm-starts each SLSQP call from the current row/time factor. It does not preserve or warm-start the active set. It also already uses `numpy.linalg.lstsq`, so merely replacing normal equations is not a complete solution.

Permanent width reduction based on a numerical singular-value threshold is not semantically neutral: it can remove feasible directions. Dormant-direction management is exact only if it preserves the product at the transformation point, keeps the full cap+1 width available, and uses no singular-value cutoff for statistical rank or acceptance.

## Changes requiring a paper-level decision

The following are outside ordinary engineering and are not proposed for the baseline: changing `B`; changing pilot or reporting caps; changing the number or identity/distribution of maintained starts; loosening stationarity, KKT, or objective-agreement tolerances; accepting one valid start; adding a singular-value rank threshold; changing the ridge or rank-zero anchor; using approximate SVD values in ratios; or replacing the pilot by another estimator. Any such change requires a new paper-level decision.
