# Mixed-state analysis

Suppose block M is already at `r_M,0` while block L is underfit. The current residual has the
schematic decomposition

`e(r)=u + X_L(M_L,0-M_L(r)) + other estimation terms`.

Therefore

`S_M^raw(r)=Z_M elementwise u
             + Z_M elementwise X_L(omitted L signal)+...`.

The second term need not be centered or small. Prediction identification prevents exact
cancellation between complete coefficient collections; it does not impose cross-block
orthogonality. Projecting onto the current M normal space or residualizing against the current
tangent removes fitted nuisance directions, not an omitted normal direction of L. Consequently a
correctly ranked M can profitably absorb part of L's omitted signal. Largest-only, add-all, and
joint-maximum greedy updates all inherit this contamination. Batch updating changes numerical
order but does not create a no-false-increment proof.

## Profiling remedy

For each block M and scalar s define

`RSS_M^prof(s)=min RSS(Theta)`

subject to `rank(M)<=s`, `rank(L)<=cap_L` for every `L!=M`, and the manuscript's literal box
constraint. Let

`Delta_M^prof(s)=RSS_M^prof(s)-RSS_M^prof(s+1)`.

When all true ranks lie within their caps, every nuisance truth is feasible in both profile
problems. Omitted signal from L is therefore not forced into M. For `s>=r_M,0`, both profile
classes contain the complete truth and the increment is a pure overfit comparison. For
`s<r_M,0`, M remains genuinely restricted regardless of nuisance profiling. This resolves the
mixed-state *definition* and makes block decisions order-invariant.

The remedy is not yet a selector theorem. Current identification and strong singular values give
an order-NT gap between any underfit profile class and truth, but they do not guarantee that every
successive marginal difference `Delta_M^prof(s)` is order NT in a general joint prediction norm.
That is exactly the localized hereditary-detectability condition retained from Phase 1.

Nuisance-at-cap profiling changes the Phase-1 greedy architecture but not the estimator, DGP, box
constraint, or downstream inference. It requires at most `sum_M(cap_M+1)` rank-at-most profile
fits rather than a full Cartesian `R_max` search. No such fit is implemented here.
