# Candidate/profile completion

Requiring the nuclear screen to contain every needed upward neighbor makes consistency depend on
an accidental screening omission. A full Cartesian `R_max` search is unnecessary with fixed caps.

For the original greedy design, the clean rule would be dynamic upward-neighbor completion:
whenever state r is reached, compute every feasible `r+e_M` not already cached. However, the
mixed-state analysis shows that these current-state gains are not protected against omitted
cross-block signal under the maintained assumptions.

For the order-invariant profiled design, the corresponding rule is dynamic profile completion:
for every required pair `(M,s)`, compute the two nested rank-at-most profile fits defining
`Delta_M^prof(s)`, irrespective of whether screening proposed their full rank vectors. Cache
shared endpoints. This needs at most `sum_M(cap_M+1)` conceptual profile endpoints, not the
Cartesian product of ranks.

The nuclear path and thresholded cap pilot remain useful starting values and screening
diagnostics, but neither controls theoretical availability. If any required profile endpoint is
numerically unresolved, the selector returns `numerically_unresolved` rather than treating the
missing comparison as acceptance or rejection. This is a design statement only; no completion
routine is implemented.
