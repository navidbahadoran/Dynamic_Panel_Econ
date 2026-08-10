# Conditional selector specification — not authorized

Phase 1.5 does **not** solve the observable envelope, so no unconditional Revision-10 selector is
specified or authorized. If the candidate Freedman/HAC envelope in
`simultaneous_envelope_derivation.md` were proved, the complete mathematical rule would be:

1. Set `L=sum_M cap_M` and use the theorem-prescribed deterministic `alpha_NT -> 0`.
2. Compute the common rank-at-most-cap residual fit and all required block profile endpoints.
3. For every block M and `s=0,...,cap_M-1`, form
   `Delta_M^prof(s)=RSS_M^prof(s)-RSS_M^prof(s+1)`.
4. Form the proved observable `vhat_M`, effective empirical curvature `muhat_M(s)`, increment
   bound `Rhat_M`, and
   `x=log((N+T)L/alpha_NT)`.
5. Set
   `Ehat_M(s)=[sqrt(2 vhat_M x)+(Rhat_M x)/3]^2/muhat_M(s)
               +2 epsilon_RSS`.
6. Select
   `rhat_M=min{s<cap_M: Delta_M^prof(s)<=Ehat_M(s)}`;
   if the set is empty, select `cap_M`.
7. Jointly compute the usual unpenalized exact-rank post-refit at `rhat`; all inference then uses
   that one full-panel fit.

This rule is block-order invariant; deterministic manuscript block order is needed only for
report serialization, not statistics. Rank zero is selected when the s=0 increment fails. Caps
are never expanded. Screening omissions are filled by dynamic profile completion.

Every required endpoint must satisfy multi-start validity, stationarity, the literal coefficient
bound, and absolute RSS suboptimality at most `epsilon_RSS,NT`. A usable truth-free numerical
condition is

`max_endpoint epsilon_RSS,NT = o_p(min{Escale_NT, NT})`,

where `Escale_NT` is the theorem's deterministic lower order for the simultaneous envelope, plus
a no-near-boundary condition that ideal null and signal comparisons are separated from their
boundaries by more than `2 epsilon_RSS,NT` with probability tending to one. Operationally, an
endpoint that does not certify its prespecified tolerance makes the selector unresolved.

The undefined theorem-prescribed `alpha_NT`, unproved `vhat_M` validity, and unproved observable
increment bound are intentionally visible. Consequently this conditional pseudocode is not an
ST1 algorithm and must not be implemented.
