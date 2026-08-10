# Phase-1.5 theorem skeleton

## Candidate theorem, conditional on unresolved lemmas

Maintain `a:stab`, `a:exog`, `a:geometry`, `a:ned`, `a:moments`, `a:signal`,
`a:identification`, `a:gram`, and `a:growth`; fixed block count and rank caps; the literal box
constraint with interior truth; and `r_0 in R_max`. Add only:

1. **N1 — uniform observable block-envelope validity.** The cap-residual, spatially weighted
   predictable-variance estimators and effective curvatures make the stated Freedman-type
   envelopes simultaneously conservative for all L profile-null gains, with envelopes `o_p(NT)`.
2. **N2 — profiled marginal detectability.** Every successive missing block rank has population
   profiled gain at least `c_gain NT`.

Require all profile objective errors to be negligible relative to the certified envelope/signal
margins and numerical failure probability to vanish.

Then the conditional profiled selector would satisfy

`P_0(rhat_ST=r_0|C)->_p 1`.

## Proof decomposition

- **Lemma A — observable normalizer validity.** Prove N1 from maintained spatial/mixing
  conditions, including cap-residual replacement and a growing matrix predictable-variance
  estimate. This is currently missing.
- **Lemma B — no false increment, including mixed states.** Nuisance-at-cap profiling makes every
  other true block feasible. Lemma A then controls any increment with `s>=r_M,0`.
- **Lemma C — missing direction detected.** N2, uniform concentration, and objective fidelity give
  `Delta_M^prof(s)/Ehat_M(s)->infinity` for all `s<r_M,0`.
- **Lemma D — finite path reaches truth.** Each block's scalar scan passes every missing rank and
  is independent of the scan order; fixed caps make the argument finite.
- **Lemma E — stop at truth.** Lemma B rejects the first unnecessary increment; for a truth-at-cap
  block the scan terminates at the cap.
- **Theorem — joint exact recovery.** Intersect the finitely many lemma events.

On exact recovery, the selected full-panel post-refit has exactly the supplied true ranks. The
existing recovery, target expansion, two-way split correction, Riesz construction, spatial
variance, and selected-rank inference results therefore transfer without substantive changes.

The skeleton stays within two primitive conditions, but N1 is not derived by the existing
manuscript. Treating it as an assumption would hide the main unresolved theorem rather than solve
it.
