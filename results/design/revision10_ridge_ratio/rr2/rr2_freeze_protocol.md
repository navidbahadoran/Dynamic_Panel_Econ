# Proposed Revision-10 theory freeze protocol (not activated)

The following objects would be frozen before implementation or finite-sample evaluation if the cap+1 pilot
gap is resolved by a new paper-level decision. RR2-PARTIAL does **not** activate this protocol.

1. The block collection and reporting caps remain fixed. For the locked three-block Monte Carlo design the
   reporting vector is `(3,3,3)`; other applications use their prespecified matrix-specific fixed caps.
2. The one joint spectral pilot uses caps `bar_r_M+1` and the literal box `||Theta||_max<=B`.
3. Input variables are mapped to stored canonical scientific units before applying the literal box.
4. The final spectrum is
   `lambda_hat_M,j=(w_M/w_A,1)^2 sigma_j(M_hat^pil)^2/(NT)`.
5. The ridge is exactly `a_NT=1/log(NT)`.
6. The zero-rank anchor is exactly `lambda_hat_M,0=1`.
7. Ratios are exactly
   `R_M,j=(lambda_hat_M,j+1+a_NT)/(lambda_hat_M,j+a_NT)`, `j=0,...,bar_r_M`.
8. Each rank is the smallest blockwise argmin. Exact numerical ties go to the smaller rank.
9. Truth at a reporting cap is assessed with the genuine cap+1 pilot singular value; cap+1 is never reported.
10. A missing, infeasible, invalid, or numerically unresolved pilot yields
    `pilot_numerically_unresolved`, no selected rank, and no selected-rank inference.
11. The final coefficient estimator is one literal unpenalized fixed-rank joint post-refit at `r_hat`.
12. The nuclear path is not part of the statistical selector. It may only supply optional deterministic warm
    starts, and its grid cannot define or alter rank.
13. Numerical diagnostics are frozen to feasibility, all start objectives, best stable objective gap,
    objective stability, stationarity/KKT residual, box activity, numerical rank, starting-value identity,
    and runtime. They do not certify the unknown global infimum.
14. Theory requires pilot operator error `O_p(sqrt(NT zeta_NT))`; the cap+1 sufficient optimization
    condition is normalized global objective gap `o_p(zeta_NT)`.
15. Rectangular asymptotics retain the exact `1/N+1/T` rate and impose no balance or relative-growth condition.
16. Revision-9 IC, candidates, threshold, penalty, and local completion have no Revision-10 statistical role.

If subsequently activated, none of these objects may change after seeing RR5 finite-sample results without reopening the paper-level
theory decision. In particular, there is no outcome-selected ridge, anchor, normalization, cap, tie rule, or
failure override.

Before activation, the paper must either extend prediction identification to the cap+1 difference class,
adopt a high-level operator-rate condition for the implemented pilot, or prove a different extra-spectrum
construction under the maintained assumptions. No RR3/RR4 implementation is authorized by this document.
