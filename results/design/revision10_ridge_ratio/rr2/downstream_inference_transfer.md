# Downstream inference transfer

Let `E_NT={r_hat=r_0}`. The final selected estimator is defined by the same literal constrained
least-squares map as the supplied-rank estimator:

`Theta_hat^sel=Theta_hat(r_hat)` and `Theta_hat^oracle=Theta_hat(r_0)`.

With the same deterministic start/tie conventions for a numerical representative, on `E_NT` the rank
restrictions, objective, box constraint, and returned estimator coincide. Every deterministic downstream map
of that estimator also coincides: aligned singular spaces, tangent spaces, empirical Riesz equations, target
plug-ins, four supplied-rank split fits, two-way correction, residuals, spatial variance, and studentized
statistics.

For any oracle statistic `Z_NT` and selected statistic `Z_NT^sel` constructed by these same maps,

`P_0(Z_NT^sel!=Z_NT | C)<=P_0(E_NT^c|C)->_p0`.

Therefore any conditional `o_p`, `O_p`, consistency, expansion, or weak-limit statement for the oracle
statistic transfers by asymptotic equivalence. In particular:

- coefficient and singular-space recovery in `thm:recovery` transfers;
- empirical tangent-space and Riesz consistency transfers;
- the target expansion in `thm:target_expansion` transfers;
- the four-fit correction `3 phi_full-phi_time_sum-phi_unit_sum` transfers;
- the distinct full/time/unit residual construction and spatial variance transfer;
- feasible conditional normality in `thm:twoway` transfers.

No substantive re-proof of those results is required. Manuscript changes are limited to replacing the
Revision-9 selector cross-reference, stating the final fixed-rank post-refit convention, and invoking the new
rank-consistency event. Rank selection remains full-panel only and is not repeated in splits.
