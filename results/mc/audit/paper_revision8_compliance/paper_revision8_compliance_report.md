# Revision-8 paper-versus-code compliance report

Audit basis: supplied paper specification, repository `main` at `25bea9c1d04ce2a72903c50ce2616ddaab1980b1`. No source/configuration/paper files were changed and no statistical Monte Carlo was run. The only numerical work was deterministic calibration feasibility using the finalized 50-draw common-random-number construction with master seed `20260807`.

## Decision

**NO-GO for a paper-labelled production Monte Carlo.** The paper's pooled-R2 target 0.50 is infeasible under its own finalized population-`c_H` construction in every requested baseline cell. In addition, current boundary-active inference does not follow the paper's suppression rule. Paper revisions must precede implementation changes.

## A. Paper/code mismatches

1. **Calibration target:** paper 0.50; current frozen/production implementation 0.65. The requested 0.50 is not merely absent: it has no admissible positive `c_xi` root in any of 16 cells.
2. **Boundary inference:** the constrained point estimate is correct and retained, but a boundary-active fit can still receive a primary interval. The paper says to suppress it.
3. **Penalty sensitivity reporting:** the baseline formula and multiplier 1 are correct, but code lists `[0.5,1.0,2.0]` as sensitivity multipliers; the paper says the sensitivities are only 0.5 and 2.0. The extra 1.0 duplicates baseline.
4. **Threshold sensitivity reporting:** the same redundant 1.0 appears in the threshold sensitivity list, while the paper names only 0.5 and 2.0 as sensitivities.

Historical preflight configs containing `4e-6` are not the production baseline. The production config inherits `ic_multiplier=1.0`, so the baseline Revision-8 rate itself matches.

## B. Internal paper inconsistencies and underspecification

### Internal inconsistency

The simultaneous claims `pi_H=0.30`, finalized population `c_H`, and average pooled R2 0.50 cannot all be satisfied by positive `c_xi` for the stated DGP/cells. The paper must change the R2 target or change a substantive DGP/calibration definition. Code must not silently choose one.

### Missing reproducibility details

- Initial values for `g_h`, `f_x`, and `x`, plus the extra primitive disturbance used by DGP3-4 at the start of burn-in.
- The enlarged conditioning field containing time-invariant unit draws.
- Numerical `B`, simulation margin, rank caps, and deterministic envelope formulas.
- Nuclear grid/stopping/weights and the fact that the path is screening only.
- Exact cap-pilot routes, rank moves, tolerances, stopping, and whether failed pilot multistart agreement is fatal.
- Exact candidate union, local-completion loop, multistart validity, and IC tie-breaking.
- Split seeds/partitions, four-fit rule, preservation of external lags, computational split-rank floor, and variance normalization/HAC definition.
- Exact target indices and weights, plus which targets are theorem validations versus stress targets.
- Rank-stress outer-product laws, iid-over-time factors, strengths, independence, common rescaling, stability, calibration seeds/draws/root rule, and zero-B normalization.
- DGP4 group assignment and explicit standardized-uniform perturbation law.

## C. R2=0.50 feasibility

The exact calibration recursion writes each draw as `y(c)=y_base+c*y_scale`, so its pooled R2 is evaluated from stored quadratic forms without simulation replications. The average was minimized over positive `c`; finite stationary points and the `c->infinity` limit were checked.

All 16 lower extrema exceed 0.50. They range from **0.5469907224** (DGP4, 100) to **0.5798064805** (DGP1, 100). The upper supremum is 1 as `c_xi -> 0`. Full cell results are in `r2_050_feasibility.csv`.

There is therefore:

- no `c_xi` solving 0.50;
- no corresponding `C_H` or `C_Theta` truth to report;
- no finite common `(B,c_B)` capable of implementing the nonexistent target calibration.

Moreover, trying to approach the asymptotic lower bounds by sending `c_xi` to infinity sends `C_H=|c_xi*c_H|*3sqrt(3)` to infinity. In the few cells with a shallow finite numerical minimum, the minimizing scale is hundreds or thousands, producing enormous H envelopes. This cannot be reconciled with a fixed finite interior coefficient box.

For reference only, not as a 0.50 solution: `C_A=0.85`; `C_B=2.0313843876` for DGP1-3 and `1.9596152423` for DGP4; and `C_H/c_xi=3.4016802571` at the rank-one population `c_H`. The existing 0.65 frozen truths fit inside production `B=10` with margin 1, but that does not validate the paper's 0.50 statement.

## D. Correct conditioning sigma-field

Use

\[
\mathcal C_{NT}=\sigma(\text{full burn+sample paths }g_a,g_b,g_h,f_x;
 \lambda_a,\lambda_b,\lambda_h,\lambda_x,\sigma_i^2,\sigma_{e,i}^2,G_i;
 c_a,c_H,c_\xi).
\]

For DGP4, `lambda_a` and `lambda_b` include the bounded group-specific perturbations. Do not include time-varying Gaussian disturbance or covariate innovations. The detailed rationale is in `conditioning_sigma_field_audit.md`.

## E-F. Target coverage and B-entry

THEOREM_COVERED: A entry; overall fixed-time A/B means; DGP1-3 fixed-time group means; DGP4 fixed-time group means and contrasts; full-panel A/B means; group-specific time averages; time-averaged group contrasts.

NOT_UNIFORMLY_COVERED: B entry; DGP1-3 fixed-time group contrasts.

No listed target was assigned NEEDS_ADDITIONAL_ASSUMPTION under the current analytic reading; instead, the two nonuniform classes should be explicitly presented as stress targets outside the theorem. If the author wants B-entry inside the theorem, the DGP needs an additional strict factor-floor restriction.

At `kappa_f_b=0.20`, `g_b` has support `[-3,3]` and `f_b` has support `[0,1.2]`. Hence **B-entry is not theorem-covered**. Empirical projection ratios cannot cure this lack of a uniform support floor.

## G. Rank-stress formulas

- `(2,1,1)`: add one independent `U[-sqrt(3),sqrt(3)]` loading times an iid-over-time factor of the same law to A raw; multiply the whole expanded A raw matrix by `E_Araw/(E_Araw+3)`; then apply realized full-horizon `c_a` stability scaling. B and H remain rank one.
- `(1,0,2)`: A remains rank one; B is exactly zero; add the same type of independent outer product to H raw and multiply the whole expanded H raw matrix by `E_Hraw/(E_Hraw+3)`. Recalculate population `c_H` for actual rescaled rank-two H. Set `c_xi=1`, declare `r2_scale_identified=false`, and report induced R2 rather than a failed target.

All cells/rank vectors receive separate actual-matrix calibration. See `rank_stress_specification.md` for the complete formula and reproducibility list.

## H-I. Boundary and cap compliance

- **Boundary:** point estimator MATCH; point retention MATCH; primary-inference suppression **MISMATCH**.
- **Selected rank at cap:** **MATCH**. The code returns `rank_at_cap` before target inference, so no primary interval is produced.

## J. Recommended paper revisions first

1. Replace the infeasible 0.50 claim with a target justified by the feasibility table, or explicitly revise the substantive calibration/DGP. The existing 0.65 choice is feasible but should be adopted only by an author decision.
2. Enlarge `C_NT` exactly as stated above.
3. State all missing initial conditions and disturbance indexing.
4. Classify B-entry and DGP1-3 fixed-time contrasts as outside uniform theorem coverage, or add assumptions/DGP restrictions that genuinely cover them.
5. Add exact rank-stress formulas and the `(1,0,2)` normalization rule.
6. State numerical B/margin/caps and all rank-selection numerical details, especially cap-pilot stability eligibility.
7. State whether boundary suppression applies to the selected full fit only or to any full/split coefficient fit.
8. State exact split and variance formulas and exact sensitivity grids.

## K. Recommended code changes only after paper lock

1. Change/regenerate calibration tables and `target_r2` only after the author selects a feasible paper target.
2. Add a boundary-active inference status that retains the point estimate but suppresses primary intervals; implement the paper's chosen full/split precedence.
3. Remove redundant 1.0 entries from sensitivity *reporting* if the paper literally restricts sensitivities to 0.5 and 2.0; keep 1.0 as baseline.
4. If the paper makes cap-pilot multistart agreement a hard condition, make failed route/confirmation agreement fatal rather than a warning.
5. Align target metadata with the paper's final theorem/stress labels.

Do not start production until the paper is internally feasible, the approved changes are implemented, and the resulting deterministic calibration/envelope checks pass.
