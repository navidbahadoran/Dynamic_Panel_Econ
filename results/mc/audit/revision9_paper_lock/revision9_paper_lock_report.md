# Revision-9 paper-lock report

Audit basis: repository `main` at `25bea9c1d04ce2a72903c50ce2616ddaab1980b1`. This task changed only files in this audit namespace. It did not alter source, estimator, rank selection, active configuration, README, frozen calibration, or manuscript, and it ran no reported/statistical Monte Carlo.

## Paper-lock result

The Revision-9 DGP candidate is calibration-feasible and verifies the proposed coefficient box.

- `pi_H=0.30` and target pooled R2 `0.65` have positive `c_xi` roots in every baseline and identified-scale rank-stress cell.
- `(1,0,2)` retains `c_xi=1`, `r2_scale_identified=false`, and reports induced R2.
- `C_Theta,max = 8.288745227963506`.
- `9-C_Theta,max = 0.711254772036494`: additional slack relative to the required `B-c_B` envelope.
- `10-C_Theta,max = 1.711254772036494`: distance from the worst-case truth to the estimation boundary.
- Therefore `B=10`, `c_B=1` is **VERIFIED**.

The worst cell is DGP4, rank-stress `(1,1,1)`, `N=T=50`.

## Calibration candidate

The candidate contains 52 cells: 16 baseline cells and 36 stress cells. It uses the current analytical population `c_H`, master seed `20260807`, 50 independent deterministic calibration draws per cell, and exact quadratic pooled-R2 evaluation. The parameter `kappa_f_b=0.15` was supplied in memory only. The active frozen calibration was not overwritten.

For rank-one H, `c_H=0.6546536707079772`. For rank-two rescaled H, `c_H=0.7301712917987002`. Identified-scale achieved R2 differs from 0.65 only by root-solver rounding. Induced R2 for `(1,0,2)` ranges from `0.5210906542635242` to `0.5393014964428011`.

The complete cell table is `final_calibration_candidate.csv`; envelope maxima and precisely named slack measures are in `final_coefficient_envelope.csv`.

## Rank-stress formulas and caps

The exact equations are in `rank_stress_exact_specification.md`, with a ready-to-paste manuscript block in `rank_stress_latex_block.tex`. In brief, `(2,1,1)` adds a unit-uniform times iid-time-uniform outer product to A raw, commonly rescales the two-term A matrix to the rank-one raw support, and then applies the realized full-horizon A stability scale. `(1,0,2)` sets B identically zero and applies the analogous two-term construction/common support rescaling to H, followed by actual-rank population `c_H` and the `c_xi=1` normalization.

The intended production and separate rank-stress configurations both use primary caps `(3,3,3)`, and every true rank lies within them. Historical preflight files use `(2,2,2)`; the paper should identify `(3,3,3)` and the maintained configuration explicitly so those artifacts cannot be mistaken for the final design.

## B-entry, initialization, and conditioning

Revision 9 gives `g_b in [-3,3]` and `f_b in [0.15,1.05]`. Hence the normalized rank-one B time singular vector satisfies the uniform leverage bound `T v_t^2 >= 1/49`. Together with positive bounded B loadings, B-entry is now compatible with target regularity.

All six requested initial conditions are zero in code: `y_i,-50`, `x_i,-50`, `g_a,-50`, `g_b,-50`, `g_h,-50`, and `f_x,-50`. **MATCH.**

The main-DGP conditioning-field mapping is exact: condition on common bounded paths and time-invariant unit design draws, while Gaussian time-varying innovations remain random. **MATCH.** For rank stress, the paper must explicitly add the new bounded factor paths and unit loadings to the corresponding two parts of the field.

## Remaining paper reproducibility items

1. State that the authoritative primary cap is `(3,3,3)` and distinguish historical `(2,2,2)` preflights and the `(4,4,4)` sensitivity.
2. Insert the exact rank-stress subsection, including iid (not AR) added factors, component strength one, common support rescaling, full-horizon A stability, almost-sure rank language, separate seed-labelled calibration, and the zero-B normalization.
3. Extend `C^MC` for rank stress with its added common factor paths and unit loadings.
4. State master seed `20260807`, 50 calibration draws, semantic separation of baseline and rank-stress calibration, and the prohibition on reusing calibration draws as reported replications.
5. State the exact finite-horizon population variance used for rank-two H and that `(1,0,2)` reports induced rather than intended R2.
6. State that `B=10`, `c_B=1`, and the verified worst envelope is `8.288745227963506` for this candidate.
7. Preserve the exact Revision-8 IC baseline; do not describe `4e-6` as the paper penalty. List only 0.5 and 2.0 as sensitivity multipliers.

This is a non-activated paper-lock candidate. Code/config alignment, if later authorized, must be a separate task after the manuscript formulas are locked.

## Repository-health validation

- `pytest -q`: the configured `results/pytest-tmp` could not be removed by the sandbox identity (`WinError 5`), after 118 tests had passed and four setups errored. The unchanged suite was rerun with a fresh task-specific `--basetemp`; **122 passed in 19.43s**.
- `ruff check src scripts tests`: **All checks passed.**
- `git diff --check`: **passed with no output.**

The base-temp override addressed only an environment ownership problem and did not alter test or source files.
