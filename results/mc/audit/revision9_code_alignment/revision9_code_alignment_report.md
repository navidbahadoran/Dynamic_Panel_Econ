# Locked Revision-9 code-alignment report

## Scope

This alignment used the locked Revision-9 TeX as the sole scientific source. It changed the DGP/defaults, maintained configurations, official frozen calibration, target metadata, boundary-inference status, manifests, README, deterministic calibration builder, and tests. It did not change the estimator objective, rank-selection objective, candidate algorithm, split correction, variance estimator, rank-stress formulas, or outlier policy. No statistical Monte Carlo was run.

## Mismatches found and resolved

Before alignment, `kappa_f_b` was 0.20; `rho_fx` resolved to 0.0 while the generator hard-coded 0.5; B-entry was non-headline; boundary activity did not suppress primary inference; sensitivities included redundant 1.0; the frozen table used the old slope factor; smoke/README examples carried non-paper numerical choices; and manifests omitted the conditioning-field classification. All are aligned as recorded in `revision9_code_compliance.csv`.

## Locked DGP and frozen calibration

The active constants are: R2 target 0.65, `pi_H=.30`, `rho_g=rho_x=rho_fx=rho_s=.5`, `delta_x=.5`, `eta_x=.3` in DGP3-4, `mu_f_a=.5`, `kappa_f_a=.10`, `mu_f_b=.6`, `kappa_f_b=.15`, DGP4 A means `.9/1.1` with scale `.08`, and B means `.8/1.2` with scale `.25`. All six initial conditions at `-50` equal zero.

The dedicated builder regenerated 52 active frozen cells using master seed `20260807` and 50 calibration-only draws per cell. Rank-one `c_H=0.6546536707079772`; rank-two-H `c_H=0.7301712917987002`. Identified `c_xi` ranges from `0.8085955961638839` to `2.4366620615516461`; achieved R2 equals 0.65 within `1e-10`. For `(1,0,2)`, `c_xi=1`, scale identification is false, and induced R2 ranges from `0.5210906542635242` to `0.5393014964428011`.

`C_Theta,max=8.288745227963506`. Hence `8.288745227963506 <= 9`, additional slack to the required envelope is `0.711254772036494`, and distance to the estimation boundary is `1.711254772036494`. `B=10,c_B=1` is verified.

## Rank stress

The current code already matched the locked Appendix. `(2,1,1)` applies `s_A=E_Araw/(E_Araw+3)` to the sum `lambda_a f_a' + tilde_lambda_A tilde_f_A'`, with iid standardized-uniform added loading and iid-over-time standardized-uniform factor, then applies realized full-horizon `c_a`. `(1,0,2)` keeps baseline A, sets B identically zero, and applies `s_H=(3sqrt(3))/(3sqrt(3)+3)` to `lambda_h g_h' + tilde_lambda_H tilde_g_H'`, with the locked design-specific population `c_H` and `c_xi=1` rule.

## Target and inference alignment

The analytic support is `g_b in [-3,3]`, so `f_b in [.15,1.05]` and B-entry is theorem-covered. DGP1-3 fixed-time group contrasts remain `weak_target_stress_outside_assumption9`; DGP4 fixed-time contrasts and all listed time-averaged contrasts are headline targets.

Boundary-active constrained points remain valid for bias/RMSE/MC SD, while primary inference is suppressed with `boundary_interiority_failure`. Selected-rank cap hits continue to suppress inference before target computation.

## Canonical dry run

```powershell
.\.venv\Scripts\python.exe scripts\run_mc.py --config configs\mc\production.toml --pooled-r2-target 0.65 --kappa-f-b 0.15 --coefficient-bound 10 --rank-caps 3,3,3 --ic-multiplier 1.0 --print-resolved-config --dry-run
```

The resolved output displays R2 `.65`, `kappa_f_b=.15`, `rho_fx=.5`, B `10`, caps `[3,3,3]`, IC multiplier `1`, IC/threshold sensitivities `[.5,2]`, gamma `.8`, epsilon `.01`, and the active frozen table. Dry-run exits before calibration, fitting, inference, or output.

## Deterministic validation

No statistical Monte Carlo was launched. Calibration-only regeneration and deterministic validation completed successfully:

- `pytest -q`: 142 passed (the command used an alternate `--basetemp` because the repository-configured Windows temporary directory was not writable in the validation sandbox);
- `ruff check src scripts tests`: all checks passed;
- `git diff --check`: passed with no whitespace errors;
- canonical Revision-9 CLI dry run: resolved the locked values and exited before any simulation work.
