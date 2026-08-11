# Revision-10 cap+1 pilot acceptance gate

## Exact Boolean logic

Source: `src/dynamic_panel_econ/rank_selection.py`, `fit_invalid_reasons` (lines 159-188) and `fit_revision10_spectral_pilot` (lines 1384-1499).

For each of three deterministic seeded starts `s`, with requested factor widths `(4,4,4)`, the code uses `require_exact_numerical_rank=False` and defines:

```text
finite_objective_s = isfinite(objective_s)
termination_s = fit.converged

if constrained_fallback_s:
    stationarity_s = stored constrained stationarity_pass
                     (factor-space box KKT <= 1e-4)
    subproblem_feasibility_s = max_constraint_violation <= 1e-8
else:
    stationarity_s = finite(projected_gradient_residual)
                     and projected_gradient_residual <= 1e-6
    subproblem_feasibility_s = True

literal_box_s = finite(max_envelope_ratio)
                and max_envelope_ratio <= 1 + 1e-8
rank_condition_s = True  # exact numerical rank explicitly disabled for pilot

valid_s = finite_objective_s and termination_s and stationarity_s
          and subproblem_feasibility_s and literal_box_s and rank_condition_s

valid = starts with valid_s, sorted by objective
two_valid = len(valid) >= 2
gap = abs(L_valid_2 - L_valid_1) / max(1, abs(L_valid_1))
stable = two_valid and gap <= 1e-6
pilot_accepted = stable
```

If `stable` is false, `RankPilotFailure` is raised before scale weights, block spectra, ridge ratios, block ranks, or the final post-refit are computed.

## Conditions

| Condition | Source | Frozen value | Role | RR5 status |
|---|---|---:|---|---|
| Three starts | `fit_revision10_spectral_pilot`, lines 1400-1405 | 3 | Purely numerical credibility | Frozen before RR5 |
| Finite objective | `fit_invalid_reasons`, line 186 | finite | Numerical | 540/540 passed |
| Solver termination | lines 167-168 | `fit.converged` | Numerical | 53/540 passed |
| Interior stationarity | `estimation.py`, lines 140-146 and 228-243 | `1e-6` | Numerical | 14/339 passed |
| Constrained KKT | `estimation.py`, lines 520-550 | `1e-4` | Numerical | 39/201 passed |
| Subproblem feasibility | `fit_invalid_reasons`, lines 173-177 | `1e-8` | Numerical enforcement of scientific box | 540/540 passed |
| Literal coefficient box | lines 181-182 | ratio `<= 1+1e-8`, with `B=10` | Scientific constraint plus numerical slack | 540/540 passed |
| Exact numerical rank | line 184; pilot call line 1419 | disabled | Correct rank-at-most semantics | Not a rejection source |
| Boundary activity | recorded only | no rejection rule | Diagnostic | Not a rejection source |
| Best-two objective agreement | line 1458 | normalized gap `<=1e-6` | Purely numerical credibility | 0/180 passed |

There is no IC, threshold rank, candidate enumeration, Revision-9 status check, or final-post-refit condition in this pilot gate.
