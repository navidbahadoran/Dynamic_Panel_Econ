# Final theorem-aligned statistical preflight

This is a 3-replication-per-cell diagnostic preflight, not publication Monte Carlo evidence.
No trimming, winsorization, magnitude filtering, or runtime filtering was applied.
The Parquet replication records are lossless. Their CSV companions omit only the two very large nested rank-selection diagnostic columns, which remain in Parquet and the dedicated rank/fit files.

## Exact commands

```text
.\.venv\Scripts\python.exe scripts\run_mc.py --config configs\mc\preflight_finalized_fixed.toml --dry-run --print-resolved-config
.\.venv\Scripts\python.exe scripts\run_mc.py --config configs\mc\preflight_finalized_selected.toml --dry-run --print-resolved-config
.\.venv\Scripts\python.exe scripts\run_mc.py --config configs\mc\preflight_finalized_fixed.toml
.\.venv\Scripts\python.exe scripts\run_mc.py --config configs\mc\preflight_finalized_selected.toml
```

## Design and matched draws

- Unique semantic DGP realizations: 24 (required 24).
- Method evaluations: 48 (24 fixed rank, 24 selected rank).
- All fixed/selected DGP realization hashes match: PASS.
- Frozen calibration SHA-256: `e8983cadc4fbca990feeba6363420542a99ee056cf445e867532e2a6ea0e7d62`.
- `B=10`, `c_B=1`; every requested frozen cell satisfies `C_Theta <= 9`.

## Fixed-rank reliability

| dgp | N | target_group | valid_full_panel_fits | constrained_fallback_rate | boundary_active_rate | four_required_split_completion_rate | point_retained_share | inference_retained_share | median_runtime_seconds |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| 1 | 50 | broad | 3 | 0.3333 | 0.3333 | 1.0000 | 1.0000 | 1.0000 | 1.5904 |
| 1 | 50 | local_plugin | 3 | 0.3333 | 0.3333 |  | 1.0000 | 1.0000 | 1.5904 |
| 1 | 100 | broad | 3 | 0.0000 | 0.0000 | 1.0000 | 1.0000 | 1.0000 | 5.9860 |
| 1 | 100 | local_plugin | 3 | 0.0000 | 0.0000 |  | 1.0000 | 1.0000 | 5.9860 |
| 2 | 50 | broad | 3 | 0.3333 | 0.3333 | 1.0000 | 1.0000 | 1.0000 | 1.9632 |
| 2 | 50 | local_plugin | 3 | 0.3333 | 0.3333 |  | 1.0000 | 1.0000 | 1.9632 |
| 2 | 100 | broad | 3 | 0.0000 | 0.0000 | 1.0000 | 1.0000 | 1.0000 | 5.4950 |
| 2 | 100 | local_plugin | 3 | 0.0000 | 0.0000 |  | 1.0000 | 1.0000 | 5.4950 |
| 3 | 50 | broad | 3 | 0.6667 | 0.6667 | 1.0000 | 1.0000 | 1.0000 | 8.9312 |
| 3 | 50 | local_plugin | 3 | 0.6667 | 0.6667 |  | 1.0000 | 1.0000 | 8.9312 |
| 3 | 100 | broad | 3 | 0.0000 | 0.0000 | 1.0000 | 1.0000 | 1.0000 | 4.4814 |
| 3 | 100 | local_plugin | 3 | 0.0000 | 0.0000 |  | 1.0000 | 1.0000 | 4.4814 |
| 4 | 50 | broad | 3 | 0.6667 | 0.6667 | 1.0000 | 1.0000 | 1.0000 | 3.7693 |
| 4 | 50 | local_plugin | 3 | 0.6667 | 0.6667 |  | 1.0000 | 1.0000 | 3.7693 |
| 4 | 100 | broad | 3 | 0.0000 | 0.0000 | 1.0000 | 1.0000 | 1.0000 | 5.8141 |
| 4 | 100 | local_plugin | 3 | 0.0000 | 0.0000 |  | 1.0000 | 1.0000 | 5.8141 |

## Selected-rank reliability and complete distributions

| dgp | N | candidate_coverage_rate | exact_rank_recovery_rate | underselection_rate | overselection_rate | rank_cap_hit_rate | selected_rank_distribution_json | point_retained_share | inference_retained_share | median_runtime_seconds |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| 1 | 50 | 1.0000 | 1.0000 | 0.0000 | 0.0000 | 0.0000 | {"[1, 1, 1]": 3} | 1.0000 | 1.0000 | 1175.1245 |
| 1 | 100 | 1.0000 | 0.0000 | 1.0000 | 0.0000 | 0.0000 | {"[1, 0, 0]": 1, "[1, 0, 1]": 2} | 1.0000 | 0.5000 | 871.7581 |
| 2 | 50 | 1.0000 | 0.6667 | 0.0000 | 0.3333 | 0.3333 | {"[1, 1, 1]": 2, "[1, 1, 2]": 1} | 0.6667 | 0.6667 | 1220.7377 |
| 2 | 100 | 1.0000 | 0.3333 | 0.6667 | 0.0000 | 0.0000 | {"[1, 0, 0]": 1, "[1, 0, 1]": 1, "[1, 1, 1]": 1} | 1.0000 | 0.6667 | 1209.5230 |
| 3 | 50 | 1.0000 | 0.6667 | 0.0000 | 0.3333 | 0.3333 | {"[1, 1, 1]": 2, "[1, 1, 2]": 1} | 0.6667 | 0.6667 | 2513.0397 |
| 3 | 100 | 1.0000 | 0.0000 | 1.0000 | 0.0000 | 0.0000 | {"[1, 0, 0]": 1, "[1, 0, 1]": 2} | 1.0000 | 0.5000 | 704.5928 |
| 4 | 50 | 1.0000 | 0.6667 | 0.0000 | 0.3333 | 0.3333 | {"[1, 1, 1]": 2, "[1, 1, 2]": 1} | 0.6667 | 0.6667 | 1111.0716 |
| 4 | 100 | 1.0000 | 0.0000 | 1.0000 | 0.0000 | 0.0000 | {"[1, 0, 0]": 1, "[1, 0, 1]": 2} | 1.0000 | 0.5000 | 1140.1848 |

## Boundary and constrained-estimator behavior

| dgp | N | method | fit_role | attempted_fits | unconstrained_outside_box_rate | constrained_fallback_rate | boundary_active_rate | constrained_solver_failure_rate | constrained_feasibility_failure_rate | constrained_optimality_failure_rate |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| 1 | 50 | fixed_rank | full_panel | 9 | 0.0000 | 0.0000 | 0.0000 | 0.0000 | 0.0000 | 0.0000 |
| 1 | 50 | fixed_rank | time_half_split | 6 | 0.1667 | 0.1667 | 0.1667 | 0.0000 | 0.0000 | 0.0000 |
| 1 | 50 | fixed_rank | unit_half_split | 6 | 0.0000 | 0.0000 | 0.0000 | 0.0000 | 0.0000 | 0.0000 |
| 1 | 50 | selected_rank | candidate_post_refit | 132 | 0.3409 | 0.3409 | 0.3258 | 0.0000 | 0.0000 | 0.0530 |
| 1 | 50 | selected_rank | cap_pilot | 298 | 0.3221 | 0.3221 | 0.3087 | 0.0000 | 0.0000 | 0.0470 |
| 1 | 50 | selected_rank | full_panel | 3 | 0.0000 | 0.0000 | 0.0000 | 0.0000 | 0.0000 | 0.0000 |
| 1 | 50 | selected_rank | time_half_split | 6 | 0.1667 | 0.1667 | 0.1667 | 0.0000 | 0.0000 | 0.0000 |
| 1 | 50 | selected_rank | unit_half_split | 6 | 0.0000 | 0.0000 | 0.0000 | 0.0000 | 0.0000 | 0.0000 |
| 1 | 100 | fixed_rank | full_panel | 9 | 0.0000 | 0.0000 | 0.0000 | 0.0000 | 0.0000 | 0.0000 |
| 1 | 100 | fixed_rank | time_half_split | 6 | 0.0000 | 0.0000 | 0.0000 | 0.0000 | 0.0000 | 0.0000 |
| 1 | 100 | fixed_rank | unit_half_split | 6 | 0.0000 | 0.0000 | 0.0000 | 0.0000 | 0.0000 | 0.0000 |
| 1 | 100 | selected_rank | candidate_post_refit | 142 | 0.1338 | 0.1338 | 0.1197 | 0.0000 | 0.0000 | 0.0070 |
| 1 | 100 | selected_rank | cap_pilot | 398 | 0.1055 | 0.1055 | 0.0804 | 0.0000 | 0.0000 | 0.0025 |
| 1 | 100 | selected_rank | full_panel | 3 | 0.0000 | 0.0000 | 0.0000 | 0.0000 | 0.0000 | 0.0000 |
| 1 | 100 | selected_rank | time_half_split | 6 | 0.0000 | 0.0000 | 0.0000 | 0.0000 | 0.0000 | 0.0000 |
| 1 | 100 | selected_rank | unit_half_split | 6 | 0.0000 | 0.0000 | 0.0000 | 0.0000 | 0.0000 | 0.0000 |
| 2 | 50 | fixed_rank | full_panel | 9 | 0.3333 | 0.3333 | 0.3333 | 0.0000 | 0.0000 | 0.0000 |
| 2 | 50 | fixed_rank | time_half_split | 6 | 0.1667 | 0.1667 | 0.1667 | 0.0000 | 0.0000 | 0.0000 |
| 2 | 50 | fixed_rank | unit_half_split | 6 | 0.1667 | 0.1667 | 0.1667 | 0.0000 | 0.0000 | 0.0000 |
| 2 | 50 | selected_rank | candidate_post_refit | 133 | 0.4286 | 0.4286 | 0.3759 | 0.0000 | 0.0000 | 0.0000 |
| 2 | 50 | selected_rank | cap_pilot | 329 | 0.4255 | 0.4255 | 0.3830 | 0.0000 | 0.0000 | 0.0030 |
| 2 | 50 | selected_rank | full_panel | 3 | 0.3333 | 0.3333 | 0.3333 | 0.0000 | 0.0000 | 0.0000 |
| 2 | 50 | selected_rank | time_half_split | 4 | 0.0000 | 0.0000 | 0.0000 | 0.0000 | 0.0000 | 0.0000 |
| 2 | 50 | selected_rank | unit_half_split | 4 | 0.0000 | 0.0000 | 0.0000 | 0.0000 | 0.0000 | 0.0000 |
| 2 | 100 | fixed_rank | full_panel | 9 | 0.0000 | 0.0000 | 0.0000 | 0.0000 | 0.0000 | 0.0000 |
| 2 | 100 | fixed_rank | time_half_split | 6 | 0.0000 | 0.0000 | 0.0000 | 0.0000 | 0.0000 | 0.0000 |
| 2 | 100 | fixed_rank | unit_half_split | 6 | 0.0000 | 0.0000 | 0.0000 | 0.0000 | 0.0000 | 0.0000 |
| 2 | 100 | selected_rank | candidate_post_refit | 125 | 0.1600 | 0.1600 | 0.1520 | 0.0000 | 0.0000 | 0.0000 |
| 2 | 100 | selected_rank | cap_pilot | 319 | 0.1223 | 0.1223 | 0.1191 | 0.0000 | 0.0000 | 0.0000 |
| 2 | 100 | selected_rank | full_panel | 3 | 0.0000 | 0.0000 | 0.0000 | 0.0000 | 0.0000 | 0.0000 |
| 2 | 100 | selected_rank | time_half_split | 6 | 0.0000 | 0.0000 | 0.0000 | 0.0000 | 0.0000 | 0.0000 |
| 2 | 100 | selected_rank | unit_half_split | 6 | 0.0000 | 0.0000 | 0.0000 | 0.0000 | 0.0000 | 0.0000 |
| 3 | 50 | fixed_rank | full_panel | 9 | 0.3333 | 0.3333 | 0.3333 | 0.0000 | 0.0000 | 0.0000 |
| 3 | 50 | fixed_rank | time_half_split | 6 | 0.1667 | 0.1667 | 0.1667 | 0.0000 | 0.0000 | 0.0000 |
| 3 | 50 | fixed_rank | unit_half_split | 6 | 0.3333 | 0.3333 | 0.3333 | 0.0000 | 0.0000 | 0.0000 |
| 3 | 50 | selected_rank | candidate_post_refit | 137 | 0.4526 | 0.4526 | 0.4453 | 0.0000 | 0.0000 | 0.0219 |
| 3 | 50 | selected_rank | cap_pilot | 403 | 0.4541 | 0.4541 | 0.4392 | 0.0000 | 0.0000 | 0.0174 |
| 3 | 50 | selected_rank | full_panel | 3 | 0.3333 | 0.3333 | 0.3333 | 0.0000 | 0.0000 | 0.0000 |
| 3 | 50 | selected_rank | time_half_split | 4 | 0.0000 | 0.0000 | 0.0000 | 0.0000 | 0.0000 | 0.0000 |
| 3 | 50 | selected_rank | unit_half_split | 4 | 0.2500 | 0.2500 | 0.2500 | 0.0000 | 0.0000 | 0.0000 |
| 3 | 100 | fixed_rank | full_panel | 9 | 0.0000 | 0.0000 | 0.0000 | 0.0000 | 0.0000 | 0.0000 |
| 3 | 100 | fixed_rank | time_half_split | 6 | 0.0000 | 0.0000 | 0.0000 | 0.0000 | 0.0000 | 0.0000 |
| 3 | 100 | fixed_rank | unit_half_split | 6 | 0.0000 | 0.0000 | 0.0000 | 0.0000 | 0.0000 | 0.0000 |
| 3 | 100 | selected_rank | candidate_post_refit | 128 | 0.1875 | 0.1875 | 0.1875 | 0.0000 | 0.0000 | 0.0000 |
| 3 | 100 | selected_rank | cap_pilot | 317 | 0.1451 | 0.1451 | 0.1451 | 0.0000 | 0.0000 | 0.0000 |
| 3 | 100 | selected_rank | full_panel | 3 | 0.0000 | 0.0000 | 0.0000 | 0.0000 | 0.0000 | 0.0000 |
| 3 | 100 | selected_rank | time_half_split | 6 | 0.0000 | 0.0000 | 0.0000 | 0.0000 | 0.0000 | 0.0000 |
| 3 | 100 | selected_rank | unit_half_split | 6 | 0.0000 | 0.0000 | 0.0000 | 0.0000 | 0.0000 | 0.0000 |
| 4 | 50 | fixed_rank | full_panel | 9 | 0.3333 | 0.3333 | 0.3333 | 0.0000 | 0.0000 | 0.0000 |
| 4 | 50 | fixed_rank | time_half_split | 6 | 0.1667 | 0.1667 | 0.1667 | 0.0000 | 0.0000 | 0.0000 |
| 4 | 50 | fixed_rank | unit_half_split | 6 | 0.5000 | 0.5000 | 0.5000 | 0.0000 | 0.0000 | 0.0000 |
| 4 | 50 | selected_rank | candidate_post_refit | 139 | 0.4317 | 0.4317 | 0.4173 | 0.0000 | 0.0000 | 0.0288 |
| 4 | 50 | selected_rank | cap_pilot | 376 | 0.4840 | 0.4840 | 0.4787 | 0.0000 | 0.0000 | 0.0053 |
| 4 | 50 | selected_rank | full_panel | 3 | 0.6667 | 0.6667 | 0.6667 | 0.0000 | 0.0000 | 0.0000 |
| 4 | 50 | selected_rank | time_half_split | 4 | 0.2500 | 0.2500 | 0.2500 | 0.0000 | 0.0000 | 0.0000 |
| 4 | 50 | selected_rank | unit_half_split | 4 | 0.5000 | 0.5000 | 0.5000 | 0.0000 | 0.0000 | 0.0000 |
| 4 | 100 | fixed_rank | full_panel | 9 | 0.0000 | 0.0000 | 0.0000 | 0.0000 | 0.0000 | 0.0000 |
| 4 | 100 | fixed_rank | time_half_split | 6 | 0.0000 | 0.0000 | 0.0000 | 0.0000 | 0.0000 | 0.0000 |
| 4 | 100 | fixed_rank | unit_half_split | 6 | 0.0000 | 0.0000 | 0.0000 | 0.0000 | 0.0000 | 0.0000 |
| 4 | 100 | selected_rank | candidate_post_refit | 137 | 0.2336 | 0.2336 | 0.1971 | 0.0000 | 0.0000 | 0.0000 |
| 4 | 100 | selected_rank | cap_pilot | 326 | 0.1411 | 0.1411 | 0.1350 | 0.0000 | 0.0000 | 0.0000 |
| 4 | 100 | selected_rank | full_panel | 3 | 0.0000 | 0.0000 | 0.0000 | 0.0000 | 0.0000 | 0.0000 |
| 4 | 100 | selected_rank | time_half_split | 6 | 0.0000 | 0.0000 | 0.0000 | 0.0000 | 0.0000 | 0.0000 |
| 4 | 100 | selected_rank | unit_half_split | 6 | 0.0000 | 0.0000 | 0.0000 | 0.0000 | 0.0000 | 0.0000 |

Successful boundary-active constrained estimates are retained in the primary results. The interior-only comparison is secondary and is stored separately.

## Target retention

| dgp | N | method | target_group | point_retained_share | inference_retained_share | interior_only_attempts | interior_only_inference_retained_share |
| --- | --- | --- | --- | --- | --- | --- | --- |
| 1 | 50 | fixed_rank | broad | 1.0000 | 1.0000 | 16 | 1.0000 |
| 1 | 50 | fixed_rank | local_plugin | 1.0000 | 1.0000 | 20 | 1.0000 |
| 1 | 50 | selected_rank | broad | 1.0000 | 1.0000 | 0 |  |
| 1 | 50 | selected_rank | local_plugin | 1.0000 | 1.0000 | 0 |  |
| 1 | 100 | fixed_rank | broad | 1.0000 | 1.0000 | 24 | 1.0000 |
| 1 | 100 | fixed_rank | local_plugin | 1.0000 | 1.0000 | 30 | 1.0000 |
| 1 | 100 | selected_rank | broad | 1.0000 | 0.5000 | 0 |  |
| 1 | 100 | selected_rank | local_plugin | 1.0000 | 0.5000 | 0 |  |
| 2 | 50 | fixed_rank | broad | 1.0000 | 1.0000 | 16 | 1.0000 |
| 2 | 50 | fixed_rank | local_plugin | 1.0000 | 1.0000 | 20 | 1.0000 |
| 2 | 50 | selected_rank | broad | 1.0000 | 1.0000 | 0 |  |
| 2 | 50 | selected_rank | local_plugin | 0.5263 | 0.5263 | 0 |  |
| 2 | 100 | fixed_rank | broad | 1.0000 | 1.0000 | 24 | 1.0000 |
| 2 | 100 | fixed_rank | local_plugin | 1.0000 | 1.0000 | 30 | 1.0000 |
| 2 | 100 | selected_rank | broad | 1.0000 | 0.6667 | 0 |  |
| 2 | 100 | selected_rank | local_plugin | 1.0000 | 0.6667 | 0 |  |
| 3 | 50 | fixed_rank | broad | 1.0000 | 1.0000 | 8 | 1.0000 |
| 3 | 50 | fixed_rank | local_plugin | 1.0000 | 1.0000 | 10 | 1.0000 |
| 3 | 50 | selected_rank | broad | 1.0000 | 1.0000 | 0 |  |
| 3 | 50 | selected_rank | local_plugin | 0.5263 | 0.5263 | 0 |  |
| 3 | 100 | fixed_rank | broad | 1.0000 | 1.0000 | 24 | 1.0000 |
| 3 | 100 | fixed_rank | local_plugin | 1.0000 | 1.0000 | 30 | 1.0000 |
| 3 | 100 | selected_rank | broad | 1.0000 | 0.5000 | 0 |  |
| 3 | 100 | selected_rank | local_plugin | 1.0000 | 0.5000 | 0 |  |
| 4 | 50 | fixed_rank | broad | 1.0000 | 1.0000 | 8 | 1.0000 |
| 4 | 50 | fixed_rank | local_plugin | 1.0000 | 1.0000 | 10 | 1.0000 |
| 4 | 50 | selected_rank | broad | 1.0000 | 1.0000 | 0 |  |
| 4 | 50 | selected_rank | local_plugin | 0.5263 | 0.5263 | 0 |  |
| 4 | 100 | fixed_rank | broad | 1.0000 | 1.0000 | 24 | 1.0000 |
| 4 | 100 | fixed_rank | local_plugin | 1.0000 | 1.0000 | 30 | 1.0000 |
| 4 | 100 | selected_rank | broad | 1.0000 | 0.5000 | 0 |  |
| 4 | 100 | selected_rank | local_plugin | 1.0000 | 0.5000 | 0 |  |

## Gram/Riesz diagnostics for theorem-covered targets

| dgp | N | method | target_group | target_support_rate | tangent_Gram_failure_rate | Riesz_failure_rate_conditional_on_support | invalid_variance_rate | minimum_Gram_eigenvalue | p90_condition_number |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| 1 | 50 | fixed_rank | broad | 1.0000 | 0.0000 | 0.0000 | 0.0000 | 0.1406 | 328.4175 |
| 1 | 50 | fixed_rank | local_plugin | 1.0000 | 0.0000 | 0.0000 | 0.0000 | 0.1406 | 328.4175 |
| 1 | 50 | selected_rank | broad | 1.0000 | 0.0000 | 0.0000 | 0.0000 | 0.1407 | 328.3216 |
| 1 | 50 | selected_rank | local_plugin | 1.0000 | 0.0000 | 0.0000 | 0.0000 | 0.1407 | 328.3216 |
| 1 | 100 | fixed_rank | broad | 1.0000 | 0.0000 | 0.0000 | 0.0000 | 0.1765 | 223.2192 |
| 1 | 100 | fixed_rank | local_plugin | 1.0000 | 0.0000 | 0.0000 | 0.0000 | 0.1765 | 223.2192 |
| 1 | 100 | selected_rank | broad | 0.5000 | 0.0000 | 0.0000 | 0.5000 | 0.1533 | 260.6199 |
| 1 | 100 | selected_rank | local_plugin | 0.5714 | 0.0000 | 0.0000 | 0.4286 | 0.1533 | 260.6199 |
| 2 | 50 | fixed_rank | broad | 1.0000 | 0.0000 | 0.0000 | 0.0000 | 0.0681 | 511.3852 |
| 2 | 50 | fixed_rank | local_plugin | 1.0000 | 0.0000 | 0.0000 | 0.0000 | 0.0681 | 511.3852 |
| 2 | 50 | selected_rank | broad | 1.0000 | 0.0000 | 0.0000 | 0.0000 | 0.0759 | 520.7137 |
| 2 | 50 | selected_rank | local_plugin | 1.0000 | 0.0000 | 0.0000 | 0.0000 | 0.0759 | 520.7137 |
| 2 | 100 | fixed_rank | broad | 1.0000 | 0.0000 | 0.0000 | 0.0000 | 0.2326 | 193.0500 |
| 2 | 100 | fixed_rank | local_plugin | 1.0000 | 0.0000 | 0.0000 | 0.0000 | 0.2326 | 193.0500 |
| 2 | 100 | selected_rank | broad | 0.6667 | 0.0000 | 0.0000 | 0.3333 | 0.2568 | 191.8123 |
| 2 | 100 | selected_rank | local_plugin | 0.7143 | 0.0000 | 0.0000 | 0.2857 | 0.2568 | 191.8123 |
| 3 | 50 | fixed_rank | broad | 1.0000 | 0.0000 | 0.0000 | 0.0000 | 0.1546 | 388.3915 |
| 3 | 50 | fixed_rank | local_plugin | 1.0000 | 0.0000 | 0.0000 | 0.0000 | 0.1546 | 388.3915 |
| 3 | 50 | selected_rank | broad | 1.0000 | 0.0000 | 0.0000 | 0.0000 | 0.2178 | 253.4542 |
| 3 | 50 | selected_rank | local_plugin | 1.0000 | 0.0000 | 0.0000 | 0.0000 | 0.2178 | 253.4542 |
| 3 | 100 | fixed_rank | broad | 1.0000 | 0.0000 | 0.0000 | 0.0000 | 0.1838 | 282.8910 |
| 3 | 100 | fixed_rank | local_plugin | 1.0000 | 0.0000 | 0.0000 | 0.0000 | 0.1838 | 282.8910 |
| 3 | 100 | selected_rank | broad | 0.5000 | 0.0000 | 0.0000 | 0.5000 | 0.2627 | 195.1272 |
| 3 | 100 | selected_rank | local_plugin | 0.5714 | 0.0000 | 0.0000 | 0.4286 | 0.2627 | 195.1272 |
| 4 | 50 | fixed_rank | broad | 1.0000 | 0.0000 | 0.0000 | 0.0000 | 0.1548 | 425.6356 |
| 4 | 50 | fixed_rank | local_plugin | 1.0000 | 0.0000 | 0.0000 | 0.0000 | 0.1548 | 425.6356 |
| 4 | 50 | selected_rank | broad | 1.0000 | 0.0000 | 0.0000 | 0.0000 | 0.1548 | 209.8419 |
| 4 | 50 | selected_rank | local_plugin | 1.0000 | 0.0000 | 0.0000 | 0.0000 | 0.1548 | 209.8419 |
| 4 | 100 | fixed_rank | broad | 1.0000 | 0.0000 | 0.0000 | 0.0000 | 0.1550 | 215.4044 |
| 4 | 100 | fixed_rank | local_plugin | 1.0000 | 0.0000 | 0.0000 | 0.0000 | 0.1550 | 215.4044 |
| 4 | 100 | selected_rank | broad | 0.5000 | 0.0000 | 0.0000 | 0.5000 | 0.2691 | 191.5438 |
| 4 | 100 | selected_rank | local_plugin | 0.5556 | 0.0000 | 0.0000 | 0.4444 | 0.2691 | 191.5438 |

## Runtime and optimization

Lossless fit-role diagnostics are in `optimization_summary.csv`; slow and extreme finite fits remain included. Selected final full-panel rows are a derived role view of the chosen candidate post-refit and are therefore intentionally non-additive with candidate rows.

## Accounting

- Target-cell accounting rows: 288; all have `R_attempted=3`.
- `R_inference <= R_point <= R_attempted`: PASS.
- `inference_valid iff primary_status=success`: PASS.
- Every method-replication contains all 18 requested target records: PASS.
- Exactly four shared split fits per broad-target replication: PASS.

## Implementation finding

The numerical estimator did not require correction. The audit identified a fit-diagnostic labeling defect: additional deterministic fixed-rank full-panel starts were emitted as generic `coefficient_fit`. The combined output relabels those records losslessly, and future instrumentation now labels all pre-split fixed starts `full_fixed_rank`. Explicit aliases `max_abs_coefficient` and `constrained_runtime_seconds` were also added.

## Medium recommendation

- Provisional `c_kappa=4e-6` remains plausible: NO.
- Recommendation for the 100-replication medium diagnostic: **NO-GO**.
- Criterion results: `{"accounting_reconciles": true, "candidate_coverage_high": true, "constrained_failures_rare": true, "fixed_rank_reliable": true, "gram_riesz_failures_rare": true, "rank_cap_hits_not_systematic": true, "selected_behavior_not_degenerate": false, "split_broad_operational": true}`.
- By-size selected-rank behavior: `{"100": {"exact_rank_recovery_rate": 0.08333333333333333, "overselection_rate": 0.0, "underselection_rate": 0.9166666666666666}, "50": {"exact_rank_recovery_rate": 0.75, "overselection_rate": 0.25, "underselection_rate": 0.0}}`.

The medium diagnostic was not launched.
