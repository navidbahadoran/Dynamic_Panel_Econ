# Final independent preflight

This run contains 12 new DGP draws (replication indices 3-5) and 24 matched method evaluations.
No medium, rank-stress, power, or production simulation was run.

## Fixed rank

```json
{
  "attempted_replications": 12,
  "full_fit_success": 11,
  "coefficient_bound_hits": 1,
  "numerical_stability_failures": 0,
  "point_retained": 198,
  "point_attempted": 216,
  "point_retained_share": 0.9166666666666666,
  "inference_retained": 134,
  "inference_attempted": 216,
  "inference_retained_share": 0.6203703703703703,
  "gram_failure_target_records": 0,
  "riesz_failure_target_records": 0,
  "split_fit_failure_target_records": 64,
  "split_fit_failure_replications": 8,
  "runtime_total_seconds": 21.37987110001268,
  "runtime_mean_seconds": 1.7816559250010566,
  "runtime_median_seconds": 1.8800524499965832
}
```

## Selected rank (`c_kappa=4e-6`)

```json
{
  "attempted_replications": 12,
  "valid_cap_pilot": 12,
  "pilot_multistart_disagreement": 1,
  "basin_confirmation_attempted": 2,
  "basin_confirmation_success": 1,
  "candidate_coverage": 12,
  "candidate_coverage_conditional_share": 1.0,
  "exact_recovery": 9,
  "underselection": 0,
  "overselection": 3,
  "rank_cap_hits": 3,
  "point_retained": 162,
  "point_attempted": 216,
  "point_retained_share": 0.75,
  "inference_retained": 114,
  "inference_attempted": 216,
  "inference_retained_share": 0.5277777777777778,
  "runtime_total_seconds": 1166.4954091999098,
  "runtime_mean_seconds": 97.20795076665915,
  "runtime_median_seconds": 97.53823699998611,
  "gram_failure_target_records": 0,
  "riesz_failure_target_records": 0,
  "split_fit_failure_target_records": 48,
  "candidate_numerically_unresolved_replications": 0,
  "selected_postrefit_stability_failures": 0,
  "rank_selection_failure_replications": 0,
  "split_fit_failure_replications": 6
}
```

Selected rank distribution:

```csv
selected_rank_vector,count,share
"[1, 1, 1]",9,0.75
"[1, 1, 2]",3,0.25
```

## Recommendation

**NO-GO** for the medium diagnostic. Although all cap pilots were valid and candidate coverage was 100%, rank-cap hits occurred in 3/12 cases, including 2/3 DGP-4 cases. That DGP-specific concentration does not satisfy the no-systematic-rank-cap gate.

The recommendation also notes the explicitly retained split-fit failures; Gram and Riesz primary failures are counted separately above.
The optimization CSV counts every executed route/candidate start, including intentionally rejected starts; replication-level numerical failures are the gate metrics reported in the JSON summaries.
