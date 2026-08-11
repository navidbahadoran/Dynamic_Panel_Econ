# Revision-10 result schema

Primary rank records expose stable machine-readable fields:

- `rank_selector_method = revision10_ridge_ratio`;
- `reporting_rank_caps`, `pilot_rank_caps`;
- `scale_weights`, `reference_weight_w_A1`, `a_NT`;
- `pilot_singular_values`, `normalized_lambda_hat`, `ridge_ratios`;
- `selected_rank_by_block`, `selected_rank_vector`;
- `minimum_ratio_by_block`, `second_smallest_ratio_by_block`, `ratio_gap_by_block`;
- `pilot_objective`, `pilot_feasibility`, `pilot_stationarity_residual`,
  `pilot_boundary_activity`;
- `pilot_start_objectives`, `pilot_start_stationarity_residuals`,
  `pilot_best_two_objective_gap`, `pilot_objective_stability_pass`;
- `rank_selection_numerically_unresolved`;
- `final_selected_rank_post_refit_status`, `final_objective`,
  `final_stationarity_residual`, `final_max_envelope_ratio`, and
  `final_objective_stability_pass`.

The complete nested diagnostic JSON also stores every block's values through cap+1, all start
records, termination states, numerical ranks, and the explicit statement that observable
diagnostics do not certify the unknown global objective gap.

Revision-10 fields coexist with historical raw Revision-9 fields. IC and candidate fields are
null/not produced on the Revision-10 path rather than being repurposed.
