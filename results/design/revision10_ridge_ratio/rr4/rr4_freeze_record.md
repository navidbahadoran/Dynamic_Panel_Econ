# RR4 implementation freeze record

- Scientific source commit: `47af37560d6e0752519d6a6267a152871d8a5157`.
- Canonical committed manuscript:
  `results/design/revision10_ridge_ratio/rr3_1/manuscript/ver8_revision10_ridge_ratio_reviewed.tex`.
- Canonical SHA-256: `d1c7224e0962df6a54fb92a550671f72b66258f30afd64d292f3bf60c7d87be9`.
- Selector method: `revision10_ridge_ratio`.
- Reporting caps: `estimation.rank_caps`.
- Pilot caps: deterministically `reporting cap + 1` for every block.
- Box and numerical fields: existing `coefficient_bound`, optimizer tolerances,
  `start_objective_stability_tol`, and maintained multistart settings.
- Formula functions: listed in `implementation_map.md` and audited in
  `selector_formula_audit.md`.
- Diagnostics: listed in `result_schema.md`.
- Failure behavior: listed in `numerical_failure_rules.md`.
- Tests: listed in `test_matrix.md`.
- Legacy disposition: `legacy_revision9_disposition.md`.

The committed manuscript blob was independently archived from the approved commit and hashed
before implementation. The pre-existing working-tree rename of that manuscript was not modified
or included in RR4.

No DGP, B, interior simulation margin, target, Riesz, split, spatial-variance, or manuscript file
is part of this implementation freeze.
