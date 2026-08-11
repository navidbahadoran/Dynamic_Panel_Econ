# RR5d engineering freeze record

Frozen implementation invariants:

- Revision-10 reporting/pilot caps `(3,3,3)`/`(4,4,4)`;
- `B=10`;
- interior projected-gradient `1e-6`;
- constrained factor-space KKT `1e-4`;
- three maintained starts and objective agreement `1e-6`;
- unchanged normalized spectrum, ridge `1/log(NT)`, rank-zero anchor and smaller-rank exact tie rule;
- unchanged literal final selected-rank post-refit;
- isolated cap+1 engineering path; supplied-rank and Revision-9 legacy solver unchanged;
- no clipping, regularization, singular-value rank threshold, factor floor or approximate statistical SVD;
- fingerprinted one-task atomic durability before legacy chunk output;
- deterministic semantic aggregation and one native thread per outer worker;
- recommended default `n_jobs=8` from deterministic computation only.

Validation must include the complete Pytest suite, Ruff over `src scripts tests`, and `git diff --check`. RR5, DGP1-DGP4 scientific Monte Carlo, rank-recovery inspection, inference and manuscript editing remain prohibited after this freeze unless separately authorized.

Final validation: `182 passed`; Ruff passed; `git diff --check` passed.
