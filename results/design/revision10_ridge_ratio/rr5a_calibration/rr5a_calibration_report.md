# RR5a rectangular calibration freeze report

## Outcome

RR5a completed and froze the 24 missing rectangular rank-stress calibration cells and did not launch RR5. The calibration, historical-invariance, and repository validation gates pass. The decision is `RR5a-CALIBRATION-FROZEN`.

The only scientific configuration change is the addition of the required records for DGPs 1–4, dimensions `(100,200)` and `(200,100)`, and true ranks `(1,1,1)`, `(2,1,1)`, and `(1,0,2)`. No selector, DGP, inference, manuscript, or existing calibration record was changed.

## Frozen protocol confirmed

- The master seed is `20260807`; every cell uses 50 calibration-only draws.
- Streams are derived by `SeedSequence` from semantic keys: master seed, `rank_stress_calibration`, DGP, `N`, `T`, the full rank tuple, and calibration replication. Adding cells therefore cannot shift any previously frozen stream.
- Each draw uses the actual rescaled rank-stress matrices with component strengths `(1,1)`.
- `c_H` is analytical: population disturbance variance and the rank-specific population variance of raw `H` are inserted into the frozen `pi_H=0.30` formula. Rank-one-H cells have `c_H=0.6546536707079772`; rank-two-H cells have `c_H=0.7301712917987002`.
- For identified cells, the raw draws are held fixed while the recursion is decomposed as `y(c_xi)=y_base+c_xi*y_scale`. The existing 80-point geometric grid brackets the first positive root, and `brentq` uses `xtol=rtol=1e-10` to target pooled R² 0.65.
- When `r_B=0`, scale is not identified. The frozen normalization `c_xi=1` is retained and the induced pooled R² is reported; it is not classified as a calibration failure.

## Calibration results

| DGP | N | T | true ranks | c_H | c_xi | pooled R² | status |
|---:|---:|---:|:---:|---:|---:|---:|:---|
| 1 | 100 | 200 | (1,1,1) | 0.654654 | 1.768546 | 0.650000 | root converged |
| 1 | 100 | 200 | (2,1,1) | 0.654654 | 0.818422 | 0.650000 | root converged |
| 1 | 100 | 200 | (1,0,2) | 0.730171 | 1.000000 | 0.540172 | normalization |
| 1 | 200 | 100 | (1,1,1) | 0.654654 | 2.016507 | 0.650000 | root converged |
| 1 | 200 | 100 | (2,1,1) | 0.654654 | 0.831243 | 0.650000 | root converged |
| 1 | 200 | 100 | (1,0,2) | 0.730171 | 1.000000 | 0.533333 | normalization |
| 2 | 100 | 200 | (1,1,1) | 0.654654 | 1.887928 | 0.650000 | root converged |
| 2 | 100 | 200 | (2,1,1) | 0.654654 | 0.823473 | 0.650000 | root converged |
| 2 | 100 | 200 | (1,0,2) | 0.730171 | 1.000000 | 0.536218 | normalization |
| 2 | 200 | 100 | (1,1,1) | 0.654654 | 1.939554 | 0.650000 | root converged |
| 2 | 200 | 100 | (2,1,1) | 0.654654 | 0.815733 | 0.650000 | root converged |
| 2 | 200 | 100 | (1,0,2) | 0.730171 | 1.000000 | 0.535364 | normalization |
| 3 | 100 | 200 | (1,1,1) | 0.654654 | 2.141694 | 0.650000 | root converged |
| 3 | 100 | 200 | (2,1,1) | 0.654654 | 0.870469 | 0.650000 | root converged |
| 3 | 100 | 200 | (1,0,2) | 0.730171 | 1.000000 | 0.535063 | normalization |
| 3 | 200 | 100 | (1,1,1) | 0.654654 | 2.241495 | 0.650000 | root converged |
| 3 | 200 | 100 | (2,1,1) | 0.654654 | 0.874451 | 0.650000 | root converged |
| 3 | 200 | 100 | (1,0,2) | 0.730171 | 1.000000 | 0.529296 | normalization |
| 4 | 100 | 200 | (1,1,1) | 0.654654 | 1.988023 | 0.650000 | root converged |
| 4 | 100 | 200 | (2,1,1) | 0.654654 | 0.842478 | 0.650000 | root converged |
| 4 | 100 | 200 | (1,0,2) | 0.730171 | 1.000000 | 0.514476 | normalization |
| 4 | 200 | 100 | (1,1,1) | 0.654654 | 2.028576 | 0.650000 | root converged |
| 4 | 200 | 100 | (2,1,1) | 0.654654 | 0.848070 | 0.650000 | root converged |
| 4 | 200 | 100 | (1,0,2) | 0.730171 | 1.000000 | 0.525465 | normalization |

Full precision, root brackets/residuals, population and realized H shares, deterministic and realized envelopes, and first/last stream states are in `rr5a_rectangular_calibration.csv`.

## Seed and integrity audits

The 1,200 new draw streams are unique. Their generated 128-bit states have no collision with the 2,600 streams represented by the existing frozen cells, the 24 locked Revision-9 preflight DGP streams, the 180 reserved RR5 DGP streams, or the explicit deterministic test streams. Per-draw evidence is in `rr5a_seed_audit.csv`.

The original frozen config SHA-256 was `a80900898ff5bfef84380bd1cbd68a27d7f22c8ae3023b8b8c64ec6cec6f471e`. Its complete bytes remain the exact prefix of the updated file, and all 52 old parsed entries are value-for-value identical. The updated SHA-256 is `51b4ea76fdd2295a2cfeb840ef7b2b3f468d32c45ad2c123523363d5b968f17e`. Per-entry hashes are in `rr5a_existing_entry_integrity.csv`.

The largest new deterministic coefficient envelope is 7.624850556803478, which satisfies the frozen simulation interior condition `C_Theta <= B-c_B = 9`.

## Scope exclusions

No cap+1 pilot, ridge ratio, selected rank, recovery classification, post-refit, Riesz object, target estimate, confidence interval, or inference result was computed or saved. RR5 seed `2026081101` remains reserved and unused for calibration.

## Validation

Before calibration, 167 tests passed, Ruff passed, and `git diff --check` passed. After the authorized historical-invariance test repair, all 167 tests pass, Ruff passes, and `git diff --check` passes.

The repaired test identifies the exact 52 historical Revision-9 keys, verifies an order-independent SHA-256 digest of every historical key/value mapping against the approved pre-RR5a Git configuration, identifies the exact 24 Revision-10 rectangular keys, and asserts that the complete frozen key set is their 76-entry union. No expected historical calibration number was changed.
