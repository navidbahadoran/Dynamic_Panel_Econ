# Locked Revision-9 no-go diagnosis

This is an offline analysis of already-computed preflight records. No simulation, fitting, rank selection, split fit, or Riesz solve was run.
Accepted source preflight commit: `9f0778d2503a66d913274ac41ae2b6e94257a5b9`.

## Candidate validity and locked IC

- N=50: truth present 12/12; valid truth post-refit 6/12; invalid 6/12.
- N=100: truth present 12/12; valid truth post-refit 12/12; invalid 0/12.
  - dgp1_N50_T50_r2026080902_truth1-1-1: `stationarity_high|objective_stability_failed`.
  - dgp2_N50_T50_r2026080902_truth1-1-1: `objective_stability_failed`.
  - dgp3_N50_T50_r2026080901_truth1-1-1: `objective_stability_failed`.
  - dgp3_N50_T50_r2026080902_truth1-1-1: `objective_stability_failed`.
  - dgp3_N50_T50_r2026080903_truth1-1-1: `stationarity_high|objective_stability_failed`.
  - dgp4_N50_T50_r2026080901_truth1-1-1: `stationarity_high|objective_stability_failed`.
- Maximum locked-IC reconstruction error over valid candidates: `5.82e-11`.
- At c=1, both sizes select (0,0,0) in all 12 saved replications.
- Best simultaneous saved-grid tradeoff occurs at c=3.4651874e-06: Exact=2/12 at N=50 and 1/12 at N=100; this is not substantial at both sizes.

## Especially requested pairwise break-even values

| N | Competitor | Available | Minimum c | Median c | Maximum c |
|--:|:--|--:|--:|--:|--:|
| 50 | (0,0,0) | 6 | 8.5939279e-06 | 1.1696778e-05 | 1.4723553e-05 |
| 50 | (0,0,1) | 6 | 8.5452077e-06 | 1.0706068e-05 | 1.2270306e-05 |
| 50 | (0,1,0) | 5 | 1.0406598e-05 | 1.3107964e-05 | 1.8122432e-05 |
| 50 | (0,1,1) | 6 | 9.6487461e-06 | 1.1878983e-05 | 1.5017899e-05 |
| 50 | (1,0,0) | 6 | 6.0424333e-06 | 6.3970637e-06 | 7.5742597e-06 |
| 50 | (1,0,1) | 2 | 5.0313584e-06 | 6.3220601e-06 | 7.6127617e-06 |
| 50 | (1,1,0) | 5 | 7.095034e-06 | 8.3641641e-06 | 8.9048643e-06 |
| 100 | (0,0,0) | 12 | 6.1930884e-06 | 6.57422e-06 | 7.4282455e-06 |
| 100 | (0,0,1) | 12 | 4.9298633e-06 | 6.1690795e-06 | 6.7720746e-06 |
| 100 | (0,1,0) | 12 | 6.7953307e-06 | 7.9322497e-06 | 9.6701513e-06 |
| 100 | (0,1,1) | 12 | 4.9927182e-06 | 6.5806757e-06 | 8.5613086e-06 |
| 100 | (1,0,0) | 12 | 3.1210451e-06 | 3.3668398e-06 | 4.3103103e-06 |
| 100 | (1,0,1) | 12 | 2.3018158e-06 | 2.739486e-06 | 3.739842e-06 |
| 100 | (1,1,0) | 12 | 4.0102077e-06 | 4.6447909e-06 | 6.0954854e-06 |

## Common-c result

- N50: common [4.5656569e-06, 5.0313584e-06], nonempty=True; maximum overlap 6/6.
- N100: common [1.2340838e-06, 2.3018158e-06], nonempty=True; maximum overlap 12/12.
- pooled_N50_N100: common [4.5656569e-06, 2.3018158e-06], nonempty=False; maximum overlap 12/18.
- Classification: **CASE C** -- truth-optimal fixed-c requirements are incompatible across sizes.

## Baseline offline rank behavior

- N=50, c=1: Exact=0, Under only=12, Over only=0, Mixed=0; distribution {"(0,0,0)": 12}.
- N=100, c=1: Exact=0, Under only=12, Over only=0, Mixed=0; distribution {"(0,0,0)": 12}.

## Supplied-rank failures and split interiority

- The four N=50 full-fit failures contain 12 starts. Every start converged, remained feasible, retained rank (1,1,1), and failed only because its residual was slightly above the locked `1e-6` stationarity threshold. The primary classification is KKT/optimality failure, not start disagreement.

| Failed semantic ID | Residual range | Maximum objective gap | Maximum coefficient | Classification |
|:--|:--|--:|--:|:--|
| dgp1_N50_T50_r2026080902_truth1-1-1 | 1.1156029e-06 to 1.1729656e-06 | 1.2848559e-08 | 5.7240132 | KKT/optimality failure |
| dgp3_N50_T50_r2026080902_truth1-1-1 | 1.009091e-06 to 1.3069844e-06 | 4.0261617e-09 | 6.8131932 | KKT/optimality failure |
| dgp3_N50_T50_r2026080903_truth1-1-1 | 1.3288969e-06 to 1.3880209e-06 | 2.874134e-09 | 7.206605 | KKT/optimality failure |
| dgp4_N50_T50_r2026080901_truth1-1-1 | 1.0839503e-06 to 1.2378604e-06 | 9.9499615e-09 | 8.2209114 | KKT/optimality failure |

- Boundary-active split fits: 9.
- Successful N=50 fixed-rank replications with zero/one/multiple active splits: 3/3/2.
- Five successful replications have at least one boundary-active split and all eight broad targets in each are suppressed by `boundary_interiority_failure`. The other three have no boundary-active split but fail the locked split numerical checks, producing `split_fit_failure`. Hence broad inference is 0/96.
- The locked records do not persist the active matrix identity or entry-level active-set count. Those fields are explicitly marked unavailable; they are not inferred.

## Fixed-positive-multiplier theorem implication

For any fixed `c_kappa > 0`, the ratio of `zeta_NT` to `kappa_NT (N+T)/(NT)` is `1/{c_kappa log(NT)} -> 0`. A fixed positive multiplier also leaves the maintained vanishing-rate condition unchanged. Therefore c=1 is a normalization rather than a constant uniquely implied by the proof. Any revision must nevertheless freeze the multiplier before a fresh independent validation experiment; post-outcome selection would be inappropriate.

## Statistical versus numerical decomposition

- Supplied rank: the N=100 estimator/inference path is fully retained; N=50 has four numerical full-fit failures and universal broad-target split-interiority/numerical failure. Gram/Riesz diagnostics are reliable wherever reached.
- Selected rank: candidate coverage is complete, but only 6/12 N=50 true-rank post-refits are numerically eligible. At N=100 all 12 are eligible. Conditional on eligibility, the locked IC penalty scale drives selection to zero rank at c=1, which removes tangent support and suppresses all selected-rank inference.
- Paper-level issue: the finite-sample normalization of the paper IC. Numerical issue for the same estimator: N=50 stationarity acceptance and split interiority reliability. Neither is changed here.

## Next paper-level decision

The author should decide whether Revision 9 should retain its c=1 normalization or be revised to pre-specify a different fixed positive multiplier and a fresh validation protocol. This audit does not select a constant. Independently, the author should decide whether the N=50 numerical/interiority evidence is acceptable for the claimed design or requires a paper-level change in supported sample sizes or numerical assumptions before any code change.

No medium or production run is authorized by this report.

## Penalty scale snapshot

| N | b_NT | kappa_base | 0->1 | 1->2 | 2->3 | truth-zero |
|--:|--:|--:|--:|--:|--:|--:|
| 50 | 15.017342 | 845107.02 | 33466.238 | 32790.152 | 32114.067 | 100398.71 |
| 100 | 19.843077 | 2833483.9 | 56386.329 | 55819.633 | 55252.936 | 169158.99 |
| 200 | 25.625536 | 8279760.7 | 82590.613 | 82176.625 | 81762.637 | 247771.84 |
| 400 | 32.526674 | 21813797 | 108932.65 | 108659.98 | 108387.3 | 326797.95 |
