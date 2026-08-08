# Offline IC break-even diagnostic

This report reconstructs the Revision-8 IC from saved preflight candidate fits. It performs no fitting, simulation, or estimator changes.

## Saved-candidate audit

- Successful replications: 8.
- Every persisted candidate record contains rank, objective, normalized residual variance, model dimension, validity, source attribution, and baseline IC.
- Candidate-record counts equal the saved final locally completed candidate counts in all eight replications.
- Replications with a valid true-rank post-refit: 7.
- Replications with an invalid true-rank post-refit: 1.
- DGP 3, replication 0 has the true rank in the screened candidate set, but it is excluded from IC competition because objective stability failed.

## Exact interval summary

- Nonempty true-rank intervals among eligible replications: 7/7.
- Median lower endpoint: 3.6332e-06.
- Median finite upper endpoint: 6.09198e-06.
- Smallest finite upper endpoint: 4.11163e-06.
- Maximum lower endpoint: 4.20548e-06 (DGP 2, replication 0).
- Minimum upper endpoint: 4.11163e-06 (DGP 3, replication 2).
- Common intersection: empty; at most 6/7 eligible intervals overlap. DGP 2, replication 0 sets the largest lower bound, while DGP 3, replication 2 sets the smallest upper bound.

The lower endpoint is imposed by higher-dimensional competitors; the upper endpoint is imposed by lower-dimensional competitors. Bounds are the intersection of all valid saved candidates, not only the all-zero model.

## Diagnostic grid

| c_kappa | true rank | underfit | overfit | selected-rank distribution |
|---:|---:|---:|---:|:---|
| 1e-08 | 1/7 (14.3%) | 1/8 (12.5%) | 6/8 (75.0%) | (1,0,1): 1, (1,1,1): 1, (1,1,2): 4, (1,2,1): 2 |
| 3e-08 | 1/7 (14.3%) | 1/8 (12.5%) | 6/8 (75.0%) | (1,0,1): 1, (1,1,1): 1, (1,1,2): 4, (1,2,1): 2 |
| 1e-07 | 1/7 (14.3%) | 1/8 (12.5%) | 6/8 (75.0%) | (1,0,1): 1, (1,1,1): 1, (1,1,2): 4, (1,2,1): 2 |
| 3e-07 | 1/7 (14.3%) | 1/8 (12.5%) | 6/8 (75.0%) | (1,0,1): 1, (1,1,1): 1, (1,1,2): 4, (1,2,1): 2 |
| 1e-06 | 1/7 (14.3%) | 1/8 (12.5%) | 6/8 (75.0%) | (1,0,1): 1, (1,1,1): 1, (1,1,2): 4, (1,2,1): 2 |
| 3e-06 | 3/7 (42.9%) | 1/8 (12.5%) | 4/8 (50.0%) | (1,0,1): 1, (1,1,1): 3, (1,1,2): 4 |
| 1e-05 | 0/7 (0.0%) | 8/8 (100.0%) | 0/8 (0.0%) | (1,0,0): 8 |
| 3e-05 | 0/7 (0.0%) | 8/8 (100.0%) | 0/8 (0.0%) | (0,0,0): 8 |
| 1e-04 | 0/7 (0.0%) | 8/8 (100.0%) | 0/8 (0.0%) | (0,0,0): 8 |
| 3e-04 | 0/7 (0.0%) | 8/8 (100.0%) | 0/8 (0.0%) | (0,0,0): 8 |
| 1e-03 | 0/7 (0.0%) | 8/8 (100.0%) | 0/8 (0.0%) | (0,0,0): 8 |

The true-rank count uses only numerically valid true-rank post-refits; its maximum is therefore seven, not eight. Underfit and overfit are separate componentwise indicators, so a mixed rank vector can satisfy both.

## Shortlist for later testing

- `c_kappa = 3e-06`: true rank in 3/7 eligible replications; underfit in 1/8 and overfit in 4/8 saved replications.

These are diagnostic shortlist values only. The active configuration remains unchanged; a later lower-cap preflight must be approved and executed before choosing a production constant.

## Literal Revision-8 penalty

The base rate is `kappa_NT = b_NT^2 log(NT)^(d_s+3)`, with `b_NT = (NT)^(1/(8+eta)) log(NT)`, `eta=4`, and `d_s=1`. The penalty table reports increments for raising one matrix from rank r to r+1 at N=T in {50,100,200,400}.

| N=T | NT | log(NT) | b_NT | base kappa | delta d1 | base increment 1 |
|---:|---:|---:|---:|---:|---:|---:|
| 50 | 2500 | 7.82405 | 15.0173 | 845107 | 99 | 33466.2 |
| 100 | 10000 | 9.21034 | 19.8431 | 2.83348e+06 | 199 | 56386.3 |
| 200 | 40000 | 10.5966 | 25.6255 | 8.27976e+06 | 399 | 82590.6 |
| 400 | 160000 | 11.9829 | 32.5267 | 2.18138e+07 | 799 | 108933 |

First-rank IC increments under each diagnostic multiplier are tabulated in `tab_ic_penalty_magnitude.csv`; the LaTeX artifact additionally provides separate panels for rank increases 0-to-1, 1-to-2, and 2-to-3.

## Zero-versus-truth break-even distribution

Among the 7 replications with valid true-rank candidates: minimum 1.05634e-05; 10th percentile 1.0925e-05; median 1.19992e-05; 90th percentile 1.25771e-05; maximum 1.32427e-05.
