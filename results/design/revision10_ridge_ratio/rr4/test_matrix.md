# RR4 deterministic test matrix

`tests/test_revision10_rank_selection.py` covers:

- exact normalized spectral formula and retention through cap+1;
- exact ridge `1/log(NT)`;
- rank-zero anchor and a selected zero rank;
- unique positive-rank selection;
- genuine cap+1 numerator at the reporting cap;
- first-index exact-tie behavior;
- separate block selection and `(1,0,2)` assembly;
- covariate and outcome-unit scale-equivariance algebra;
- uncentered full-sample RMS weights, `w_H=1`, and no reference floor;
- natural pilot numerical-rank collapse below cap+1;
- nonfinite, stationarity-failed, and feasibility-failed pilot fixtures;
- pilot-to-ratios-to-one-final-post-refit call sequence;
- no final/fallback call after unresolved pilot;
- deterministic serial/spawn-process semantic equality;
- Revision-10 default routing and exclusion of legacy tuning fields.

`tests/test_audit.py` explicitly routes historical candidate/IC and cap-hit regression tests to
`select_ranks_revision9`. Existing supplied-rank and parallel tests remain part of the full suite.

No RR4 test inspects rank recovery in DGP 1--4 or selects/tunes a scientific constant.
