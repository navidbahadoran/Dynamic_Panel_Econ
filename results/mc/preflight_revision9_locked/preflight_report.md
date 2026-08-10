# Locked Revision-9 statistical preflight

This is the single authorized 24-realization, 48-method-evaluation preflight. No penalty or threshold sensitivities were run and no tuning was performed.

## Execution accounting

- Fixed-rank evaluations: 24; selected-rank evaluations: 24.
- Unique matched semantic DGP realizations: 24; SHA-256 match: True.
- Frozen calibration SHA-256: `a80900898ff5bfef84380bd1cbd68a27d7f22c8ae3023b8b8c64ec6cec6f471e`.
- Global deterministic envelope: `8.288745227963506` <= `9`.

## Supplied-rank estimator

- Replication-level execution success: 0.833.
- N=50: 8/12 successful; N=100: 12/12 successful.
- N=50 target retention: point 144/216, inference 80/216; N=100: point 216/216, inference 216/216.
- Split fits: N=50 32/48; N=100 48/48. Every successful replication supplied exactly four split coefficient fits.
- Fixed split boundary-active fits: 9; fixed full-panel boundary-active fits: 0.
- Gram/Riesz: minimum empirical tangent-Gram eigenvalue 0.113834, maximum condition number 411.884, maximum Riesz residual 9.97317e-09; no invalid finite variances or SEs.
- Retention accounting reconciled: True.
- Runtime seconds: N=50 total 98.452, median 1.804, max 38.812; N=100 total 79.047.

## Selected-rank procedure

- Rank-result numerical completion: 1.000.
- Candidate coverage: 1.000; exact recovery: 0.000; cap-hit rate: 0.000.
- Complete selected-rank distribution: (0,0,0) in 24/24 replications; underselection 24/24 and overselection 0/24.
- P(true rank absent from candidates)=0; P(selected rank != truth | truth in candidates)=1.
- Selected split fits: N=50 48/48; N=100 48/48.
- Point retention is 432/432; inference retention is 0/432 because the selected zero tangent spaces make the requested targets unsupported.
- Runtime seconds: N=50 total 32180.154, median 2591.494, max 5041.653; N=100 total 34048.977, median 2377.046, max 7815.426.
- Assessment: **PAPER RANK-SELECTOR FINITE-SAMPLE NO-GO**.
- The complete distribution and candidate-level IC decomposition are stored separately.

## Medium-diagnostic decision

**NO-GO.** This decision is diagnostic and does not authorize tuning or a medium run.

## Reporting corrections

- Failure placeholder target rows now inherit semantic IDs and DGP hashes from the authoritative attempt ledger.
- Candidate Q_hat is reported as twice the numerical half-loss, matching the information-criterion implementation; invalid candidates retain infinite IC.
- The decision helper now treats exact rank recovery as a required selector diagnostic, in addition to numerical completion, coverage, and cap behavior.

## Commands

```powershell
.\.venv\Scripts\python.exe scripts\run_mc.py --config configs\mc\preflight_revision9_locked_fixed.toml --print-resolved-config --dry-run
.\.venv\Scripts\python.exe scripts\run_mc.py --config configs\mc\preflight_revision9_locked_selected.toml --print-resolved-config --dry-run
.\.venv\Scripts\python.exe scripts\run_mc.py --config configs\mc\preflight_revision9_locked_fixed.toml
.\.venv\Scripts\python.exe scripts\run_mc.py --config configs\mc\preflight_revision9_locked_selected.toml
.\.venv\Scripts\python.exe scripts\audit_revision9_locked_preflight.py --fixed-root results\mc\preflight_revision9_locked\fixed_rank\f87b622f889053fe --selected-root results\mc\preflight_revision9_locked\selected_rank\aa152561964a7ec3 --output-root results\mc\preflight_revision9_locked
```

Rank summary rows: 8; distribution rows: 8; Gram/Riesz rows: 264; split-summary rows: 32; boundary-summary rows: 64.
