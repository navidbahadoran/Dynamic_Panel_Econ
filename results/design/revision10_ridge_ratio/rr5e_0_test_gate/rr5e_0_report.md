# RR5e.0 test-gate report

Date: 2026-08-11

## Frozen-state gate

- `HEAD`, local `main`, and `origin/main` were identical at `caa8572556a8d892ae07936fc4e6460c73819285` before diagnosis.
- The worktree was clean.
- No scientific file was touched.

## Diagnosis

All eight pytest errors were setup failures while pytest tried to recursively clear its configured `results/pytest-parallel-audit` basetemp. The final `os.rmdir` operation raised `PermissionError: [WinError 5]`. `.pytest_cache` independently produced permission warnings.

The paths were directories, untracked, ignored, and test-only. Their ACLs exposed full control only to `SYSTEM`, `Administrators`, and `OWNER RIGHTS`. There were no Python/pytest processes before cleanup, so the failure was a stale ACL/permission artifact rather than a live Windows handle.

## Repair and lifecycle conclusion

Only `results/pytest-parallel-audit` and `.pytest_cache` were removed. No source/test repair was warranted: after cleanup, the full suite passed twice, no workers survived either run, directory remove/recreate/remove succeeded, and a generated manifest could be renamed and restored. The evidence therefore does not establish an RR5d resource-lifecycle bug.

## Validation

- first `pytest -q`: `182 passed in 42.04s`;
- second `pytest -q`: `182 passed in 42.13s`;
- post-suite Python/pytest process count: zero;
- manifest rename out/back: passed;
- basetemp/cache recursive deletion: passed;
- Ruff: passed;
- `git diff --check`: passed.

## Scope confirmation

- no estimator, solver, selector, calibration, DGP, tolerance, start, manuscript, or scientific evidence changed;
- no test was weakened, skipped, marked xfail, or otherwise bypassed;
- no Monte Carlo or DGP generation ran;
- `RR5E_MASTER_SEED` was not inspected, reserved, derived from, or used;
- RR5e was not launched.

Decision: **RR5e.0-PASS**.
