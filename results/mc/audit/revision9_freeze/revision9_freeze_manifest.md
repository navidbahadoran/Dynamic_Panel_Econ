# Revision 9 freeze manifest

Status: **CLOSED NO-GO benchmark**. Frozen from published `main` at
`b4bad7bdca3fc9d68ae3befa007d7fd69dccf549`.

This freeze preserves, without reinterpretation or tuning, the locked preflight, CASE C fixed-`c_kappa` diagnosis, CASE P3 cap-pilot diagnosis, and parallel-execution validation. Revision 9 selected `(0,0,0)` in 24/24 locked realizations although the true rank was in 24/24 candidate sets. The N=50 and N=100 truth-optimal fixed-multiplier intervals do not intersect. The thresholded cap pilot recovered the truth in 0/24 and was over-ranked in 23/24. At least one nuclear-path proposal was the truth in 24/24; screening and candidate coverage were therefore not the primary failure.

The supplied-rank estimator was clean at N=100. N=50 retains four narrow full-panel stationarity failures, nine boundary-active split fits, and broad inference failure. It is provisionally retained as a **small-sample stress design**, not deleted.

Outer parallel execution is scientifically equivalent across one and four workers. The deterministic scaling audit recommends the machine-specific setting `--n-jobs 12`.

No Revision-9 result is reopened here. In particular, this freeze does not replace `c_kappa`, retune `tau_NT`, promote either preliminary estimator to the final selector, modify code or manuscript, or authorize Monte Carlo.

## Evidence commits

- Locked preflight: `9f0778d2503a66d913274ac41ae2b6e94257a5b9`.
- CASE C diagnosis: `0606b3a7653581f18571dfb993f51b0e8b92efb5`.
- Parallel validation: `f3147a8ca9ca8161e1ccc07eef37049b27912a2e`.
- CASE P3 diagnosis: `b4bad7bdca3fc9d68ae3befa007d7fd69dccf549`.

The detailed ledger is `revision9_evidence_ledger.csv`; the permitted use of this evidence is governed by `revision9_to_revision10_firewall.md`.
