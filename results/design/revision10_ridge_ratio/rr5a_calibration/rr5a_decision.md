# RR5a decision

Decision: `RR5a-CALIBRATION-FROZEN`.

All 24 authorized rectangular calibration entries are present. The frozen seed, 50-draw protocol, analytical `c_H`, identified-cell pooled-R² target, common random numbers, and zero-slope `c_xi=1` normalization were preserved. All 52 prior calibration entries remain byte/value unchanged, no rank-selection evidence was generated, and no scientific source or manuscript was modified.

The authorized historical-invariance test now verifies the exact 52-entry Revision-9 subset and its frozen key/value digest, the exact 24-entry Revision-10 rectangular subset, and the complete 76-entry union. All 167 tests pass; Ruff and `git diff --check` pass.

This decision freezes calibration only. It does not launch or approve execution of RR5 within this task.
