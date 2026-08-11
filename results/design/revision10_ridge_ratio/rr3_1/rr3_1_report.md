# RR3.1 reviewed-manuscript lock report

## Decision

`RR3.1-PASS`

The author-reviewed Revision-10 manuscript is installed as the canonical scientific source for RR4. Historical Revision-9 and RR1/RR2/RR2.5/RR2-FREEZE/RR3 records were not changed.

## Source and canonical hash

- Supplied source SHA-256: `82714cdbee6c99f65b1fe35193b1a8910238057ffa1050826e1123c949a968c5`
- Canonical source: `manuscript/ver8_revision10_ridge_ratio_reviewed.tex`
- Canonical SHA-256: `d1c7224e0962df6a54fb92a550671f72b66258f30afd64d292f3bf60c7d87be9`

The hash change is limited to the requested Assumption 8 clarification, consistent use of $\calD_{\max}^{+}$ and $c_+$ in Appendix A.7, exact unresolved-status text, and replacement of unnecessary nuclear-norm phrasing by the standard fixed-rank operator--Frobenius inequality.

## Validation

- Three-pass noninteractive `pdflatex` compilation succeeded and produced 75 pages.
- All 29 citations resolve; there are no undefined citations or duplicate literal labels.
- The Revision-10 rank table and ridge-ratio diagnostic table resolve.
- Nine undefined references remain only because four pre-existing numerical MC table fragments and five pre-existing MC figure fragments were not supplied and were intentionally represented by empty validation stubs.
- No input file is missing from the validation build.
- There are no overfull boxes. Five underfull-box notices occur in unchanged empirical/bibliographic prose.
- Representative pages covering the abstract, selector, Assumptions 7--8, theorem/corollary, implementation, Monte Carlo design, table shells, conclusion/references, and the full A.7 proof were visually inspected without an RR3.1 layout defect.
- `git diff --check` passed before commit.
- No source/helper code changed, so pytest and Ruff were not run under the manuscript-only protocol.

No Monte Carlo, model fitting, selector implementation, DGP change, or tuning was performed.
