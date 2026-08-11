# RR4 report

RR4 implements the frozen Revision-10 blockwise ridge-ratio selector as a separate primary path.
The implementation uses one joint at-most-cap+1 literal box-constrained pilot, exact fitted-value
RMS normalization, the fixed ridge and anchor, separate smallest-index block argmins, and one final
literal selected-rank fixed-rank post-refit.

Revision-9 IC/path/candidate logic remains preserved behind `revision9_ic` and is absent from the
primary option set. Unresolved pilots have no fallback and suppress primary selected-rank point and
inference output. Supplied-rank mode is unchanged.

Validation completed with `167 passed` under `pytest -q`. `ruff check src scripts tests` and
`git diff --check` passed. The serial/spawn-process test compares the actual Revision-10 pilot,
spectra, normalized values, ratios, selected ranks, final post-refit objective and coefficients;
identities/ranks/statuses are exact and numerical outputs agree at `rtol=atol=1e-12`.

The exact RR4 commit, push, and synchronization result are recorded in the task handoff.

No Monte Carlo experiment, DGP preflight, empirical fit, rank-recovery inspection, or scientific
tuning was performed in RR4.
