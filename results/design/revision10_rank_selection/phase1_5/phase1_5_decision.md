# Phase-1.5 decision

Classification: **ST3 — THEORY COST TOO HIGH** for the requested self-tuning route at present.

Progress is material:

- raw, normal, operator, local-gain, and post-refit objects are separated exactly;
- the local gain has a symbolic score-squared-over-curvature upper bound;
- mixed-state contamination is shown to invalidate the original greedy proof;
- nuisance-at-cap block profiling supplies an order-invariant mixed-state remedy;
- fixed caps reduce simultaneous multiplicity to `L=sum_M cap_M`;
- a precise rectangular-Freedman/spatial-variance candidate envelope is identified.

The remaining envelope is not a small missing constant. It requires a new theorem for an
observable, uniformly conservative growing operator-valued predictable variance based on
cap-pilot residuals under conditional spatial mixing/NED, plus residual replacement and
truncation. Revision 9 proves only an operator-score `O_p` rate and scalar fixed-direction HAC.
Neither finite union bounds nor ordinary residual MSE fills this gap. A bootstrap would require a
comparably substantial new validity theory.

The design still uses only two *candidate* new primitive conditions, but N1 cannot credibly be
assumed away. Therefore there is no exact tuning-free boundary, no authorized selector, and no
Phase-2/3 theorem or implementation authorization. No iid, Gaussian, independence, or
homoskedastic shortcut is accepted.
