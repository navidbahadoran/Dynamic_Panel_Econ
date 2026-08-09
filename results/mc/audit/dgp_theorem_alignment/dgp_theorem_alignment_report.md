# DGP/theorem alignment report

## Outcome

The proposed calibration design is theorem-compatible at the DGP layer and has
a finite deterministic common coefficient envelope.  It is prepared but not
activated.  The medium diagnostic remains cancelled, and no estimator,
rank-selection, inference, or production configuration was changed.

## Frozen design

The candidate table contains 52 ex-ante cells:

- baseline DGPs 1--4 at N=T in {50,100,200,400};
- rank-stress DGPs 1--4 at N=T in {50,100,200} for `(1,1,1)`, `(2,1,1)`, and
  `(1,0,2)`.

It uses analytical population `c_h` and a separate 50-draw calibration-only
experiment for frozen `c_xi(DGP,N,T,rank design)`.  Identified cells hit 0.65 in
the calibration experiment to root tolerance.  For `(1,0,2)`, `c_xi=1` and the
induced-R2 ranges are:

| DGP | Induced R2 range across N=T=50,100,200 |
|---:|---:|
| 1 | 0.528614--0.539301 |
| 2 | 0.530349--0.531442 |
| 3 | 0.532223--0.535142 |
| 4 | 0.521091--0.536590 |

Population `pi_H` is exactly 0.30 in every cell.  Finite independent calibration
draws fluctuate around that target, but the fluctuation no longer changes H.

## Envelopes and proposed box

Every row's `C_A`, `C_beta`, `C_H`, and `C_Theta` is finite and deterministic.
The maximum is `8.410761115894578`, attained by DGP 4, N=T=50, rank-stress
`(1,1,1)`, whose frozen `c_xi` is `2.472531361047685`.

The current pair `B=9,c_B=1` fails because its allowed truth envelope is 8.
The reported common proposal is `B=10,c_B=1`.  This covers every candidate cell
with deterministic margin `1.589238884105422`.  Neither B nor the candidate
table has been activated in an execution config.

## DGP 4 group structure

The calibration change scales only H and u.  It has no effect on A, B, the
stability rescaling `c_a`, or any realized A/B group difference for a fixed
semantic primitive draw.

With loading means 0.9 and 1.1 and `E[f_a]=0.5`, the raw population A group
means at a date are 0.45 and 0.55, with group-2-minus-group-1 difference 0.10.
The common realized stability scale multiplies both raw group targets; its
realized postscale difference is unchanged by this calibration redesign.

With slope loading means 0.8 and 1.2 and `E[f_b]=0.6`, the population B group
means are 0.48 and 0.72, with difference 0.24.  These are unchanged.

The maintained `kappa_f_b=0.20` gives slope-factor support `[0,1.2]`, so B-entry
still lacks a deterministic positive lower time-leverage bound.  That separate
limitation was not altered or combined with this task.

## Activation and next task

The candidate DGP is ready for author review.  After approval, activation would
require pointing the intended execution configs to the frozen table and setting
the common coefficient bound to the approved value.  Only after that should the
literal box-constrained estimator alignment proceed.  No medium diagnostic
should restart before both approvals and the subsequent estimator correction.

