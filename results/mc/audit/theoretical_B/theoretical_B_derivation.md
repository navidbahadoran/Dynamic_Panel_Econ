# Deterministic coefficient-envelope audit

## Decision

**NOT VERIFIED.**  Under the maintained calibration mapping there is no finite
deterministic upper bound for `c_h`.  Consequently there is no finite uniform
bound for `H_0 = c_xi c_h H_raw`, no finite common `C_Theta`, and no analytical
verification that the current `B = 9` and `c_B_sim = 1` satisfy
`C_Theta <= B - c_B_sim = 8`.

This is a support calculation, not an empirical calculation.  No Monte Carlo
draws, empirical maxima, quantiles, or calibration simulations were used.
Following the task's mandatory STOP rule, no optimizer was changed and no
coefficient-bound replay was run.

## Bounded-AR primitive

For the recursion

`g_t = rho g_{t-1} + sqrt(1-rho^2) z_t`, `g_0 = 0`,
`|z_t| <= sqrt(3)`,

the infinite-horizon deterministic envelope is

`G(rho) = sqrt(3) sqrt(1-rho^2) / (1-|rho|)`.

All maintained coefficient factors use `rho = 0.5`, so `G(0.5) = 3`.
Thus

- `|f_a| <= |0.5| + |0.1|*3 = 0.8`;
- `|f_b| <= |0.6| + |0.2|*3 = 1.2`;
- `|H_raw| <= sqrt(3)*3 = 3 sqrt(3) = 5.196152422706632`.

## Autoregressive coefficient A_0

For DGPs 1--3, `|lambda_a| <= 1 + 0.1 sqrt(3)`, hence

`C_A_raw = (1 + 0.1 sqrt(3))*0.8 = 0.9385640646055102`.

For DGP 4, the group means remain 0.9 and 1.1 and
`|lambda_a| <= max(0.9,1.1) + 0.08 sqrt(3)`, hence

`C_A_raw = (1.1 + 0.08 sqrt(3))*0.8 = 0.9908512516844082`.

Writing `m = max|A_raw|`, the code applies
`c_a = min(1, 0.85/m)`.  For every `m > 0` this gives exactly

`max|A_0| = c_a m = min(m,0.85) <= 0.85`.

Therefore `C_A = 0.85` for all baseline and rank-stress designs.  The baseline
implementation does not guard the degenerate `m=0` case before dividing by
`m`; the rank-stress implementation does.  The event has probability zero
under the continuous draws, but a literal total deterministic mapping should
use the rank-stress guard (`c_a=1` when `m=0`).  This observation does not cause
the failed H-envelope result.

## Slope coefficient B_0

For DGPs 1--3,

`C_BETA = (1 + 0.4 sqrt(3))*1.2 = 2.031384387633061`.

For DGP 4,

`C_BETA = (max(0.8,1.2) + 0.25 sqrt(3))*1.2`
`= 1.959615242270663`.

For rank stress `(1,0,2)`, the slope matrix is identically zero and its bound
is zero.  The other maintained rank vectors retain the corresponding baseline
bound.

## Higher-rank deterministic rescaling

Each added rank component has entry envelope `3*|strength|`.  The code uses the
common scale

`s_r(C) = C / (C + sum_{j=1}^{r-1} 3|strength_j|)`.

Consequently the triangle inequality gives
`s_r(C)[C + sum 3|strength_j|] = C`: the rank-two A and H stress matrices have
the same prespecified raw envelopes as the rank-one matrices.  With maintained
rank-two strength one, the scales are:

- A, DGPs 1--3: `0.2383010785682161`;
- A, DGP 4: `0.2482806772781126`;
- H, every DGP: `0.6339745962155613`.

This rescaling preserves the raw envelope, but it cannot bound the later
calibration multiplier.

## Calibration multipliers and H_0

For cells with positive slope rank, the successful `c_xi` root is searched on
the deterministic grid interval `[1e-4, 1e3]` and solved inside the first
bracket.  Hence successful identified calibrations satisfy `c_xi <= 1000`.
When the slope rank is zero, the code fixes `c_xi = 1`.  The diagnostic
evaluation at `1e8` is not a selected calibration value.

The other multiplier is

`c_h = sqrt([pi_H/(1-pi_H)] * var(u_tilde)/var(H_raw))`, with `pi_H=0.30`.

This mapping has no finite deterministic upper bound over admissible primitive
draws:

1. `u_tilde` contains Gaussian innovations, so its empirical calibration
   variance has unbounded support.
2. The bounded H primitives can be arbitrarily close to zero across the finite
   calibration panel, so `var(H_raw)` has no deterministic positive lower
   bound.  At zero the calibration rejects the draw; arbitrarily near zero it
   still produces arbitrarily large `c_h`.

It follows that

`C_H = 3 sqrt(3) * sup|c_xi c_h| = infinity`

for every maintained design with positive H rank.  The failure is caused by
`c_h`, even though `c_xi` itself has the deterministic restriction above.

## Overall envelope and status of B

`C_Theta = max(C_A, C_BETA, C_H) = infinity` under the current calibration
mapping.  Equivalently, no finite uniform `C_Theta` can presently be
established.  The minimum theoretical requirement is a finite fixed `B`
satisfying `B >= C_Theta + c_B`; because `C_Theta` is unavailable, there is no
finite numerical minimum to report and merely increasing B cannot verify the
assumption.

`B` is not estimated from a sample.  It is an ex-ante fixed parameter-space
constant.  There is no unique B implied by the theorem.  A finite choice can be
justified only after the maintained DGP/calibration construction itself yields
a finite deterministic coefficient envelope.  No calibration cap or revised
DGP is proposed here because the task expressly forbids inventing one.

