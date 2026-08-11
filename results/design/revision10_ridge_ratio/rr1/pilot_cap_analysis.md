# Spectral-pilot cap analysis

## Mechanical truncation

A rank-at-most-`bar_r_M` matrix has `sigma_(bar_r_M+1)=0` mechanically. The ratio at the reporting
cap would therefore use an artificial zero numerator and could create a cap minimum regardless of
the data. A separate spectral pilot is required.

## Enlarged fixed cap

Define pilot caps `bar_r_M^pil=bar_r_M+1` for every block, while retaining the reporting set
`{0,...,bar_r_M}`. Let `Theta_hat^(cap+1)` minimize the same normalized loss and literal box
constraint over these fixed enlarged caps.

The Revision-9 cap-pilot proof extends verbatim:

1. `r_0 in R_max` makes truth feasible under the enlarged caps.
2. Every difference matrix has rank at most `bar_r_M+1+r_M,0`, still a fixed constant.
3. The existing uniform score and empirical prediction lower bounds apply to every fixed
   rank-capped difference collection.
4. If the computed normalized-loss objective gap relative to the enlarged-cap infimum is
   `delta_NT=o_p(zeta_NT)`, the basic inequality gives
   `||Theta_hat^(cap+1)-Theta_0||_F=O_p(sqrt(NT zeta_NT))`.
5. Therefore, uniformly over the fixed block set,
   `||M_hat^(cap+1)-M_0||_op=O_p(sqrt(NT zeta_NT))`.

No DGP, dependence, signal, or relative-growth assumption changes.

## Alternatives

- **A: cap+1 rank-at-most pilot — preferred.** Its operator rate follows directly from the
  existing proof and it exposes exactly the one extra singular value needed at the reporting cap.
- **B: nuclear/full-spectrum pilot.** It supplies a full spectrum, but Revision 9 does not prove
  the required uniform operator rate for every nuclear-path estimate. Using it would add a pilot
  theorem.
- **C: existing cap pilot.** It has the rate but mechanically lacks the extra cap singular value.
  A separate generic pilot would be acceptable only if it independently satisfies the same
  operator rate.

The cap+1 pilot is a theoretical object here; it is not implemented.
