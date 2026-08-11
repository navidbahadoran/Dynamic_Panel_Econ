# RR1 decision

Classification: **RR1-PASS**.

The ridge-ratio selector is consistent under the existing rectangular-panel assumptions, conditional on the
same cap-pilot objective-accuracy condition already used by Revision 9. The proof is internal to this paper:

- the enlarged fixed cap+1 pilot has operator error `O_p(sqrt(NT zeta_NT))`;
- Weyl gives order-one normalized positive singular values and `O_p(zeta_NT)` post-rank values;
- the maintained growth condition implies `zeta_NT log(NT)->0` for every admitted rectangular sequence;
- `a_NT=1/log(NT)` therefore separates signal, ridge, and noise without a selected constant;
- the anchor `lambda_M,0=1` resolves rank zero;
- the extra pilot singular value resolves truth at the reporting cap;
- fixed blocks and caps give joint consistency by a finite union.

There are zero new substantive DGP assumptions. Pu et al.'s independence and `n=o(T^(1/2))` restrictions
are not used. The cap+1 pilot and its existing-style numerical objective condition are localized
rank-selection changes, not new DGP restrictions.

RR1-PASS authorizes theorem/design work only. It does not authorize implementation, simulation, or a
manuscript edit.
