# Manuscript assumption cost

## Required change

The minimal paper change is one rank-selection-specific extension to part (i) of
`a:identification`:

> For construction of the spectral rank pilot, part (i) additionally holds with a fixed constant `c_+>0`
> on the fixed pilot class
> `D_max^+`, obtained by increasing each reporting rank cap by one while retaining the same coefficient box.

Define `R_max^+` and `D_max^+` immediately before or after the rank-selection theorem/corollary, and state the
block difference bound `rank(Delta_M)<=2bar_r_M+1`.

The proof of `prop:uniform_concentration` should add that its score and empirical prediction conclusions hold
on `D_max^+` when the enlarged clause is invoked; only constants depending on the fixed enlarged ranks change.
The cap+1 pilot corollary can then cite that version and display the approximate-optimization basic inequality.

## Classification and containment

This is a localized identification-domain extension. It is mathematically stronger than Revision 9, so it
must be disclosed. It is not a new stochastic/dependence/growth assumption and does not change the rate.

The supplied-rank recovery, tangent-space, Riesz, target-expansion, split-correction, and variance theorems
continue to use the original `R_max` identification domain. The enlarged clause is used only to construct and
analyze the spectral pilot. The reported rank set remains `R_max`.

## No other paper cost

There is no change to stability, exogeneity, spatial geometry, NED/mixing, moments, signal strength,
incoherence, coefficient support, loading/factor Gram conditions, target regularity, or rectangular growth.
Neither `B` nor `c_B` is changed. A global joint-regressor lower-eigenvalue primitive may be mentioned as a
sufficient condition, but it need not replace the minimal restricted-domain formulation.
