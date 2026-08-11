# Route A: supplied-rank main paper

## Assessment

Route A makes the supplied-rank estimator the only operative estimator in the theory, Monte Carlo, and empirical application. This is the most direct response to RR5 and RR5e: the supplied-rank recovery and inference theory remains intact, while the implemented cap+1 least-squares pilot failed the maintained acceptance gate in 360/360 fresh attempts. The evidence concerns computation of that pilot, not statistical performance of ridge ratios.

The theoretical cost is low and mostly subtractive. The supplied-rank theorems retain their original identification domain and assumptions. The cap+1 pilot rate, enlarged pilot class, ridge ratios, rank-zero anchor, and joint rank-consistency theorem are no longer needed for the operative method. Removing them reduces the number of high-level numerical premises.

The manuscript cost is moderate because Revision 10 currently gives rank selection a visible role. The benefit is a cleaner match between theorem, computation, and evidence. Referee risk shifts from “the advertised selector cannot be computed” to the manageable question “where do ranks come from?” The honest answer is that ranks are specification inputs, prespecified over a small set, and substantive conclusions are checked over neighboring rank vectors. The paper must not call this rank-robust inference.

Empirical practicality is good when the block ranks are small: report estimates over a prespecified grid, use fit and singular-value diagnostics descriptively, and identify conclusions stable across neighboring specifications. This is familiar specification-sensitivity practice, although it cannot be marketed as automatic consistent rank selection.

## Ridge-ratio theorem disposition

Remove the ridge-ratio theorem and its cap+1 pilot from the submitted paper rather than retaining a theorem for an operative procedure the project cannot compute reliably. Preserve the theorem, proof, implementation, and RR5/RR5e evidence in the repository as research history. A later paper or revision could reconsider it only after an independently justified computational pilot exists.

## Costs and risks

- Theory cost: low; deletion and theorem renumbering.
- Manuscript cost: moderate; contribution and simulation narrative must be rewritten.
- Referee concern: ranks are inputs, so the paper must show a disciplined ex-ante sensitivity protocol.
- Empirical clarity: high if the supplied-rank grid is prespecified and compact.
- Submission timing: shortest credible route.

## Recommendation

Recommended. It aligns the paper's claims with feasible computation, protects the intact inference contribution, avoids post-hoc numerical tuning, and removes unnecessary pilot-only notation.
