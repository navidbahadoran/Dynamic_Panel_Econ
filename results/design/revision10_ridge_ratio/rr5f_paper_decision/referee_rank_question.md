# Referee question: “How are the ranks chosen in practice?”

## Strongest honest Route A response

“The block ranks are specification inputs, not estimated nuisance parameters in our inferential theorem. Before examining target estimates, we define a small economically and computationally plausible rank grid for the separately ranked coefficient matrices. We report the headline supplied-rank specification together with neighboring-rank estimates, fit diagnostics, and singular-value diagnostics. Conclusions are described as stable only when they remain so across this prespecified sensitivity set. Our contribution is valid target-specific inference conditional on supplied ranks; we do not claim rank-robust inference or a consistent automatic rank selector.”

This answer is sufficient if the paper's title, abstract, introduction, theorems, Monte Carlo, and empirical tables consistently present target-specific inference—not factor-number estimation—as the contribution. The sensitivity grid must be prespecified and small enough to be interpretable, not selected after seeing favorable target estimates.

## Strongest honest Route A+ response

The Route A answer must remain the operative answer. One may add: “An appendix gives a conditional asymptotic rank-consistency result for a spectral pilot satisfying a stated operator rate. The currently implemented nonconvex pilot did not meet our numerical acceptance criteria, so it is not used for empirical rank choice.” This is honest but invites the further question of why the auxiliary theorem belongs in the paper.

## Routes B and C

Neither provides a current practical answer. Under B, ranks would depend on a newly changed numerical rule requiring independent validation. Under C, ranks would depend on a new pilot whose statistical definition and theory do not yet exist. Until those phases are complete, neither route may be described as operational.
