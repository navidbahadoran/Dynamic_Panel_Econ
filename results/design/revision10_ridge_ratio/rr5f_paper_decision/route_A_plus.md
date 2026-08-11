# Route A+: supplied-rank main method plus auxiliary ridge-ratio theorem

## Assessment

Route A+ is mathematically coherent. Supplied-rank estimation and inference would be the only headline method, while an appendix theorem would state joint ridge-ratio consistency conditional on a cap+1 pilot satisfying the proved operator rate and normalized global objective-gap condition. The paper would disclose that the current nonconvex constrained least-squares implementation was not used for headline inference because RR5 and RR5e produced no acceptable pilot.

This route is potentially publishable if the auxiliary theorem is sharply quarantined: no finite-sample success claim, no empirical rank choice based on it, and no implication that its numerical premise was verified. Mathematically, the theorem remains valid; computationally, its current pilot is inaccessible under the maintained diagnostic standard.

The difficulty is intellectual economy. The theorem requires the pilot-only enlarged class, Assumption 8 extension, `R_max^+`, `D_max^+`, the cap+1 objective-gap premise, spectral normalization, ridge, anchor, and ratio proof. Retaining all this asks readers and referees to absorb a second method that does not produce the paper's reported estimates. A predictable referee question is why a computationally unusable procedure occupies theorem space in an inference paper.

## Costs and risks

- Theory cost: already paid mathematically, but exposition and qualification remain substantial.
- New assumptions: the localized pilot identification extension remains.
- Implementation risk: high for the auxiliary procedure; irrelevant to the operative estimator only if separation is explicit.
- Referee risk: medium-high because the computational premise lacks supporting finite-sample evidence.
- Empirical clarity: lower than Route A; readers may confuse auxiliary rank consistency with implemented rank handling.
- Submission timing: longer than A because the manuscript must maintain two carefully separated tracks.

## Recommendation

Not recommended for the present submission. The theorem is valid but distracts from the feasible contribution and retains notation whose only operative purpose has failed numerically. Preserving it in the repository is preferable to presenting it as an auxiliary result without a credible implemented pilot.
