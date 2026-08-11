# Route B: modify the numerical acceptance rule

## What each relaxation would address

### Objective-agreement relaxation

Only 5/180 RR5e pilots had two individually valid starts; none of those five met the frozen objective-agreement tolerance. Relaxing objective agreement could affect at most those five pilots. It does not address the 150 pilots with zero valid starts or the 25 with one. It therefore cannot rescue the design's dominant failure by itself.

### One-valid-start acceptance

Allowing one valid start could make at most 30/180 RR5e pilots eligible before other diagnostics: 25 with one valid start and 5 with two. It leaves 150/180 with no valid start. It also removes the main empirical check that the nonconvex fit is objective-stable across starts. This is a paper-level procedural change, not an engineering detail.

### Stationarity/KKT relaxation

This is the only listed change that directly targets the dominant failure: 505/540 starts failed the frozen stationarity/KKT condition, and 150/180 pilots had zero valid starts. Relaxing these criteria could increase valid-start counts. It would also redefine what the paper accepts as a solution to the nonconvex cap+1 problem and weaken the numerical evidence supporting the theorem's approximate-global-optimization premise.

## Scientific assessment

All three changes would be chosen after observing two complete 0/180 pilot-acceptance experiments. Without an independent numerical principle, certified error analysis, or theorem connecting a revised diagnostic to the required objective gap, they carry severe post-hoc-tuning risk. No tolerance value is justified here.

Any revision would require a new paper-level numerical protocol fixed before data generation, an explanation of why it remains credible for the theoretical objective-gap premise, deterministic solver validation, and fresh independent Monte Carlo seeds. Previous RR5/RR5e evidence could diagnose the old rule but could not validate the new one.

## Recommendation

Not recommended. Objective and one-start relaxations do not address the dominant failure, while stationarity/KKT relaxation would change the paper's numerical standard precisely because the frozen procedure failed.
