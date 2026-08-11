# Revision-10 Monte Carlo table design

No Monte Carlo was run and no numerical entry was created in RR3.

The main rank-recovery shell is `tables/mc/mc_rank_revision10.tex`. It reports DGP, `N`, `T`, true rank, pilot validity, exact recovery, under-only, over-only, mixed, numerically unresolved, cap-selected count, and modal selected rank. The categories are mutually exclusive among resolved outcomes: exact equals truth; under-only is coordinatewise weakly below with at least one strict inequality; over-only is symmetric; mixed has both an under- and an overselected coordinate.

The appendix shell is `tables/mc/app_mc_ridge_ratio_diagnostics.tex`. It reports blockwise normalized pilot singular values through cap+1, all ridge ratios including rank zero, the smallest and second-smallest ratios and their gap, pilot validity/boundary/objective/stationarity diagnostics, and final-post-refit validity. It contains no candidate, IC, penalty, threshold, or local-completion columns.

The existing DGPs and true-rank vectors `(1,1,1)`, `(2,1,1)`, and `(1,0,2)` remain unchanged. Balanced `N=T` cells remain a simulation design; unequal-`N,T` cells continue to assess rectangular growth. Reporting caps are `(3,3,3)` and spectral-pilot caps are `(4,4,4)`. The latter are never reportable ranks. There is no ridge multiplier or ridge sensitivity experiment.
