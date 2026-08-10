# Revision-9 thresholded rank-pilot diagnosis

This is an offline audit of 24 already-computed locked selected-rank records. It ran no DGP,
screening, fitting, IC selection, split fitting, Riesz solve, inference, or Monte Carlo.

## Cap-pilot ranks and numerical status

The pooled thresholded-rank distribution is: (1,0,1): 1, (1,1,2): 1, (1,1,3): 2, (1,2,3): 5, (1,3,3): 5, (2,1,2): 1, (2,1,3): 2, (2,2,2): 1, (2,2,3): 4, (2,3,3): 1, (3,1,3): 1. Exact recovery is 0/24; Over only
is 23/24; Under only is 1/24;
Mixed is 0/24. At N=50 the counts are 0 Exact, 1 Under only, 11 Over
only, 0 Mixed. At N=100 they are 0 Exact, 0 Under only, 12 Over only, 0 Mixed. Every DGP has zero
exact recovery.

All 24 pilots are recorded as converged, objective-stable, and stationarity/KKT-valid under the
locked interior (`1e-6`) or constrained (`1e-4`) criterion. 13/24 are boundary-active by
the saved coefficient-envelope diagnostic. This is therefore not a recorded pilot
numerical-validity failure; the accepted rank-adaptive pilot outputs are mostly too large.

## Candidate coverage source

Truth is in the candidates 24/24. Direct cap-pilot attribution is 0/24. A nuclear-path
proposal equals truth in 24/24. Truth is also a one-coordinate neighbor of the cap
pilot in 2/24 and of at least one nuclear proposal in 18/24. The
coverage attribution is nuclear-direct only in 6/24 and multiple-source in
18/24. The complete per-replication mechanisms are in
`candidate_coverage_source.csv`. The direct pilot did not account for any of the observed 24/24
coverage.

## Singular-value threshold margins

The locked rank rows save `tau_NT` and the thresholded vector but do not save the accepted pilot's
matrix-specific singular values. The fit-diagnostic table retains only aggregate `sigma_1` and
`sigma_r` for individual route fits; it cannot identify sigma_2 for each accepted A, B, and H.
Therefore sigma_1/tau, sigma_2/tau, signal/noise margins, and near-threshold counts are unavailable
without rerunning a prohibited fit. `cap_pilot_singular_value_margins.csv` records this limitation
for all 72 matrix-realization pairs rather than fabricating values.

## Pilot versus final IC

The pilot is incorrect and the final IC is incorrect in 24/24. The other three requested cells
(correct pilot/incorrect IC, incorrect pilot/correct IC, correct pilot/correct IC) are all zero.
Thus the IC did not destroy correct information already present in the pilot: the pilot overfit
while the locked IC selected `(0,0,0)` in every realization.

## Nuclear path

At least one path proposal equals `(1,1,1)` in 24/24. This makes the nuclear path an effective
screening source in these records, not a demonstrated final estimator. Modal proposals and full
saved proposal sequences are reported in `nuclear_path_rank_summary.csv`. Pooled over all 528
path positions, the modal rank is ["(1,1,3)"] with frequency
144.

## Theory and evidence classification

The maintained `o_p(tau_NT)` pilot operator-error condition, strong-factor singular values of
order `sqrt(NT)`, `tau_NT=o(sqrt(NT))`, and fixed matrix count directly imply joint thresholded
rank consistency by Weyl's inequality. The complete argument is in
`direct_pilot_rank_theory.md`.

The evidence classification is **CASE P3**: thresholded cap-pilot recovery is poor (0/24 exact),
so removing the IC would not solve the locked finite-sample problem. This does not contradict the
asymptotic proof; it shows that its separation regime is not visible in these two saved sizes.
Missing singular spectra prevent a margin diagnosis but do not obscure the observed 0/24 rank
recovery.

## Paper-level options

1. Retaining the current IC preserves Revision 9 but must acknowledge its finite-sample NO-GO.
2. Using pilot ranks has a direct asymptotic proof, but the locked evidence rejects it as the
   current finite-sample solution; it should not be implemented on this record.
3. A new rank selector requires new theory and a separately locked validation design. It is the
   only option aimed at resolving both observed failures, but is a future paper-level project.

Recommended next decision: retain the NO-GO while the author decides whether to keep Option 1 or
undertake Option 3. Do not adopt Option 2 from these data and do not search for another fixed
`c_kappa`.

## Separate N=50 supplied-rank issue

The existing locked audit reports four narrow full-fit stationarity failures, nine boundary-active
split fits, broad-target inference failure at N=50, and clean N=100 supplied-rank performance.
That supports reporting N=50 as a small-sample stress design rather than a headline fixed-rank
design, without deleting it or changing B, tolerances, or the estimator. This issue is distinct
from selected-rank failure, which occurs at both N=50 and N=100.
