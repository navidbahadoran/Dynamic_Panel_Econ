# Direct thresholded-pilot rank consistency

## Audit basis and conclusion

The maintained assumptions supplied for this audit are sufficient to imply consistency of the
thresholded pilot rank vector. This is an asymptotic implication; it does not make the observed
locked pilot ranks accurate at N=50 or N=100.

For each coefficient matrix M in the fixed collection
`{A^(1),...,A^(p),B^(1),...,B^(K),H}`, let M0 have fixed rank r_M and let M_tilde be the saved
rank-cap pilot conceptually covered by the maintained pilot lemma. The assumptions used are:

1. `max_M ||M_tilde-M0||_op = o_p(tau_NT)` (the maintained pilot operator-norm error condition);
2. `tau_NT=o(sqrt(NT))` and `tau_NT>0`;
3. for every positive-rank M0, `sigma_{r_M}(M0) >= c_M sqrt(NT)` with fixed `c_M>0` with
   probability tending to one, equivalently the maintained strong-factor singular-value condition;
4. p, K, the number of matrices, and all true ranks are fixed.

## Proof

Fix one matrix M. Weyl's singular-value perturbation inequality gives, for every j,

`|sigma_j(M_tilde)-sigma_j(M0)| <= ||M_tilde-M0||_op`.

Write `e_NT=||M_tilde-M0||_op`. From assumption 1,
`P(e_NT<tau_NT/2)->1`. If `r_M>0`, assumptions 2-3 imply
`sigma_{r_M}(M0)/tau_NT -> infinity` in probability, so
`P(sigma_{r_M}(M0)>2 tau_NT)->1`.

On the intersection of these events, for every `j<=r_M`,

`sigma_j(M_tilde) >= sigma_{r_M}(M0)-e_NT > 3 tau_NT/2 > tau_NT`.

Thus no true positive singular value is thresholded away. For every `j>r_M`,
`sigma_j(M0)=0`, and Weyl gives

`sigma_j(M_tilde) <= e_NT < tau_NT/2 < tau_NT`.

Thus no population-zero singular value exceeds the threshold. If `r_M=0`, only the second
argument is needed and the thresholded rank is zero. Hence the thresholded rank of M_tilde equals
`r_M` with probability tending to one.

The same argument applies separately to every A^(ell), B^(k), and H. Because their number is
fixed, a union bound over the finitely many failure events yields joint recovery of the complete
rank vector with probability tending to one.

## Manuscript implications

The existing pilot operator-error lemma and the manuscript's Weyl perturbation step supply the
entire mathematical argument above. The manuscript itself is not stored in this repository, so
this offline audit cannot responsibly assign exact lemma/proposition numbers or quote their
wording. Before editing, the author should map the maintained pilot lemma and strong-factor
assumption to their Revision-9 labels.

If the pilot became the final rank estimator, the algorithm definition, rank-selection
consistency theorem, and every downstream oracle/inference theorem that currently invokes the
IC-selected rank would need corresponding revision. IC-specific consistency arguments,
`kappa_NT` separation/rate material, candidate post-refit global-gap assumptions used only for IC
selection, and local-completion consistency material could potentially be removed. Candidate
fits might remain computational diagnostics, but they would no longer define the final rank.
No manuscript change is made by this audit.
