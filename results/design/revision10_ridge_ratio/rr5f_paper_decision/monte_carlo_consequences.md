# Monte Carlo consequences

## Route A recommendation

Headline inference Monte Carlo should use supplied true ranks. This isolates the estimator and inference theory actually claimed by the paper.

Recommended design organization for the next separately authorized phase:

- headline cells: `N=T` in `{100,200,400}`;
- `N=T=50`: retain as a labeled small-sample numerical stress design;
- rectangular cells: retain in an appendix and report 100×200 separately from 200×100;
- headline outcomes: estimator bias/RMSE, standard-error accuracy, coverage, rejection, numerical validity, and target-specific diagnostics under supplied true ranks;
- nearby supplied-rank sensitivity: a prespecified, limited set, reported separately from the theorem-aligned true-rank design.

Rank-selection recovery tables should be removed from headline Monte Carlo. RR5/RR5e should not be converted into selector performance tables because no pilot was accepted. At most, a short numerical-development note may state that an investigated automatic pilot was not used. The submitted paper need not carry the full failed-pilot evidence unless transparency or a referee request warrants an appendix diagnostic.

Replace rank-selection tables with supplied-rank sensitivity tables showing how point estimates, uncertainty, fit, boundary activity, and numerical validity change over neighboring inputs. Do not summarize sensitivity as rank-robust inference.

## Other routes

A+ uses the same headline supplied-rank design and may include no positive rank-selection table. B or C would require a completely frozen new protocol and independent validation before any selector Monte Carlo; neither can reuse RR5/RR5e as validation.

No Monte Carlo is authorized by this document.
