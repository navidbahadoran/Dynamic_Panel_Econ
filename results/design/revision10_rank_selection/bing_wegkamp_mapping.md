# Bing-Wegkamp principle and limits

Primary source: Xin Bing and Marten H. Wegkamp, “Adaptive estimation of the rank of the coefficient matrix in high-dimensional multivariate response regression models,” *Annals of Statistics* 47(6), 2019, DOI 10.1214/18-AOS1774 ([author PDF](https://pi.math.cornell.edu/~marten/bw19.pdf); [supplement](https://pi.math.cornell.edu/~marten/supplement_bw19.pdf)).

## Literal construction

Their model is `Y=XA+E` with `Y` n by m, `rank(X)=q`, and one coefficient matrix A of scalar rank r. Baseline theory takes the entries of E iid Gaussian with variance sigma squared; special extensions retain independent/iid errors with sub-Gaussian or finite-fourth-moment restrictions.

Let P project onto the column space of X and `(PY)_k` be its rank-k truncated SVD. For fixed lambda they minimize

`sigma_hat_k^2 = ||Y-(PY)_k||_F^2 / (nm-lambda k)`

over `0<=k<=K_lambda`, where `K_lambda=floor((nm-1)/lambda) min q min m`. The unique minimizer is the last k whose kth squared singular value is at least `lambda sigma_hat_k^2`; the criterion falls to that k and rises thereafter. Thus zero rank is selected when the first step does not improve the criterion.

For STRS, with Z a q by m iid standard Gaussian matrix and `S_j=E[d_j^2(Z)]`, they start with `lambda_hat_0=2(1+epsilon)S_1` and minimize the residual criterion. If the selected rank is positive, they update

`lambda_hat_(t+1)=nm / ((1-epsilon) R_hat_t/U_hat_t + k_hat_t)`,

where `R_hat_t=(n-q)m+sum_{j=2 k_hat_t+1}^{q min m} S_j` and `U_hat_t=max{S_1,S_(2 k_hat_t+1)+S_(2 k_hat_t+2)}`, then minimize only over ranks at least the current rank. Lambda decreases, rank does not decrease, and their high-probability argument keeps every iterate no larger than truth. They stop when the rank is unchanged. Signal recovery uses the separation between a nonzero singular value of XA and the leading projected-noise singular value plus the residual threshold. Their implementation may evaluate the S_j expectations by a separate Gaussian Monte Carlo; that device is part of their model-specific calibration, not evidence for this panel.

## Principle retained

The useful principle is to start below truth, compare genuine signal to a contemporaneous data-scaled noise boundary, relax the boundary only while maintaining no-overfit control, move monotonically upward, and allow rank zero to terminate. The variance level is not inserted as an externally tuned constant.

## Why the formula is not transplanted

- This paper jointly ranks multiple A^(ell), B^(k), and H matrices; there is no single scalar rank or single exchangeable noise scale.
- The dynamic-panel post-refit is a joint, nonlinear, box-constrained low-rank optimization. It lacks the nested truncated-SVD identity and the unique one-dimensional residual-criterion minimizer on which STRS monotonicity rests.
- Lagged outcomes and predetermined regressors make block scores different and data dependent. Conditional temporal/spatial mixing and NED replace iid Gaussian errors; projected-noise singular values do not have the Z calibration used by STRS.
- Rank-vector neighbors need not be globally nested in fitted values after all blocks are jointly refitted. Choosing one coordinate can change later gains, creating an order issue absent from a scalar path.

Consequently neither their residual denominator nor their Gaussian S_j update has a theorem here without a new panel-specific derivation. No official author implementation was needed for, or imported into, this design.
