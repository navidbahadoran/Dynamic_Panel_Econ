# Primary and supporting literature mapping

## Pu et al. (2025)

Primary reference: Dan Pu, Kuangnan Fang, Wei Lan, Jihai Yu, and Qingzhao Zhang,
“Reduced Rank Spatio-Temporal Models,” *Journal of Business & Economic Statistics* 43(1),
98–109, DOI [10.1080/07350015.2024.2326142](https://doi.org/10.1080/07350015.2024.2326142)
([author-hosted PDF](https://kuangnanfang.com/zb_users/upload/2024/04/202404161713254298791730.pdf)).

- **Low-rank object.** Their spatial dynamic model is
  `Y_t=B_1Y_t+B_2Y_(t-1)+epsilon_t`, with `B_1=U_1V` and `B_2=U_2V`; both square
  coefficient matrices share one d-dimensional weight/loading matrix V.
- **Spectral object.** For autocovariances `Gamma_k=cov(Y_(t+k),Y_t)`, they use
  `M=sum_(k=1)^k0 Gamma_k Gamma_k'`, whose rank is d, and its sample counterpart from sample
  autocovariances.
- **Criterion.** If `lambda_hat_i` are decreasing eigenvalues of the sample M, they select
  `argmin_(1<=i<=dmax) (lambda_hat_(i+1)+c)/(lambda_hat_i+c)`. Their method does not include
  rank zero.
- **Ridge condition.** Condition C7 requires `c=o(1)` and `n^2/(cT)=o(1)`. Their practical
  recommendation `(n log n)^2/(20T)` contains a simulation-motivated numerical factor and is not
  imported here.
- **Signal and theorem.** Condition C6 requires the smallest nonzero singular value of their
  transition object `S_n20` to be bounded below. Under C1–C7, Theorem 2 proves
  `P(d_hat=d)->1` as n,T diverge with `n=o(T^(1/2))`.
- **Disturbances.** C3 assumes the `epsilon_it` are independent, mean zero, heteroskedastic across
  i within fixed variance bounds, and have moments through order six. Their spatial dependence is
  primarily in the coefficient dynamics, not the paper's conditionally spatially dependent
  innovations.
- **Computation.** QMLE is optimized by gradient descent with Armijo line search and multiple
  randomized initial values, retaining the maximum objective. Rank selection computes sample
  autocovariances, their spectral matrix, and the ratios.

The publisher describes the supplement as containing derivative formulae, auxiliary lemmas,
proofs of Theorems 1–2, the optimization algorithm, covariate extensions, and further simulations.
The supplement payload was protected by the publisher's browser challenge in this environment;
the theorem, conditions, criterion, and algorithm summarized above are stated explicitly in the
primary article and publisher record. No claim below depends on an inaccessible supplementary
step.

## Barigozzi–He–Li–Trapani and related spectral work

- [Robust Tensor Factor Analysis](https://arxiv.org/abs/2303.18163) studies a Tucker tensor factor
  model and a modified mode-wise eigenvalue ratio. Its ridge is tied to an upper rate for noise
  eigenvalues; the authors prove factor-number consistency without restrictions on relative
  divergence rates and permit heavy tails. Their iterative implementation deliberately retains
  extra eigenvectors to avoid spectral truncation.
- [Statistical Inference for Large-dimensional Tensor Factor Model by Iterative
  Projections](https://arxiv.org/abs/2206.09800) allows weak cross-sectional and temporal
  correlation in tensor idiosyncratic components and uses iterative projected covariance
  eigenvalue ratios.
- [Barigozzi and Trapani's nonstationary factor-dimension procedure](https://arxiv.org/abs/1806.03647)
  derives eigenvalue-divergence tests with no relative N/T rate restriction under weak factor-model
  dependence.

These papers support the spectral-separation principle, finite extra-spectrum requirement, and
possibility of rectangular asymptotics. They concern covariance eigenvalues of observed factor or
tensor processes, not singular values of jointly estimated dynamic-panel coefficient matrices.
Their assumptions and ridge formulae are therefore not transplanted. The RR1 proof uses only the
Revision-9 pilot operator rate and signal bounds.
