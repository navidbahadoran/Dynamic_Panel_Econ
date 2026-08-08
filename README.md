# Low-Rank Heterogeneous Dynamic Panels: Monte Carlo Replication

This repository implements the Monte Carlo design and estimator for a dynamic panel with
separately low-rank autoregressive coefficients, covariate slopes, and interactive effects.
It is designed as referee-facing research code: configurations are immutable at run time,
randomness is reproducible across worker counts, failures remain in the output, and reported
estimates always come from unpenalized supplied-rank least squares.

## 1. Scope

The repository currently implements the Monte Carlo experiment. It does not implement or
compare an alternative estimator, and it contains no empirical data or fabricated empirical
results. The supplied older archive was used only as a guide to practical command-line and
output organization; none of its cross-fitting, buffering, fold, or one-step econometric logic
is used here.

## 2. Econometric model

For the baseline design,

\[
y_{it}=a_{it}y_{i,t-1}+\beta_{it}x_{it}+H_{0,it}+u_{it}.
\]

The matrices `A`, `B`, and `H` have their own ranks and their own factors and loadings. The
internal `Coefficients` representation also supports lists of multiple lags and covariates.
No factors or loadings are shared across coefficient matrices.

## 3. Method-to-code map

| Responsibility | Location |
|---|---|
| Coefficient algebra, fitted operator, exact adjoint | `src/dynamic_panel_econ/core.py` |
| DGPs 1--4, burn-in, DGP 4 truths and gap pilot | `src/dynamic_panel_econ/dgp.py` |
| Deterministic `c_H`/`c_xi` calibration | `src/dynamic_panel_econ/calibration.py` |
| Joint supplied-rank ALS and nuclear screening | `src/dynamic_panel_econ/estimation.py` |
| Screened candidates, post-refits, IC, local completion | `src/dynamic_panel_econ/rank_selection.py` |
| Linear target directions and exact truths | `src/dynamic_panel_econ/targets.py` |
| Matrix-free Riesz and two-way split correction | `src/dynamic_panel_econ/inference.py` |
| Bartlett spatial HAC | `src/dynamic_panel_econ/spatial.py` |
| Seeds, chunks, multiprocessing, manifests | `src/dynamic_panel_econ/monte_carlo.py` |
| Aggregation and paper-facing tables | `src/dynamic_panel_econ/reporting.py` |

## 4. Repository structure

```text
configs/mc/                  smoke, pilot, production, and rank-stress TOML files
scripts/                     validation, simulation, aggregation, and table CLIs
src/dynamic_panel_econ/      tested production package
tests/                       mathematical identity and integration tests
results/mc/                  trackable generated calibration/results/tables
```

Manuscripts are not part of the repository. Generated table fragments under `results/mc` are
deliberately not ignored.

## 5. System requirements

- Python 3.11 or newer (the package is also tested here with Python 3.13)
- NumPy, SciPy, Pandas, PyArrow, and threadpoolctl
- Windows PowerShell, macOS, or Linux
- Enough memory for several `N x T` arrays per worker; production `400 x 400` jobs should be
  piloted before selecting the worker count

The nuclear path, repeated joint post-refits, matrix-free Riesz iterations, and four split fits
for every broad target are the principal costs.

## 6. Installation

PowerShell:

```powershell
cd D:\Programming\Dynamic_Panel_Econ
py -3.11 -m venv .venv
.\.venv\Scripts\Activate.ps1
python -m pip install --upgrade pip
pip install -e ".[dev]"
```

POSIX shells:

```bash
cd /path/to/Dynamic_Panel_Econ
python3.11 -m venv .venv
. .venv/bin/activate
python -m pip install --upgrade pip
pip install -e '.[dev]'
```

## 7. Quick start

```powershell
pytest -q
ruff check src scripts tests
python scripts\validate_dgp.py --config configs\mc\smoke.toml
python scripts\report_dgp4_truths.py --config configs\mc\production.toml --draws 100
python scripts\run_mc.py --config configs\mc\smoke.toml --n-jobs 1 --overwrite
python scripts\aggregate_mc.py --config configs\mc\smoke.toml
python scripts\make_mc_tables.py --config configs\mc\smoke.toml
python scripts\run_mc.py --config configs\mc\rank_stress_smoke.toml --overwrite
python scripts\aggregate_mc.py --config configs\mc\rank_stress_smoke.toml
python scripts\make_mc_tables.py --config configs\mc\rank_stress_smoke.toml
```

The smoke run uses two `20 x 20` replications per DGP. It exercises rank selection, a local
fixed-time target, a broad corrected target, diagonal/spatial variance, aggregation, and tables.
It is a software check, not evidence about statistical performance.

## 8. Tests and code quality

`pytest -q` covers reproducible DGP arrays, burn-in and lag alignment, spatial and predetermined
dependence, exact DGP 4 truths, true ranks, the fitted-map adjoint identity, joint ALS behavior,
nuclear screening at `lambda_max`, target reconstruction, tangent projections, empirical Riesz
equations, and dense-versus-lag-sum Bartlett calculations. Run `ruff check src scripts tests`
before sharing changes.

## 9. Monte Carlo DGPs

All designs generate response periods `t=-49,...,T` from `y[i,-50]=0` and discard the first 50.
The returned lag at observed time `t` is the actually simulated `y[i,t-1]`.

- DGP 1 uses independent Gaussian heteroskedastic disturbances.
- DGP 2 adds the `rho_s=0.5` one-dimensional spatial recursion in O(`NT`) work.
- DGP 3 adds `eta_x * u_tilde[i,t-1]` to the covariate. It never uses the current disturbance.
- DGP 4 retains DGP 3 and uses deterministic equal groups with separate `A` and `B` loading
  means. It remains rank one.

The common factors use bounded standardized-uniform innovations. The raw autoregressive matrix
is multiplied by one common
`c_a=min(1,0.85/max(abs(A_raw)))`; no entrywise clipping is used. Every replication stores the
pre-scaling maximum, `c_a`, and the final maximum.

The simulation coefficient box is fixed at `B=9` with deterministic interior margin
`c_B_sim=1`. The bounded-AR support calculation is applied to every calibrated DGP/cell before
tasks are dispatched. Higher-rank stress blocks are multiplied by one deterministic common
factor per matrix so their support envelope equals the corresponding rank-one envelope; they are
never clipped entrywise. Manifests and rows retain theoretical and realized blockwise envelopes.

## 10. DGP 4 group truths and the gap pilot

Every DGP 4 row and the standalone truth-report command store realized fixed-time and time-average
G1 and G2 means and G2-minus-G1 contrasts for both `A` and `B`, after all scaling. They also store the raw pre-stability `A`
contrasts and `c_a`.

`auto_adjust_group_gap=false` by default. The separate deterministic pilot evaluates only the
prespecified gap grid. If enabled, it widens the two `A` loading means symmetrically around their
original center, chooses the smallest gap meeting the configured floor in every design cell,
and writes the frozen means to `resolved_config.json` and `run_manifest.json` before production
draws. It never uses estimates, bias, coverage, power, or production replications.

The validation command prints the exact post-scaling DGP 4 group means and differences and writes
`dgp4_group_truths.csv` plus `dgp4_group_gap_pilot.json`.

## 11. Calibration and an identified feasibility issue

For each DGP/cell, deterministic calibration draws are separate from production seeds. The code
uses the prescribed

\[
c_H^2=\frac{\pi_H}{1-\pi_H}
      \frac{\operatorname{Var}(\widetilde u)}{\operatorname{Var}(H^{raw})},
\qquad \pi_H=0.30,
\]

then brackets `c_xi>0` and calls `scipy.optimize.brentq` for the requested average pooled R-squared.
All trial values reuse common random numbers. Conditional on a raw draw, the response is affine in
`c_xi`, so calibration precomputes its two dynamic-recursion components. Output records the minimum
grid R-squared, approximate large-`c_xi` floor, feasibility, and successful bracket.

There is an important reproducible feasibility problem in the supplied statistical specification:
with `H0=c_xi*c_H*H_raw` and `u=c_xi*u_tilde`, the lagged-outcome and persistent interactive-effect
terms can keep pooled R-squared above 0.50 even as `c_xi` becomes arbitrarily large. In deterministic
DGP 1 checks its numerical floor is about 0.58. Therefore the mandated 0.50 root need not exist.
The code does not invent a replacement scaling or report a closest value: it raises a transparent
calibration failure. The production configuration preserves 0.50; the smoke configuration uses
the explicit feasible value 0.70 solely to test the full software pipeline. Before production,
the econometric specification must resolve whether a different object is intended to be scaled or
a feasible R-squared target should be prespecified.

Calibration output includes targets, achieved moments, coefficient ranges, stability summaries,
and primitive conditional variance summaries.

Rank-stress calibration is separate for every DGP, cell, and true-rank vector and uses the actual
rescaled stress matrices. A further structural feasibility issue is recorded, not overridden:
when the true vector is `(1,0,2)`, the slope block is zero and both the interactive effect and
primitive disturbance scale with `c_xi`. The resulting outcome is homogeneous in `c_xi`, so its
pooled R-squared cannot be tuned by `c_xi`; at the current cells it is about 0.52--0.54 rather than
0.65. The strict runner therefore refuses the common-target rank experiment until the statistical
specification supplies a compatible resolution.

## 12. Targets and exact truth

The mathematical indices are `i0=floor(N/4)` and `t0=floor(T/2)`, converted once from one-based
notation to zero-based Python indices. Every target is a `Coefficients` direction `D`; truths and
plug-ins are always computed as `<D,Theta>`.

Implemented targets include entries, fixed-time cross-sectional means, DGP 4 fixed-time group
means and contrasts, full-panel means, and time-averaged group means and contrasts. Entries and
fixed-time targets use the full-panel plug-in. Full-panel and time-averaged targets use the
two-way split correction.

## 13. Estimator pipeline

1. Solve the convex nuclear-norm path from `lambda_max` downward. The nuclear-plus-box prox uses
   inner Dykstra iterations; one SVT followed by clipping is not treated as exact.
2. Fit one joint rank-at-most-cap pilot. Starting from nuclear-path singular spaces, a
   rank-adaptive outer routine tests one-coordinate normal-gradient increases and numerically
   redundant decreases, jointly refits after every proposed move, and stops when no admissible
   move improves the objective beyond tolerance.
3. Threshold path and cap-pilot singular values at `sqrt(NT)/log(NT)`.
4. Add valid one-coordinate neighbors only.
5. Jointly refit every candidate by unpenalized supplied-rank ALS from two deterministic starts;
   use a third deterministic randomized start only when needed. A candidate is numerically
   unresolved unless at least two valid starts have normalized objective gap at most
   `start_objective_stability_tol`; an unresolved candidate cannot minimize the IC.
6. Minimize the specified IC and locally complete omitted one-coordinate neighbors.
7. Use the selected unpenalized full-panel fit for target estimation and Riesz inference.

Nuclear estimates are screening and warm starts only. They are never reported estimates.

## 14. Rank selection and stress design

The dimension is the sum of `r*(N+T-r)` over matrices. Revision 8 uses
`b_NT=(NT)^(1/(8+eta))*log(NT)` and
`kappa_NT=b_NT^2*log(NT)^(spatial_dimension+3)`. The lattice designs set
`spatial_dimension=1`, so the logarithmic exponent is four. The penalty fixes
`eta_for_penalty=4.0` in configuration and never estimates it from a simulated sample. Ties prefer
smaller dimension and then lexicographic rank order. Raw records retain candidate counts, IC gaps,
cap hits, candidate sources, convergence, and selected ranks.

`configs/mc/rank_stress.toml` prespecifies `(1,1,1)`, `(2,1,1)`, and `(1,0,2)`. The generic stress
generator adds independent bounded loading-factor components, makes rank-zero blocks exactly zero,
applies deterministic common-factor envelope rescaling, and then applies the same autoregressive
stability scaling. Calibration is performed separately on each actual true-rank design.

Rank diagnostics are stored once per replication under `rank/`, not repeated on every target row.
They include candidate coverage of the actual true vector, exact/under/over/zero-rank recovery,
cap margins, candidate validity, and nuclear-path proposals. Prespecified stability checks cover
the dense `sqrt(0.8)` nuclear grid, threshold multipliers `0.5,1,2`, IC multipliers `0.5,1,2`,
larger caps, and best one-coordinate target changes. They are diagnostics, not rank-robust inference.

## 15. Riesz computation

For every selected matrix, the tangent projector is

\[
P_T(Z)=UU'Z+ZVV'-UU'ZVV'.
\]

The ambient matrix-free normal operator applies `P_T A* A P_T`. Conjugate gradients are primary,
with MINRES and LSMR fallbacks. Convergence is judged against the original empirical Riesz
equation. The target-specific diagnostic is named
`riesz_target_rayleigh_quotient`; its numerator and denominator both use the projected solution
returned by the solver. It is not described as a minimum eigenvalue.

A separate matrix-free Lanczos diagnostic uses exactly `r(N+T-r)` nonredundant orthonormal tangent
coordinates per block. It reports estimated smallest and largest tangent-Gram eigenvalues and a
condition-number estimate without introducing ambient normal-space zero eigenvalues.

## 16. Two-way split-panel correction

Ranks are selected once on the full panel. Independent balanced time and unit partitions then
produce exactly four more fits with those ranks fixed. Unit splits are stratified within target
groups. Time splits retain observed lagged outcomes even when the lag time belongs to the other
half. Restricted directions add exactly to the original target.

The estimate is `3*phi_full - sum(phi_time_halves) - sum(phi_unit_halves)`. The corrected score uses
three distinct fitted objects at every observation:

```text
3 * Psi_full * residual_full
- Psi_time * residual_time
- Psi_unit * residual_unit
```

There are no buffers, folds, rank reselection, cross-fitting, or same-panel one-step correction.
For every positive supplied split rank, the numerical record also contains `sigma_1`, `sigma_r`,
and `sigma_r/sigma_1`. A prespecified numerical relative-rank floor flags computational collapse;
this is separate from the statistical threshold and never selects ranks.

## 17. Variance and spatial HAC

DGP 1 uses the diagonal score sum. DGPs 2--4 use Bartlett spatial weights with
`h_N=ceil(c_sp*log(NT))`. The O(`T*N*h_N`) lag-sum implementation never constructs a dense matrix
in production. There is no temporal HAC.

## 18. Parallel execution

Set `parallel_level` to `replications` or `none`. The production runner parallelizes replications
and keeps candidate fitting serial inside workers. `threadpoolctl` caps
BLAS threads to avoid oversubscription. Do not nest process pools.

```powershell
python scripts\run_mc.py --config configs\mc\pilot.toml --n-jobs 4
```

## 19. Random seeds

SHA-256 converts semantic identifiers into stable 32-bit seed components passed to NumPy
`SeedSequence`. Keys include the master seed, DGP, cell, replication, and stage. Python's randomized
`hash()` and worker order never enter a persistent seed.

## 20. Resume and restart

Each run lives under `results/mc/<name>/<config-hash>/`. Chunks have deterministic replication
ranges, write first to `.partial`, and are atomically renamed. `--resume` skips completed chunks;
`--overwrite` explicitly recomputes them. A changed statistical or numerical setting changes the
hash. Workers return records to the parent and never append concurrently to one file.

```powershell
python scripts\run_mc.py --config configs\mc\pilot.toml --resume --n-jobs 4
```

## 21. Outputs and schemas

- `resolved_config.json`: complete immutable settings and hash
- `group_gap_pilot.json`: prespecified candidates and frozen DGP 4 means
- `calibration.parquet`: cell calibration and achieved moments
- `raw/*.parquet`: one immutable replication chunk per deterministic range
- `rank/*.parquet`: compact replication-level rank-selection diagnostics
- `summary/mc_summary.{parquet,csv}`: target-level performance statistics
- `summary/rank_summary.{parquet,csv}`: true-rank-conditioned selection statistics
- `tables/tab_mc_*.{tex,csv,parquet}`: deterministic paper-facing fragments
- `run_manifest.json`: code identifier, failure vocabulary, requested count, frozen group means

Truth may vary by replication. Bias, RMSE, tests, and coverage therefore use each row's exact truth.

## 22. Reproducing tables

```powershell
python scripts\run_mc.py --config configs\mc\pilot.toml --n-jobs 4
python scripts\aggregate_mc.py --config configs\mc\pilot.toml
python scripts\make_mc_tables.py --config configs\mc\pilot.toml
```

Substitute `production.toml` only after resolving the 0.50 calibration feasibility issue and after
timing the pilot. Builders create `tab_mc_dgp1` through `tab_mc_dgp4`, `tab_mc_main_summary`,
`tab_mc_coeff_summary`, `tab_mc_bias_correction`, `tab_mc_dgp4_groups`, `tab_mc_rank`,
`tab_mc_computation`, and `tab_mc_spatial_sensitivity` in matching LaTeX, CSV, and Parquet formats.
Performance tables use named mathematical target panels rather than an internal target-string column.

## 23. Failure rules and diagnostics

No requested replication disappears. It produces successful target rows or a failure row with a
standardized status such as calibration failure, invalid cap pilot, no stable pair of start
objectives, no valid candidate,
nonconvergence, high stationarity residual, active coefficient bound, cap hit, Riesz failure,
split rank loss, split target-support loss, split Riesz target instability, nonpositive variance, or an
unexpected exception with its type and message. Aggregation reports requested and successful
counts plus failure counts. A selected rank at its cap never produces a primary interval, and an
invalid candidate cannot minimize the IC. Raw/rank records include objectives, iteration counts, ranks, IC gaps,
Riesz diagnostics, split assignments/counts, cutoff, stage runtimes, seeds, config hash, and commit.

## 24. Performance notes

Start with `smoke.toml`, then time `pilot.toml`. Reduce worker count when memory pressure or BLAS
oversubscription appears; do not silently reduce replications. Profile nuclear SVDs, joint ALS,
Riesz products, and split fits before changing algorithms. Candidate caps and numerical tolerances
are explicit configuration choices.

## 25. Empirical reuse/status

The operator, supplied-rank estimator, target, Riesz, and spatial modules are designed for later
unemployment and housing applications. No empirical application or data construction is claimed
in this repository at present.

## 26. Git and manuscript safety

`.gitignore` excludes common manuscript directories, root `ver*.tex`/`ver*.pdf`, revision drafts,
and private audit notes. It does not globally ignore TeX, PDF, data, results, tables, or figures.
The table scripts never read, edit, or copy a manuscript.

## 27. Citation and license

See `LICENSE`. Until the final article citation is available, cite this repository as:

```text
Low-Rank Heterogeneous Dynamic Panels: Monte Carlo Replication Package,
version 0.1.0, forthcoming article citation.
```

## 28. Troubleshooting

- **Calibration root not bracketed:** inspect the requested R-squared and the feasibility note in
  Section 11. Do not silently accept a boundary value.
- **ALS nonconvergence:** inspect objective history, stationarity residual, rank support, and bound
  ratio; increase configured sweeps only as a prespecified numerical change.
- **Riesz failure:** inspect the original-equation residual, target-specific Rayleigh quotient,
  and the separately reported tangent-coordinate Gram spectrum. Small split panels can lack target
  support even when the full panel is identified.
- **Memory pressure:** lower `n_jobs` and keep one BLAS thread per worker.
- **Interrupted run:** rerun with `--resume`; incomplete `.partial` files are never treated as
  completed chunks.
- **Output already exists:** choose `--resume` or the explicit destructive recomputation option
  `--overwrite`.
