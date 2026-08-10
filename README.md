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

The nuclear path, repeated joint post-refits, matrix-free Riesz iterations, and exactly four
split coefficient fits per replication are the principal costs. Split coefficient optimization
is not repeated when the target direction changes.

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

### Fixed-rank and selected-rank modes

The two primary methods share DGP seeds, targets, Riesz inference, split correction, variance
estimation, and failure rules. Fixed-rank mode bypasses every rank-selection operation:

```powershell
python scripts\run_mc.py --dgp-grid 1,2,3,4 --balanced-grid 50,100,200,400 `
  --replications 1000 --rank-mode fixed --fixed-ranks 1,1,1 `
  --pooled-r2-target 0.65 --workers 8 --output-root results/mc/production_fixed --resume
```

Selected-rank mode runs the locked Revision-9 selection pipeline. The paper fixes the baseline
IC multiplier at exactly `1.0`:

```powershell
python scripts\run_mc.py --dgp-grid 1,2,3,4 --balanced-grid 50,100,200,400 `
  --replications 1000 --rank-mode selected --rank-caps 3,3,3 `
  --kappa-f-b 0.15 --coefficient-bound 10 --ic-multiplier 1.0 `
  --pooled-r2-target 0.65 --frozen-calibration configs/mc/frozen_dgp_calibration.toml --workers 8 `
  --output-root results/mc/production_selected --resume
```

The canonical non-executing paper-validation command is:

```powershell
python scripts\run_mc.py --config configs\mc\production.toml --pooled-r2-target 0.65 `
  --kappa-f-b 0.15 --coefficient-bound 10 --rank-caps 3,3,3 `
  --ic-multiplier 1.0 --print-resolved-config --dry-run
```

Run `python scripts/run_mc.py --help` for the complete interface. Configuration precedence is
CLI over TOML over package defaults. `--dry-run --print-resolved-config` resolves and displays a
design without calibration, fitting, inference, or output. Every executed run writes
`command.txt`, `resolved_config.toml`, and `manifest.json`/`run_manifest.json` with the Git commit
and CLI arguments.

`--pooled-r2-target` is solely a DGP calibration parameter. It is not estimator tuning, rank
selection tuning, or an estimated-model fit target. The official designs read `c_h` and `c_xi`
from `configs/mc/frozen_dgp_calibration.toml`; they never recalibrate from a production draw.
The independently calibrated, frozen identified-cell target is 0.65. In rank-stress cells with
`r_B=0`, `c_xi=1` remains the normalization, the requested R-squared is not imposed, and the
induced R-squared is reported.

### Lossless Monte Carlo accounting and reporting

Every attempted replication receives one mutually exclusive `primary_status`. Point performance
uses every finite numerically valid estimate; inference performance additionally requires a
finite positive standard error and variance. No estimate is trimmed or winsorized. The optional
median/MAD `extreme_estimate_flag` is descriptive and never changes bias, RMSE, coverage, power,
or rejection inclusion.

The canonical reporting command is:

```powershell
python scripts\report_mc.py --input-run results/mc/production_fixed/RUN `
  --input-run results/mc/production_selected/RUN `
  --output-dir results/mc/method_comparison --tables all --figures all
```

Individual artifacts use `--table optimization`, `--table failure_accounting`,
`--figure rmse`, `--figure coverage`, `--figure runtime`, or `--figure power_A`. Reporting refuses
to overwrite an existing artifact unless `--overwrite` is explicit. Failure and retained-share
denominators always use attempted replications; coverage and rejection use valid-inference
replications. Runtime summaries retain all completed runtimes and report mean, median, p10, p90,
and p95.

The smoke run uses one `50 x 50` replication per DGP, the smallest configured cell for which the
prespecified 0.65 calibration is feasible in all four designs. It exercises rank selection,
entry and fixed-time targets, theorem-applicability tagging, several broad corrected targets,
diagonal/spatial variance, aggregation, and tables. It is a software check, not evidence about
statistical performance.

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

The simulation and estimator coefficient box is the fixed parameter-space constant `B=10`, with
`c_B=1`. It is not estimated or changed by replication. The common analytical DGP envelope is
`C_Theta,max=8.288745227963506`, so its distance to the boundary B is
`1.711254772036494`, while its additional slack relative to the required envelope `B-c_B=9` is
`0.711254772036494`. Higher-rank stress blocks are multiplied by one deterministic common factor
per matrix so their support envelope equals the corresponding rank-one envelope; truths are never
clipped entrywise. Manifests and rows retain theoretical and realized blockwise envelopes.

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

The official DGP uses an ex-ante frozen table. Its `c_h` is calculated analytically from population
moments using

\[
c_H^2=\frac{\pi_H}{1-\pi_H}
      \frac{\operatorname{Var}(\widetilde u)}{\operatorname{Var}(H^{raw})},
\qquad \pi_H=0.30,
\]

For all four DGPs, `Var_population(u_tilde)=1`. Rank-one H has population variance one and
`c_h=0.6546536707079772`; rank-two H stress has population variance `0.803847577293368` and
`c_h=0.7301712917987002`. Thus population `pi_H=0.30` without using a realized disturbance
variance. A separate 50-draw calibration-only experiment brackets `c_xi>0` and calls
`scipy.optimize.brentq`; its constants are frozen by DGP, N, T, and rank design before production.
Startup validates the requested cells, analytical `c_h`, coefficient envelopes, and table hash.

Locked Revision 9 prespecifies `pooled_r2_target=0.65`. Smoke, pilot, Riesz-diagnostic, and
production configurations use this value; calibration draws are independent of reported
replications and the active frozen table is never recalculated by a production replication.

Calibration output includes targets, achieved moments, coefficient ranges, stability summaries,
and primitive conditional variance summaries.

Rank-stress calibration is separate for every DGP, cell, and true-rank vector and uses the actual
rescaled stress matrices. A further structural feasibility issue is recorded, not overridden:
when the true vector is `(1,0,2)`, the slope block is zero and both the interactive effect and
primitive disturbance scale with `c_xi`. The outcome is homogeneous in `c_xi`. These cells set
`c_xi=1` as an explicit normalization, record `r2_scale_identified=false`, set requested
R-squared to null, and report the induced pooled R-squared. They are not calibration failures.

## 12. Targets and exact truth

The mathematical indices are `i0=floor(N/4)` and `t0=floor(T/2)`, converted once from one-based
notation to zero-based Python indices. Every target is a `Coefficients` direction `D`; truths and
plug-ins are always computed as `<D,Theta>`.

Implemented targets include entries, fixed-time cross-sectional means, DGP 4 fixed-time group
means and contrasts, full-panel means, and time-averaged group means and contrasts. Entries and
fixed-time targets use the full-panel plug-in. Full-panel and time-averaged targets use the
two-way split correction.

Every target row stores the true ratio `||P0 D||/||D||`; entry rows also store scaled unit and
time leverage. Fixed-time G2-minus-G1 contrasts in DGPs 1--3 are retained only as
`weak_target_stress_outside_assumption9` and excluded from headline theorem-coverage tables. The
same contrasts remain theorem-covered in DGP 4. `B_entry` is theorem-covered under the locked
factor `f_b=0.6+0.15*g_b`, whose deterministic support is `[0.15,1.05]`.

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

The dimension is the sum of `r*(N+T-r)` over matrices. Locked Revision 9 uses
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
the dense `sqrt(0.8)` nuclear grid, threshold sensitivities `0.5,2`, IC sensitivities `0.5,2`,
larger caps, and best one-coordinate target changes. They are diagnostics, not rank-robust inference.
Each sensitivity uses the same rank-at-most-cap pilot, stable candidate post-refits, and local
candidate-completion algorithm, changing only its named tuning input.

## 15. Riesz computation

For every selected matrix, the tangent projector is

\[
P_T(Z)=UU'Z+ZVV'-UU'ZVV'.
\]

Each fitted sample prepares one reusable nonredundant tangent-coordinate Riesz system. Its
matrix-free Gram operator and optional Lanczos spectrum are cached across targets. Conjugate
gradients are primary, with MINRES and LSMR fallbacks. Only the projected target right-hand side
and solution are target-specific. The diagnostic is named
`riesz_target_rayleigh_quotient`; its numerator and denominator both use the projected solution
returned by the solver. It is not described as a minimum eigenvalue.

A separate matrix-free Lanczos diagnostic uses exactly `r(N+T-r)` nonredundant orthonormal tangent
coordinates per block. It reports estimated smallest and largest tangent-Gram eigenvalues and a
condition-number estimate without introducing ambient normal-space zero eigenvalues.
If the target tangent norm is below `target_support_tolerance`, it is classified
`target_unsupported_selected_rank` before any Riesz iteration. When Gram diagnostics are enabled,
eigensolver failure or an estimated minimum eigenvalue below
`tangent_gram_min_eigenvalue_floor` suppresses the primary interval under its own status.

## 16. Two-way split-panel correction

Ranks are selected once on the full panel. Independent balanced time and unit partitions then
produce exactly four more coefficient fits with those ranks fixed. These fits, residuals,
partitions, rank diagnostics, and Riesz systems are prepared once per replication and reused for
every broad target. Unit splits are stratified within deterministic groups. Time splits retain
observed lagged outcomes even when the lag time belongs to the other half. Restricted directions
add exactly to the original target.

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

The numerical fixed-rank routine solves the paper's literal problem over the supplied-rank
factorization and the reconstructed entrywise box. It first runs unconstrained ALS. A solution
strictly inside `B` is the fast-path constrained solution. Otherwise, alternating row and time
convex quadratic subproblems impose `-B <= loading*factor' <= B` directly. Successful boundary
activity is retained as a valid constrained point estimate and remains in bias, RMSE, and
Monte Carlo-SD summaries. It receives `boundary_interiority_failure` for primary theorem-based
inference, so it is excluded from coverage and rejection denominators. Other failures distinguish
constrained solver, feasibility, optimality/KKT, and nonfinite solutions.

The maintained slope factor has `mu_f_b=0.6` and `kappa_f_b=0.15`. Since the bounded AR support is
`g_b in [-3,3]`, the factor satisfies `f_b in [0.15,1.05]`; this gives the deterministic positive
time-leverage floor used for theorem-covered `B_entry` inference.

## 18. Parallel execution

Set `parallel_level` to `replications` or `none`. The production runner parallelizes replications
as one global outer pool over `DGP x N x T x semantic replication`; it keeps the nuclear path,
cap pilot, candidate post-refits, local completion, and split fits serial inside each worker.
There is no inner process pool. On Windows the pool uses explicit `spawn` semantics and only
module-level, pickle-safe worker functions.

Use the single `n_jobs` TOML value or the `--n-jobs` CLI option. Precedence is CLI, then TOML,
then the conservative default of one. The resolved configuration and manifest record both
`requested_n_jobs` and `effective_n_jobs`; the latter is capped at the number of available outer
tasks. A request is never silently replaced by a hidden worker count.

```powershell
python scripts\run_mc.py --config configs\mc\pilot.toml --n-jobs 4
python scripts\run_mc.py --config configs\mc\production.toml --n-jobs 12
```

The correct value is machine dependent. `n_jobs=1` is the deterministic debugging setting. When
more than one worker is active, `threadpoolctl` explicitly limits detected BLAS, LAPACK, and
OpenMP libraries to the configured `blas_threads` value (normally one), preventing multiplication
such as 12 Python workers times 8 native threads. Configuration and frozen calibrations are
initialized once per worker. Panel/design arrays are generated and owned by one worker rather
than cached and copied from the parent.

On Windows, logical processor count is available in Task Manager under **Performance > CPU** or
with `Get-CimInstance Win32_ComputerSystem | Select-Object NumberOfLogicalProcessors`. Task Manager
also reports total CPU utilization and process/system memory under **Performance** and **Details**.
Monitor both during a computational smoke: increasing workers is useful only while wall time falls
and memory retains safe headroom.

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
- `summary/target_regularity.{parquet,csv}`: true projection and entry-leverage diagnostics
- `tables/tab_mc_*.{tex,csv,parquet}`: deterministic paper-facing fragments
- `run_manifest.json`: code identifier, failure vocabulary, requested count, frozen group means

Truth may vary by replication. Bias, RMSE, tests, and coverage therefore use each row's exact truth.

## 22. Reproducing tables

```powershell
python scripts\run_mc.py --config configs\mc\pilot.toml --n-jobs 4
python scripts\aggregate_mc.py --config configs\mc\pilot.toml
python scripts\make_mc_tables.py --config configs\mc\pilot.toml
```

Substitute `production.toml` only after one fresh locked-Revision-9 preflight passes and the
worker count is fixed ex ante. Builders create `tab_mc_dgp1` through `tab_mc_dgp4`, `tab_mc_main_summary`,
`tab_mc_coeff_summary`, `tab_mc_bias_correction`, `tab_mc_dgp4_groups`, `tab_mc_rank`,
`tab_mc_computation`, `tab_mc_target_regularity`, and `tab_mc_spatial_sensitivity` in matching LaTeX, CSV, and Parquet formats.
Performance tables use named mathematical target panels rather than an internal target-string column.

## 23. Failure rules and diagnostics

No requested replication disappears. It produces successful target rows or a failure row with a
standardized status such as calibration failure, invalid cap pilot, no stable pair of start
objectives, no valid candidate,
nonconvergence, constrained solver/feasibility/KKT failure, cap hit, Riesz failure,
split rank loss, unsupported selected-rank targets, full/split tangent-Gram numerical failures,
split Riesz target instability, nonpositive variance, or an
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
