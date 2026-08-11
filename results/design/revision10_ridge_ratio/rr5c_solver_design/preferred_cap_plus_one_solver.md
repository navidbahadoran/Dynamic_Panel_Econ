# Preferred semantically exact cap+1 solver

## Decision-ready architecture

Use a deterministic **gauge-aware, safeguarded alternating solver with QP-specific constrained updates**. It retains the three maintained starts, the `(cap+1)` factor widths, `B=10`, normalized least-squares objective, update order, and all frozen acceptance thresholds. It changes no coefficient-product feasible set and adds no penalty or rank threshold.

### 1. Maintained start construction

Construct the same three seeded Gaussian starts in the same order and with the same RNG calls. Record their initial factors and products before any numerical transformation. Start identity is therefore unchanged.

### 2. Interior path

1. Perform the existing loading-row then time-row alternating least-squares sweep.
2. Solve every joint least-squares system by rank-revealing SVD/QR with the frozen `rcond=1e-10`; use deterministic signs and pivot tie rules. Never form normal equations.
3. After each *completed* sweep, reconstruct each product and apply the balanced thin-SVD gauge
   \[
   M=P\Sigma Q',\quad U=P\Sigma^{1/2},\quad V=Q\Sigma^{1/2},
   \]
   zero-padded to width `cap+1`. This is used only on the interior path, whose frozen stationarity diagnostic is product based.
4. If a zero singular direction has both factor columns zero, preserve the zero product while placing a deterministic orthonormal complement on one side and zero on the other. Alternate the chosen side by a fixed documented rule. This keeps the full width available; it is not rank selection and uses no threshold beyond exact/numerical zero handling required by the factorization routine.
5. Apply a monotone safeguard: retain a snapshot of the factors and objective before every half/full sweep; commit the candidate only if its recomputed frozen objective is finite and nonincreasing within the existing roundoff allowance. Otherwise restore both factors and objective from the snapshot. This repairs the current latent mismatch in which an increasing interior sweep can leave post-sweep factors paired with the prior stored objective. A rejected accelerated step falls back to the unaccelerated accepted ALS iterate, not to a different estimator.
6. When ordinary ALS stalls above the frozen product projected-gradient target, take deterministic quotient-aware Gauss-Newton directions in the current fixed-width factorization. Remove only gauge-null components from the *direction*, not product directions. Use a deterministic Armijo/backtracking line search on the same objective and require the product to remain strictly interior. If no admissible descent step exists, retain the ALS iterate.
7. Stop and validate with the unchanged objective rule and the unchanged `1e-6` product tangent projected-gradient calculation. If the product is not strictly inside `B-1e-8`, use the constrained path.

The Gauss-Newton refinement is an optimizer for the same factorized least-squares objective, not a new estimator. Its role is to address the 325/339 interior starts that ended above the frozen stationarity threshold despite finite objectives.

### 3. Constrained path

1. Start from the same uniformly scaled interior result and the same balanced square-root SVD refactorization used now.
2. Preserve the loading-row, blockwise QR, time-row update order and full widths.
3. For each row/time problem, form the identical convex quadratic objective and identical product-box linear constraints.
4. Apply an invertible, deterministic column equilibration only inside the subproblem: `x=D z`, transform `(X,C)` to `(X D,C D)`, solve for `z`, then map back to the maintained factor coordinate `x`. Do not retain the transformed gauge globally.
5. Replace generic SLSQP with a deterministic dense primal-dual active-set QP solver specialized to the quadratic/linear-inequality problem. Warm-start both the primal point and active set from the previous sweep. Use QR/SVD null-space solves and deterministic constraint ordering/ties. Continue until the existing subproblem feasibility and KKT acceptance rules pass; do not change their tolerances or maximum-iteration rule without a separate engineering amendment.
6. Cache the block-diagonal constraint matrix and symbolic work for all loading rows within a half-sweep, and analogously for time rows. Cache only values mathematically invariant within that half-sweep.
7. Commit an update only if it passes the existing finite, feasibility, nonincrease, and subproblem KKT rule. Preserve the prior feasible coordinate if a trial fails; record the failure exactly.
8. At each full sweep, compute the same maximum factor-space box KKT residual in the maintained factor coordinates. Stop only under the existing objective and `1e-4` KKT rules.

Internal equilibration improves the linear algebra but maps the solution back before the frozen KKT evaluator. Arbitrary non-orthogonal global gauge balancing is intentionally excluded from the constrained baseline because the finite normalized factor-KKT residual is coordinate dependent, even though exact KKT zero is invariant.

### 4. Multistart and outputs

Run all three maintained starts. Apply the unchanged per-start validity logic. Require at least two valid starts and the unchanged normalized best-two objective gap `<=1e-6`. Compute exact full singular values only after pilot acceptance. Preserve `rank_selection_numerically_unresolved` with no fallback rank if the gate fails.

### 5. Determinism, resources, Windows safety

- Stable semantic start IDs and existing seeds; no new stochastic choices.
- Fixed unit/time/constraint order and deterministic tie-breaking.
- One native BLAS/OpenMP thread per worker and no nested pool.
- Per-worker reusable workspaces bounded by the small joint factor width; no `NT x NT` objects.
- Spawn-safe module-level worker functions and immutable worker-local configuration.
- Telemetry is observational and excluded from scientific outputs/hashes where timestamps differ.

## Why this is preferred

It targets both observed failure modes: product-based interior stationarity and expensive constrained coordinate optimality. It uses product-preserving gauge control where the frozen residual is gauge invariant, while avoiding a hidden change to the constrained finite residual. It naturally permits lower-rank products through zero singular values, retains every cap+1 direction, and keeps the literal product box exact. The main implementation cost is the deterministic QP active-set kernel and its tests; memory remains modest because subproblems have only the sum of block widths as variables.
