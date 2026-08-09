# Initial-condition audit

## Code behavior

The generator initializes all recursive state before the first generated burn-in observation as follows:

- `y_i,-50 = 0` through `previous = zeros(N)` in `_simulate`.
- `g_a,-50 = g_b,-50 = g_h,-50 = 0` because every call to `_ar_uniform` starts from `previous = 0`.
- `f_x,-50 = 0` in the recursion used for the covariate common factor (the code variable is the zero-mean AR factor itself).
- `x_i,-50 = 0` through `previous_x = zeros(N)`.

The first stored AR factor is generated from that zero initial state. The code generates one additional primitive disturbance at `t=-50`; DGPs 3-4 use it as `tilde_u_i,t-1` when constructing the first covariate value after initialization.

## Missing from the paper

The paper states the initial values for `y`, `g_a`, and `g_b`, but it must also state:

1. `g_h,-50 = 0`;
2. `f_x,-50 = 0`;
3. `x_i,-50 = 0`;
4. the indexing rule that a primitive `tilde_u_i,-50` is drawn and enters the first DGP-3/4 covariate recursion;
5. whether the 50-period burn-in means the first retained observation follows 50 generated transition steps (the code does).

These are paper underspecifications, not code mismatches.
