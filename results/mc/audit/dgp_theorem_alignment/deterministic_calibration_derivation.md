# Deterministic Monte Carlo calibration derivation

## Population variance of the primitive disturbance

The heteroskedastic scale satisfies `sigma_i^2 ~ U[0.5,1.5]`, so
`E[sigma_i^2]=1`.  DGP 1 uses `u_tilde_it=sigma_i epsilon_it`, with
`Var(epsilon_it)=1`, and therefore

`Var_population(u_tilde_it)=1`.

For DGPs 2--4, `z_0t=epsilon_0t` and

`z_it=rho_s z_(i-1,t)+sqrt(1-rho_s^2) epsilon_it`.

Because the initial value has variance one, induction gives `Var(z_it)=1` for
every unit.  Independence of `sigma_i` and `z_it` again gives
`Var_population(u_tilde_it)=E[sigma_i^2] Var(z_it)=1`.  Spatial covariance
changes dependence, not this marginal variance.  DGP 4 groups do not enter the
disturbance construction, so there is no group adjustment.

## Population variance of H_raw

The bounded time factor follows

`g_j=0.5 g_(j-1)+sqrt(1-0.5^2) z_j`, `g_-1=0`, `Var(z_j)=1`.

Thus `Var(g_j)=1-0.5^(2(j+1))`.  With burn-in `b=50`, the average marginal
variance across the T observed dates for rank one is

`V_H,1(T) = 1 - (1/T) sum_(j=b+1)^(b+T) 0.5^(2j)`.

Since the loading is independent, mean zero, and variance one,
`Var(lambda_i g_t)=V_H,1(T)`.  The correction from one is below displayed
double precision for every maintained T, so the numerical table reports 1.

For rank-two H stress, the added component is the product of independent
variance-one bounded loading and factor draws.  Its maintained strength is one.
The deterministic support-preserving scale is

`s_H=(3 sqrt(3))/(3 sqrt(3)+3)=0.6339745962155613`.

Hence

`V_H,2(T)=s_H^2 [V_H,1(T)+1]`

and its maintained numerical value is `0.803847577293368`.

## Deterministic c_h

For H rank `r_H`, the replacement is

`c_h(r_H,T)=sqrt([pi_H/(1-pi_H)] * 1/V_H,r_H(T))`, `pi_H=0.30`.

The maintained numerical values are:

- H rank one: `c_h=0.6546536707079772`;
- H rank two: `c_h=0.7301712917987002`.

These values depend only on prespecified primitive distributions, burn-in, T,
rank, and deterministic stress strengths.  At displayed precision they are
common across N, T, and DGP.  They do not use any realized calibration or
production disturbance.

The population H share is exactly

`c_h^2 V_H / (c_h^2 V_H + 1) = 0.30`.

The independent finite calibration samples have realized shares from
`0.2917529981` to `0.3108818271`; those fluctuations are diagnostics and do not
feed back into `c_h`.

## Frozen c_xi

Closed-form population pooled R-squared is impractical here because the dynamic
outcome includes random heterogeneous A and B matrices, predetermined x in
DGPs 3--4, burn-in, and spatial dependence.  The candidate therefore uses the
permitted alternative: a separate calibration experiment with seed 20260807
and 50 draws for each DGP x N x T x rank design.  The outcome recursion is
affine in `c_xi`; its pooled centered sum of squares is evaluated exactly as a
quadratic, and the positive root targets 0.65.  The resulting constants are
then frozen in a TOML table before any production draw.

For positive slope rank, frozen `c_xi` ranges from `0.8231633625` to
`2.4725313610`.  Calibration-sample achieved R-squared equals 0.65 to numerical
root tolerance.  For `(1,0,2)`, scale is unidentified, `c_xi=1`, and induced
R-squared is retained rather than forced to 0.65.

Cell-specific constants are proposed because the maintained estimand is a
finite-panel pooled R-squared and the candidates vary materially across N, T,
and rank-design calibration samples.  A DGP-level common constant would be
simpler, but it would deliberately give different finite-cell R-squared values
and would require a separate design decision.  The full cell table is explicit
and fixed ex ante, so cell specificity does not introduce replication-level
randomness.

## Coefficient envelopes

The prior analytical results remain:

- `C_A=0.85` in every design;
- `C_beta=2.031384387633061` in DGPs 1--3 with positive slope rank;
- `C_beta=1.959615242270663` in DGP 4 with positive slope rank;
- `C_beta=0` for `(1,0,2)`.

For each frozen row,

`C_H=3 sqrt(3) |c_h c_xi|`,

and `C_Theta=max(C_A,C_beta,C_H)`.  The common maximum is
`C_Theta_max=8.410761115894578`.

Therefore `B=9,c_B=1` does not satisfy `C_Theta_max<=8`.  The smallest simple
integer proposal under `c_B=1` is `B=10`, for which the deterministic interior
margin is `10-C_Theta_max=1.589238884105422`.  This is a reported candidate,
not an activated estimator/configuration change.

