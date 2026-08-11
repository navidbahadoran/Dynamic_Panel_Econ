# Rectangular ridge-rate audit

Let `m=NT`, `h_NT=1/N+1/T`, and use the manuscript definitions

`b_NT=m^(1/(8+eta)) log(m)`,

`zeta_NT=b_NT^2 log(m)^(d_s+2) h_NT
        =m^(2/(8+eta)) log(m)^(d_s+4) h_NT`.

The maintained growth condition is

`G_NT=m^(4/(8+eta)) log(m)^(d_s+6) h_NT ->0`.

With `a_NT=1/log(m)`, both dimensions diverge under the maintained panel asymptotics, so `m->infinity`
and `a_NT->0`. Moreover,

`zeta_NT/a_NT
 =m^(2/(8+eta)) log(m)^(d_s+5) h_NT
 =G_NT/[m^(2/(8+eta)) log(m)] ->0`.

This is an exact identity. It never replaces `1/N+1/T` by one of its terms and never invokes a ratio between
`N` and `T`. It also implies `zeta_NT->0`.

## Unequal-growth examples admitted by the same condition

Let integer parts be understood.

1. `T=N^(1+eta/8)`. The polynomial exponent in `G_NT` is
   `4(2+eta/8)/(8+eta)-1=-eta/[2(8+eta)]<0`; logarithms are dominated.
2. `N=T^(1+eta/8)`. This is the transposed unequal sequence and has the same negative exponent.
3. `T=N{log N}^10`. The polynomial factor is `N^(-eta/(8+eta))`; every fixed logarithmic power is
   dominated, so the condition holds.

The ridge proof applies unchanged to all three and, more generally, every rectangle satisfying `G_NT->0`.

## Hidden-balance checklist

- Pilot error retains `sqrt(NT zeta_NT)` exactly.
- Spectra are divided by `NT`, not `N^2` or `T^2`.
- The fitted-value scale weights use the actual `NT` observations.
- The ratio proof uses only `zeta_NT=o(a_NT)`.
- The finite union uses fixed blocks and caps, not dimension balance.
- Downstream transfer is an event argument and adds no rate.

No line requires `N=T`, `N/T->c`, or either polynomial ordering prohibited by the task.
