# Final rectangular-rate check

Let `m=NT`, `h_NT=1/N+1/T`, and retain

`b_NT=m^(1/(8+eta))log(m)`,

`zeta_NT=m^(2/(8+eta))log(m)^(d_s+4)h_NT`.

The maintained growth condition is

`G_NT=m^(4/(8+eta))log(m)^(d_s+6)h_NT->0`.

With the frozen ridge `a_NT=1/log(m)`, `N,T->infinity` implies `a_NT->0`, and exactly

`zeta_NT/a_NT
 =m^(2/(8+eta))log(m)^(d_s+5)h_NT
 =G_NT/[m^(2/(8+eta))log(m)]->0`.

The cap+1 class changes only constants depending on fixed ranks. It changes neither `zeta_NT` nor this
identity. The theorem retains `1/N+1/T` throughout and imposes no `N=T`, `N/T->c`, `N asymp T`, or hidden
polynomial ordering. Balanced Monte Carlo cells, if later used, are computational design choices rather than
theory restrictions.
