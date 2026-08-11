# Rectangular ridge-rate derivation

Write `m=NT` and `h_NT=1/N+1/T`. The manuscript defines

`b_NT=m^(1/(8+eta)) log(m)`

and

`zeta_NT=b_NT^2 log(m)^(d_s+2) h_NT
         =m^(2/(8+eta)) log(m)^(d_s+4) h_NT`.

Its maintained growth condition is

`G_NT=m^(4/(8+eta)) log(m)^(d_s+6) h_NT -> 0`.

For `a_NT=1/log(m)`, clearly `a_NT->0` because both N and T diverge. Moreover,

`zeta_NT/a_NT
 =m^(2/(8+eta)) log(m)^(d_s+5) h_NT
 =G_NT/[m^(2/(8+eta)) log(m)]`.

The numerator tends to zero and the denominator tends to infinity. Hence

`zeta_NT/a_NT -> 0`.

This algebra uses neither `N=T`, a finite N/T limit, nor a polynomial ordering between N and T.
It holds for every rectangular sequence admitted by the existing panel-growth assumption. It also
shows `zeta_NT->0`, as required by the pilot perturbation argument.
