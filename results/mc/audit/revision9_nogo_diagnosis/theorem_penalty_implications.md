# Theorem implication of a fixed positive multiplier

The maintained rate definitions are

```text
zeta_NT = b_NT^2 log(NT)^(d_s+2) (1/N + 1/T)
kappa_NT = c_kappa b_NT^2 log(NT)^(d_s+3).
```

Therefore

```text
zeta_NT / {kappa_NT (N+T)/(NT)} = 1 / {c_kappa log(NT)} -> 0
```

for every fixed `c_kappa > 0`. Multiplication by a fixed positive constant also does not alter
the maintained requirement `kappa_NT (N+T)/(NT) -> 0`. With `eta=4`, `d_s=1`, and balanced
`N=T=m`, that term is proportional to `c_kappa m^(-2/3) log(m^2)^6`, which converges to zero.

Thus `c_kappa=1` is a normalization, not a constant uniquely implied by the proof, and any fixed
positive multiplier preserves these asymptotic rate statements. That does not authorize choosing
a multiplier after inspecting production outcomes. If the paper is revised, the multiplier and
its rationale must be fixed before a fresh, independent validation experiment; the revised paper,
configuration, seeds, decision rule, and permitted sensitivities must all be locked in advance.
