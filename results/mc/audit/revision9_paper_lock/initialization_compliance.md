# Monte Carlo initialization compliance

| State at period `-50` | Revision-9 intended value | Current code value | Result |
|---|---:|---:|---:|
| `y_i,-50` | 0 | 0 (`previous=zeros(N)` in outcome recursion) | MATCH |
| `x_i,-50` | 0 | 0 (`previous_x=zeros(N)`) | MATCH |
| `g_a,-50` | 0 | 0 (`_ar_uniform` starts from `previous=0`) | MATCH |
| `g_b,-50` | 0 | 0 | MATCH |
| `g_h,-50` | 0 | 0 | MATCH |
| `f_x,-50` | 0 | 0 (the zero-mean AR covariate factor uses `_ar_uniform`) | MATCH |

**Overall: MATCH.** The code additionally draws `tilde_u_i,-50` so DGPs 3--4 can use the correct lagged disturbance in the first covariate transition after initialization.
