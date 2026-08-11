# Complete cap+1 pilot operator-rate proof

Assume the localized identification extension and the enlarged uniform concentration result. Let

`Theta_hat^pil in M_cap+1`

be feasible and let

`delta_NT=L_NT(Theta_hat^pil)-inf_(Theta in M_cap+1)L_NT(Theta)=o_p(zeta_NT)`.

Set `Delta_hat=Theta_hat^pil-Theta_0`, `n=NT`,

`Gamma_hat(Delta,Delta)=sum_it Y_it(Delta)^2`,

and `S(Delta)=sum_it u_it Y_it(Delta)`.

Truth belongs to `M_cap+1`, so

`inf_(Theta in M_cap+1)L_NT(Theta)<=L_NT(Theta_0)`

and therefore

`L_NT(Theta_hat^pil)-L_NT(Theta_0)<=delta_NT`.

Because `y_it=Y_it(Theta_0)+u_it`, the complete squared-loss expansion is

`L_NT(Theta_hat^pil)-L_NT(Theta_0)
 =(1/(2n)) sum_it[{u_it-Y_it(Delta_hat)}^2-u_it^2]
 =(1/(2n))Gamma_hat(Delta_hat,Delta_hat)-(1/n)S(Delta_hat)`.

Consequently,

`(1/2)Gamma_hat(Delta_hat,Delta_hat)
 <=S(Delta_hat)+n delta_NT`.

The enlarged concentration bounds give, jointly,

`|S(Delta_hat)|<=C_n sqrt(n zeta_NT)||Delta_hat||_F`, with `C_n=O_p(1)`,

and

`Gamma_hat(Delta_hat,Delta_hat)
 >=(c_+/2)||Delta_hat||_F^2-C n zeta_NT`.

Substitution gives

`(c_+/4)||Delta_hat||_F^2
 <=C_n sqrt(n zeta_NT)||Delta_hat||_F
   +(C/2)n zeta_NT+n delta_NT`.

Since `delta_NT=o_p(zeta_NT)`, the last term is `o_p(n zeta_NT)`. With probability tending to one it can be
absorbed into a random `O_p(n zeta_NT)` term. For `x=||Delta_hat||_F`, `A_n=O_p(sqrt(n zeta_NT))`, and
`D_n=O_p(n zeta_NT)`, the inequality is

`(c_+/4)x^2<=A_n x+D_n`.

The positive root bound

`x<=4A_n/c_++2sqrt(D_n/c_+)`

therefore yields

`||Theta_hat^pil-Theta_0||_F=O_p(sqrt(NT zeta_NT))`.

Finally, for every block,

`||M_hat^pil-M_0||_op
 <=||M_hat^pil-M_0||_F
 <=||Theta_hat^pil-Theta_0||_F`,

so the fixed block collection satisfies

`max_M ||M_hat^pil-M_0||_op=O_p(sqrt(NT zeta_NT))`.

This proves the generic RR2 pilot requirement. The proof retains the approximate-global-objective term and
does not infer globality from observable numerical diagnostics.
