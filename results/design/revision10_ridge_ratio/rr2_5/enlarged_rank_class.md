# Enlarged cap+1 pilot class

Let the fixed coefficient-block collection be

`J={A^(1),...,A^(P),B^(1),...,B^(K),H}`.

For each block `M`, retain the reporting cap `bar_r_M` and the true-rank restriction
`0<=r_M,0<=bar_r_M`. The reported-rank set remains

`R_max=product_M {0,1,...,bar_r_M}`.

Only for the spectral pilot, define

`R_max^+=product_M {0,1,...,bar_r_M+1}`

and

`D_max^+={Theta-Theta_0: Theta in union_(r in R_max^+) M(r), ||Theta||_max<=B}`.

Equivalently, the pilot parameter space is

`M_cap+1={Theta: rank(M)<=bar_r_M+1 for every M, ||Theta||_max<=B}`.

Truth is feasible because `r_M,0<=bar_r_M<bar_r_M+1` and
`||Theta_0||_max<=B-c_B<B`.

For `Delta=Theta-Theta_0` in `D_max^+`, rank subadditivity gives, block by block,

`rank(Delta_M)
 <=rank(Theta_M)+rank(M_0)
 <=bar_r_M+1+r_M,0
 <=2bar_r_M+1`.

Write

`s_M^+=bar_r_M+1+r_M,0` and `S_+=sum_(M in J) s_M^+`.

The number of blocks, every `bar_r_M`, every `s_M^+`, and `S_+` are fixed constants. Also
`||Delta_M||_max<=2B` and `||Delta||_F<=2B sqrt(|J|NT)`. The extra rank changes neither the reported set nor
the coefficient box. It exists only to expose the `(bar_r_M+1)`st pilot singular value.
