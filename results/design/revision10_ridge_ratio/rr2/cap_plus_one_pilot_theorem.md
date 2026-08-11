# Cap-plus-one pilot audit and unresolved implementation corollary

## Generic pilot condition

The rank theorem is algorithm-free. Its statistical pilot requirement is

`max_M ||M_hat^pil-M_0||_op=O_p(sqrt(NT zeta_NT))`,

conditionally on the common-shock sigma-field. No rank theorem premise refers to a nuclear path, starting
value, or optimization algorithm.

## Proposed joint cap+1 pilot

Let

`M_cap+1={Theta: ||Theta||_max<=B, rank(M)<=bar_r_M+1 for every M}`

and let the computed joint pilot be feasible and satisfy

`delta_NT=L_NT(Theta_hat^pil)-inf_{Theta in M_cap+1} L_NT(Theta)=o_p(zeta_NT)`.

Truth is feasible because `r_M,0<=bar_r_M`. For `Delta=Theta-Theta_0`, every block difference can have

`rank(Delta_M)<=bar_r_M+1+r_M,0<=2bar_r_M+1`,

a fixed bound. Fixedness is enough for the covering and score arguments, but it is **not** enough for the
maintained prediction-identification assumption.

## The precise gap

Assumption `a:identification`, lines 940--953 of the supplied manuscript, states its global lower bound for
differences `Delta=Theta-Theta_0` with `Theta` ranging over ranks in `R_max`. When truth is at a reporting cap,
this controls differences of rank at most `2bar_r_M`. The proposed enlarged pilot permits a candidate of rank
`bar_r_M+1`, whose difference from cap-rank truth can have rank `2bar_r_M+1`. Such a candidate is outside the
stated identification domain.

This is not repaired merely because all ranks are fixed. A linear fitted-value map can have a null direction
of rank `2bar_r_M+1` while being bounded below on the compact set of unit-norm reporting-cap differences.
Choose cap-rank truth `M_0` and a cap+1-rank `M_1` so that `D=M_1-M_0` has rank `2bar_r_M+1`, and let the fitted-
value map have kernel spanned by `D` but no nonzero intersection with the lower-rank reporting difference
class. Then the maintained restricted identification can hold, while `M_0` and `M_1` have identical loss.
An exact global cap+1 optimizer may consequently contain an order-`sqrt(NT)` extra singular value. An
`o_p(zeta_NT)` objective gap does not prevent it.

## Conditional basic-inequality derivation

Because truth is feasible,

`L_NT(Theta_hat^pil)<=L_NT(Theta_0)+delta_NT`.

Writing `Delta_hat=Theta_hat^pil-Theta_0` and expanding squared loss gives

`(2NT)^(-1)||Y(Delta_hat)||_2^2
 <=(NT)^(-1)<u,Y(Delta_hat)>+delta_NT`.

The existing uniform score bound extends to this fixed-rank class. **If**, in addition, the prediction lower
bound were valid on the enlarged cap+1 difference class, it would give

`<u,Y(Delta)> = O_p(sqrt(NT zeta_NT)) ||Delta||_F`,

and existing conditional prediction identification plus empirical concentration gives

`||Y(Delta)||_2^2 >= c||Delta||_F^2-O_p(NT zeta_NT)`

uniformly. Multiplying the basic inequality by `NT`, absorbing fixed constants, and using
`NT delta_NT=o_p(NT zeta_NT)` yields

`c||Delta_hat||_F^2
 <=O_p(sqrt(NT zeta_NT))||Delta_hat||_F+O_p(NT zeta_NT)`.

The quadratic inequality implies

`||Delta_hat||_F=O_p(sqrt(NT zeta_NT))`.

Hence, for the fixed block collection,

`max_M ||M_hat^pil-M_0||_op
 <=||Delta_hat||_F
 =O_p(sqrt(NT zeta_NT))`.

The algebra shows that `o_p(zeta_NT)` is the correct normalized objective-gap scale, but the emphasized
enlarged-class curvature premise is not supplied by Revision 9. Hence this is not an RR2 corollary under the
unchanged assumptions.

## Resolution paths for a later paper-level decision

Exactly one of the following would be needed before a freeze:

1. strengthen `a:identification` to the cap+1 difference class (a substantive assumption change);
2. assume the generic pilot operator rate directly for the implemented cap+1 pilot (a high-level statistical
   pilot assumption, not the requested derived corollary); or
3. replace the cap+1 construction by a different extra-spectrum device and prove its rate under the existing
   assumptions.

RR2 does not choose among these paths because the task prohibits a new substantive assumption and requires
the cap+1 corollary to be derived.

## Numerical versus theoretical statements

- Statistical requirement: the displayed operator rate.
- Sufficient optimization condition: normalized global cap+1 objective gap `o_p(zeta_NT)`.
- Observable diagnostics: feasibility, all start objectives, best stable objective gap, stationarity/KKT
  residual, box activity, numerical rank, deterministic starts, and runtime.

Those diagnostics make the theoretical optimization condition credible; they do not mathematically certify
the unknown nonconvex global infimum. The condition is not weakened for computational convenience.
