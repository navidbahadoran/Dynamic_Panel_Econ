# Boundary and selected-cap inference compliance

## Coefficient-bound behavior

| Question | Current code | Paper compliance |
|---|---|---|
| Is the constrained point estimator computed with the literal box? | Yes. Fixed-rank fitting uses the coefficient bound; when the unconstrained iterate is outside, the constrained solver enforces the linearized entrywise box and checks feasibility/KKT diagnostics. Nuclear screening uses a nuclear-plus-box proximal map. | MATCH for the constrained point estimator. |
| Are boundary-active point estimates retained? | Yes. `boundary_active` is recorded, and finite point estimates remain eligible for bias/RMSE accounting. | MATCH with the paper's distinction between point estimation and interiority diagnostics. |
| Is primary inference suppressed when the bound is active? | No. `classify_inference_status` has no full-fit or split-fit `boundary_active` failure rule. A boundary-active fit can remain `success`, obtain a standard error, and enter coverage/rejection calculations if the other diagnostics pass. | MISMATCH. |

The implementation defines boundary activity numerically when the fitted maximum absolute coefficient is within ten constraint tolerances of `B`. It records boundary flags for full and split fits, but the flags are diagnostic only.

Recommended order: first make the paper explicit about whether *any* active full/split bound suppresses an interval or only the selected full fit does. Then change status precedence and tests. Preserve the finite constrained point estimate while setting primary inference invalid.

## Selected-rank cap behavior

The selected-rank rule matches the paper. In selected-rank mode, if any selected component equals its imposed cap, `selected_rank_at_cap=true`, the replication receives `rank_at_cap`, and execution returns a rank record plus a failure record before full/split target inference. Therefore no primary interval is reported.

This is exact for the baseline imposed cap. Larger-cap sensitivity is diagnostic and does not override the primary early stop.

## Related numerical caveat

Candidate post-refits enforce two-valid-start objective agreement, with a third start when needed; unresolved candidates receive infinite IC and cannot win. The cap pilot records route agreement and attempts confirmation of the best basin, but if both route agreement and confirmation fail, it currently returns the best valid route with a disagreement warning. If the paper intends objective stability to be a hard eligibility condition for the cap pilot itself, that point is underspecified in the supplied paper specification and should be resolved in the paper before code is changed.
