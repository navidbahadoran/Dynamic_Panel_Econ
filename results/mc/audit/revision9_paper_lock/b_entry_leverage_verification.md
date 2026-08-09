# Revision-9 B-entry leverage verification

For the bounded AR recursion

\[
g_{b,t}=0.5g_{b,t-1}+\sqrt{0.75}\,v_{b,t},
\qquad v_{b,t}\in[-\sqrt3,\sqrt3],
\]

the infinite-horizon deterministic envelope is

\[
|g_{b,t}|\leq
\frac{\sqrt{0.75}\sqrt3}{1-0.5}=3.
\]

Thus `support(g_b,t)` is contained in `[-3,3]` uniformly in the horizon. Under Revision 9,

\[
f_{b,t}=0.6+0.15g_{b,t}\in[0.15,1.05].
\]

The lower bound 0.15 is deterministic and uniform; it does not rely on sampled minima.

For rank-one `B=lambda_b f_b'`, the normalized time singular vector is `v=f_b/||f_b||`. At every paper entry date,

\[
T v_t^2=\frac{T f_{b,t}^2}{\sum_{s=1}^T f_{b,s}^2}
\geq \frac{0.15^2}{1.05^2}=\frac1{49}.
\]

This supplies the required uniform time-side leverage floor. The unit loadings also have positive deterministic floors: `1-0.4sqrt(3)>0` in DGPs 1--3 and `0.8-0.25sqrt(3)>0` in the lower-mean DGP4 group. Therefore the Revision-9 B-entry target is compatible with the paper's rank-one target-regularity assumption.

**Result: MATCH / THEOREM-COMPATIBLE.**
