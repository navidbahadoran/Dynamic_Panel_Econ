# Exact rank-stress DGP in the current implementation

This specification transcribes the current code without changing its formulas. Let `L=T+50` be the full burn-in plus observed horizon and let all displayed stress matrices first be constructed on `N x L`; the final `T` columns are retained. All added variables below are independent of the baseline DGP primitives and, block by block, of one another. Current component strengths equal one.

## Baseline `(1,1,1)`

The rank-stress generator first draws the main rank-one DGP:

\[
A^{(1),raw}_{it}=\lambda_{a,i}f_{a,t},\qquad
B^{(1)}_{it}=\lambda_{b,i}f_{b,t},\qquad
H^{(1),raw}_{it}=\lambda_{h,i}g_{h,t}.
\]

The first-term factors obey the main recursions

\[
g_{m,t}=0.5g_{m,t-1}+\sqrt{0.75}\,v_{m,t},\quad
v_{m,t}\stackrel{iid}{\sim}U[-\sqrt3,\sqrt3],\quad g_{m,-50}=0,
\]

for `m=a,b,h`, with `f_a=0.5+0.1g_a` and, for Revision 9, `f_b=0.6+0.15g_b`. The H loading satisfies `lambda_h,i iid U[-sqrt(3),sqrt(3)]`; A/B loadings are exactly those of the relevant main DGP.

For rank one no stress rescaling is applied. The A stability scale is

\[
c_a=\min\left\{1,\frac{0.85}{\max_{i,t\le L}|A^{(1),raw}_{it}|}\right\},
\quad A^{(1)}=c_aA^{(1),raw}.
\]

The final H is `H_0=c_xi*c_H*H^(1),raw`. These matrices have rank `(1,1,1)` almost surely.

## `(2,1,1)`

Only A receives a second term. Draw

\[
\widetilde\lambda_{A,i}\stackrel{iid}{\sim}U[-\sqrt3,\sqrt3],\qquad
\widetilde f_{A,t}\stackrel{iid}{\sim}U[-\sqrt3,\sqrt3].
\]

The second factor is iid over time; it has no AR recursion. Before common rescaling, its coefficient is one, so the two A terms are

\[
\lambda_{a,i}f_{a,t}
\quad\text{and}\quad
\widetilde\lambda_{A,i}\widetilde f_{A,t}.
\]

Let

\[
E_{A,raw}=L_{A,\max}(0.5+3\cdot0.1),
\]

where `L_A,max=1+0.1sqrt(3)` in DGPs 1--3 and `L_A,max=1.1+0.08sqrt(3)` in DGP4. Since the added outer product has support envelope 3, the code uses one common deterministic multiplier

\[
s_A=\frac{E_{A,raw}}{E_{A,raw}+3},\qquad
A^{(2),raw}_{it}=s_A\left(\lambda_{a,i}f_{a,t}
+\widetilde\lambda_{A,i}\widetilde f_{A,t}\right).
\]

It then recomputes stability over the full `N x L` realization:

\[
c_a=\min\left\{1,\frac{0.85}{\max_{i,t\le L}|A^{(2),raw}_{it}|}\right\},
\qquad A^{(2)}=c_aA^{(2),raw}.
\]

There is no clipping. B and H remain the rank-one baseline matrices. The same positive common multiplier applies to both A terms, so their pre-rescaling relative coefficient remains `1:1`; their support envelopes differ because the first term has envelope `E_A,raw` and the added term has envelope 3. The resulting ranks are exactly `(2,1,1)` almost surely.

## `(1,0,2)`

A remains the baseline rank-one A, including its standard full-horizon `c_a` scale. B is identically zero:

\[
B_{it}\equiv0.
\]

For H, draw

\[
\widetilde\lambda_{H,i}\stackrel{iid}{\sim}U[-\sqrt3,\sqrt3],\qquad
\widetilde g_{H,t}\stackrel{iid}{\sim}U[-\sqrt3,\sqrt3].
\]

The first H factor `g_h` is the bounded AR(1) above; the second factor is iid over time and has no recursion. With `E_H=3sqrt(3)`, the code sets

\[
s_H=\frac{3\sqrt3}{3\sqrt3+3},\qquad
H^{(2),raw}_{it}=s_H\left(\lambda_{h,i}g_{h,t}
+\widetilde\lambda_{H,i}\widetilde g_{H,t}\right).
\]

The population observed-entry variance used for calibration is

\[
V_{H,2}(T)=s_H^2\left\{
\frac1T\sum_{k=51}^{50+T}(1-0.5^{2k})+1
\right\},
\]

and

\[
c_H=\left[\frac{0.30}{0.70}\frac{1}{V_{H,2}(T)}\right]^{1/2}.
\]

The final interactive effect is `H_0=c_H*H^(2),raw`, because the established zero-B normalization is exactly `c_xi=1`. The scale equation is marked `r2_scale_identified=false`; 0.65 is not an intended target for this vector, and the pooled R2 induced by the normalization is reported. The resulting ranks are exactly `(1,0,2)` almost surely.

## Calibration separation

Every `DGP x N x T x true-rank` stress cell uses its own deterministic calibration draws with semantic seed label `rank_stress_calibration`, master seed `20260807`, and 50 draws in this candidate. No calibration draw is a reported replication.
