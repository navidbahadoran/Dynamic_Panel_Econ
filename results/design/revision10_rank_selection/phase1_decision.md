# Phase 1 decision

Classification: **ST2**.

A credible working selector exists: a conservative, monotone-upward, matrix-specific, self-normalized one-coordinate post-refit-gain procedure starting at the zero rank vector. It keeps all working Revision-9 estimation and inference layers and could transfer downstream inference on exact rank recovery.

One important theoretical package remains unresolved: construct and prove a tuning-free, uniformly conservative, rate-sharp observable envelope for the optimized block normal-space score under the paper's conditional temporal/spatial dependence, including mixed states where some blocks are underfit and others are at truth. A second localized requirement is hereditary one-coordinate detectability.

Therefore Phase 2/3 theorem drafting and implementation are **not authorized**. The next theory task is to resolve the envelope and induction without exceeding two new primitive conditions or strengthening the DGP to iid/Gaussian/independent/homoskedastic errors. Failure on that standard is THEORY COST TOO HIGH and triggers reconsideration of Route B.
