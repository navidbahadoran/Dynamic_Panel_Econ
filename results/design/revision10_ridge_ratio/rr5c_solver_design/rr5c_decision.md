# RR5c decision

## RR5c-PASS-SEMANTIC-ENGINEERING

A plausible semantically exact route exists. The preferred architecture keeps the mathematical cap+1 pilot, coefficient box, normalized objective, cap widths, maintained starts and all acceptance thresholds unchanged. It targets the established numerical failures with:

- balanced product-preserving gauges and quotient-aware refinement on the interior path, whose stationarity diagnostic is product invariant;
- reversible internal preconditioning and a deterministic QP-specific active-set solver on the constrained path, mapped back before the frozen factor-KKT evaluation;
- monotone safeguards, deterministic condition telemetry and exact full-width/lower-rank handling;
- mandatory non-scientific solver tests, task-atomic resumability equivalence and resource-instrumented scaling benchmarks.

The dominant RR5 problem was not objective-agreement tolerance: 177/180 pilots never produced two individually valid starts. No implementation bug, wrong scale, feasibility failure or statistical rank failure was established. RR5 produced no accepted normalized spectrum or selected rank, so it is not evidence against the ridge-ratio rule.

Changing `B`, caps, starts, thresholds, objective agreement, ridge, rank-zero anchor, singular-value rules or estimator remains prohibited without a new paper-level decision. RR5c does not establish that semantically exact engineering will succeed; it establishes a specific, testable route that must be implemented and validated before any scientific rerun.

No Monte Carlo, DGP generation, tuning, solver implementation, benchmark, source/config/test/manuscript edit, or RR5 rerun occurred in RR5c.
