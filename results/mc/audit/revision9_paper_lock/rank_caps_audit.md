# Rank-cap audit

## Intended maintained experiments

| Experiment | Maintained configuration | Primary caps | True ranks inside cap set? |
|---|---|---:|---:|
| Main rank-one DGPs | `configs/mc/production.toml` | `(3,3,3)` | Yes: `(1,1,1) <= (3,3,3)` coordinatewise. |
| Separate rank-stress experiment | `configs/mc/rank_stress.toml` and `configs/mc/rank_stress_medium.toml` | `(3,3,3)` | Yes: `(1,1,1)`, `(2,1,1)`, and `(1,0,2)` are all coordinatewise inside. |

The common code default is also `(3,3,3)`. Both maintained configurations specify the vector explicitly, so there is no baseline-versus-rank-stress conflict.

Both maintained designs also declare a larger-cap sensitivity `(4,4,4)`. This is not the primary cap.

## Repository-wide conflict disclosure

Several historical preflight files specify `(2,2,2)`, including the finalized fixed/selected preflight and replay configurations. Those files conflict with the maintained `(3,3,3)` cap if treated as active designs. Their names and run purposes identify them as historical preflight artifacts, not the intended production or rank-stress configurations. No file was changed and this audit does not silently promote either value: the paper should name the authoritative configuration and state `(3,3,3)` explicitly.
