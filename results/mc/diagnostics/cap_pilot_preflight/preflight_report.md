# Cap-pilot numerical-start preflight

Decision: **NO-GO** for medium diagnostics. The strict cap-pilot success rate was 66.7% (8/12), below the requested nearly-all gate.

Candidate coverage was 66.7%; all 12 true coefficient envelopes satisfied the deterministic interior condition. No medium or production run was started.

| dgp | replication | cap_pilot_success | attempted_routes | stable_routes | best_two_gap | thresholded_cap_rank | candidate_covers_111 | selected_rank | runtime_seconds | truth_envelope | cap_envelope_ratio |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| 1 | 0 | True | 6 | 4 | 0 | [1, 2, 1] | True | [0, 0, 0] | 129.48 | 4.13386 | 0.911999 |
| 1 | 1 | False | 6 | 1 | 0.206195 |  | False |  | 101.922 | 3.2028 |  |
| 1 | 2 | False | 6 | 1 | 0.0949079 |  | False |  | 81.9669 | 4.59421 |  |
| 2 | 0 | True | 6 | 5 | 3.6602e-06 | [1, 1, 2] | True | [0, 0, 0] | 103.438 | 4.32234 | 0.9225 |
| 2 | 1 | False | 6 | 1 | 0.129366 |  | False |  | 86.2245 | 4.3362 |  |
| 2 | 2 | True | 6 | 5 | 0 | [1, 1, 2] | True | [0, 0, 0] | 107.633 | 4.16729 | 0.710846 |
| 3 | 0 | True | 6 | 4 | 2.18702e-06 | [1, 1, 1] | True | [0, 0, 0] | 106.897 | 5.54148 | 0.875952 |
| 3 | 1 | True | 6 | 3 | 0 | [1, 1, 1] | True | [0, 0, 0] | 108.988 | 4.45622 | 0.677699 |
| 3 | 2 | True | 6 | 3 | 0 | [1, 1, 2] | True | [0, 0, 0] | 108.494 | 5.03029 | 0.99602 |
| 4 | 0 | False | 6 | 1 | 1.3047 |  | False |  | 77.536 | 4.14581 |  |
| 4 | 1 | True | 6 | 4 | 3.70858e-06 | [1, 1, 1] | True | [0, 0, 0] | 118.444 | 4.93909 | 0.99194 |
| 4 | 2 | True | 6 | 3 | 0 | [1, 1, 1] | True | [0, 0, 0] | 100.896 | 5.30044 | 0.809674 |

Full route objectives and machine-readable columns are in `tab_cap_pilot_preflight.csv` and `tab_cap_pilot_preflight.parquet`.
