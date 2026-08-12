# Table and figure reproducibility audit

The N=100 raw bundles and aggregated replication records contain the fields required for the currently planned outputs:

| Planned output | Required recorded fields | Status |
|---|---|---|
| Bias/RMSE figures | target, replication ID, replication-specific truth, estimate, point-valid flag | Complete |
| Coverage figures | truth, corrected estimate, standard error, interval coverage, inference-valid flag | Complete |
| SE/SD figures | replication-specific error, standard error, common inference-valid set | Complete |
| Retention/failure figures | attempted, point/inference flags, primary status, failure detail, boundary/interiority fields | Complete |
| Target tables | target, applicability, truth, estimate, SE, coverage, rejection, correction components | Complete |
| Optimization/retention tables | fit type/start, convergence, stationarity, box/KKT, rank and split diagnostics | Complete |
| Runtime tables | task, full-fit, split-fit, aggregate inference, worker/resource and I/O timing | Complete with limitation below |

No planned scientific statistic is missing. Raw target rows also retain semantic IDs, realization hashes, seed information, coefficient envelopes, Riesz/Gram diagnostics, split assignments and diagnostics, variance estimates, and extreme-error flags sufficient to reconstruct the reported tables without rerunning the cell.

Timing limitation: the frozen instrumentation records aggregate inference time but does not separately clock the Riesz solve and variance construction. Those two phase times therefore cannot be separated. The Monte Carlo was not rerun and scientific code was not changed to add presentation-only clocks.
