# Numerical-fix validation

All four observations are exact-draw diagnostic replays, not new Monte Carlo replications.
No independent preflight, medium diagnostic, rank-stress, power, or production run was launched.

## Cap-pilot confirmation

```csv
semantic_replication_id,diagnostic_replay,c_kappa,original_route_count,original_route_objectives,original_best_objective,original_second_best_objective,original_stability_gap,confirmation_objectives,confirmation_objective_gaps,number_confirmation_valid,number_confirmation_matching_best,confirmed_best_basin,final_pilot_acceptance_basis,pilot_now_passes,replay_primary_status
dgp1_N50_T50_r00002_truth1-1-1,True,4e-06,6,"[1.4840205145923895, 1.4840205145923895, 1.4840205145923895, 1.3374322823826592, 1.4840205145923895, 1.4840098626330025]",1.3374322823826592,1.4840098626330025,0.1095962630640354,"[1.3374138854304949, 1.3374138557251798, 1.3374063971096273]","[1.3755427027340009e-05, 1.3777637733248563e-05, 1.935445508001242e-05]",3,0,False,failure,False,rank_pilot_failure
dgp3_N50_T50_r00000_truth1-1-1,True,4e-06,6,"[1.9866879291269903, 1.9866879291269903, 1.9866879291269903, 1.556584043245411, 1.9866946605935665, 1.5474148138510786]",1.5474148138510786,1.9866879291269903,0.2838754749818408,"[1.5473969703647856, 1.5473965285000546, 1.547395555428718]","[1.1531159022938604e-05, 1.1816709301423414e-05, 1.2445546073508337e-05]",3,0,False,failure,False,rank_pilot_failure
dgp3_N50_T50_r00002_truth1-1-1,True,4e-06,6,"[1.9637266912133622, 1.9637266912133622, 1.9637266912133622, 1.9637266912133622, 1.7152297315950344, 1.907626622252391]",1.907626622252391,1.9637266912133624,0.0294083067968155,"[1.9076121085299138, 1.9076121349091328, 1.90761253013998]","[7.608261652401447e-06, 7.5944333598395495e-06, 7.387248765878271e-06]",3,3,True,confirmed_best_basin,True,success
```

## Fixed-rank replay

```csv
semantic_replication_id,diagnostic_replay,objective_stability_pass,acceptance_basis,selected_envelope,stable_interior_solution_found,bound_failure_remains,replay_primary_status
dgp4_N50_T50_r00001_truth1-1-1,True,False,failure,14.532951298155336,False,True,coefficient_bound_hit
```

## Split-status correction

```csv
inconsistent_before,inconsistent_after,R_point_change,R_inference_change
20,0,0,-20
```

The 20 inconsistencies are Case A: at least one of the four required split coefficient fits failed strict numerical validity. Finite point estimates remain retained; inference does not.

`c_kappa=3e-6` is discontinued from further preflight consideration. `c_kappa=4e-6` remains the only candidate for one future independent validation and is not designated a production constant.
