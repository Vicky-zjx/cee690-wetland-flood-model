# Direct wetland coverage-bound checks

A/B are synthetic response classes. All values below come from the public `analysis(v,c)` function and official `ode4u`; no optimizer is run. Negative physical/scaled inequality residuals satisfy their corresponding limit. Feasibility uses the existing declared numerical tolerance. The integration-status residual is -1 for a valid completed run. The maximum residual also includes the storage constraints, so a value close to zero can indicate a nearly empty storage rather than an upstream safety limit.

| Pair | Covers | Step (s) | Status | Physical river peak (m) | Max scaled residual | Max balance error (m³) | Violated constraints |
|---|---|---:|---|---:|---:|---:|---|
| A/A | 0, 0 | 120 | valid_feasible | 97.594269322 | -4.04208e-37 | 0.0341611 | None |
| A/A | 0, 0 | 60 | valid_feasible | 97.594277061 | -4.04207e-37 | 0.00854072 | None |
| A/A | 0, 1 | 120 | valid_feasible | 97.583811314 | -4.04208e-37 | 0.0500543 | None |
| A/A | 0, 1 | 60 | valid_feasible | 97.583811302 | -4.04207e-37 | 0.0125136 | None |
| A/A | 1, 0 | 120 | valid_feasible | 97.577488545 | -4.04208e-37 | 0.0554297 | None |
| A/A | 1, 0 | 60 | valid_feasible | 97.577488533 | -4.04207e-37 | 0.0138574 | None |
| A/A | 1, 1 | 120 | valid_feasible | 97.567207102 | -4.04208e-37 | 0.078958 | None |
| A/A | 1, 1 | 60 | valid_feasible | 97.567207086 | -4.04207e-37 | 0.0197391 | None |
| A/B | 0, 0 | 120 | valid_feasible | 97.594269322 | -4.04208e-37 | 0.0341611 | None |
| A/B | 0, 0 | 60 | valid_feasible | 97.594277061 | -4.04207e-37 | 0.00854072 | None |
| A/B | 0, 1 | 120 | valid_feasible | 97.588363257 | -4.04208e-37 | 0.0530525 | None |
| A/B | 0, 1 | 60 | valid_feasible | 97.588363245 | -4.04207e-37 | 0.0132631 | None |
| A/B | 1, 0 | 120 | valid_feasible | 97.577488545 | -4.04208e-37 | 0.0554297 | None |
| A/B | 1, 0 | 60 | valid_feasible | 97.577488533 | -4.04207e-37 | 0.0138574 | None |
| A/B | 1, 1 | 120 | valid_feasible | 97.571673287 | -4.04208e-37 | 0.083532 | None |
| A/B | 1, 1 | 60 | valid_feasible | 97.571676590 | -4.04207e-37 | 0.0208826 | None |
| B/A | 0, 0 | 120 | valid_feasible | 97.594269322 | -4.04208e-37 | 0.0341611 | None |
| B/A | 0, 0 | 60 | valid_feasible | 97.594277061 | -4.04207e-37 | 0.00854072 | None |
| B/A | 0, 1 | 120 | valid_feasible | 97.583811314 | -4.04208e-37 | 0.0500543 | None |
| B/A | 0, 1 | 60 | valid_feasible | 97.583811302 | -4.04207e-37 | 0.0125136 | None |
| B/A | 1, 0 | 120 | valid_feasible | 97.580990162 | -4.04208e-37 | 0.0595123 | None |
| B/A | 1, 0 | 60 | valid_feasible | 97.581000502 | -4.04207e-37 | 0.014878 | None |
| B/A | 1, 1 | 120 | valid_feasible | 97.570618540 | -4.04208e-37 | 0.0823172 | None |
| B/A | 1, 1 | 60 | valid_feasible | 97.570621591 | -4.04207e-37 | 0.0205789 | None |
| B/B | 0, 0 | 120 | valid_feasible | 97.594269322 | -4.04208e-37 | 0.0341611 | None |
| B/B | 0, 0 | 60 | valid_feasible | 97.594277061 | -4.04207e-37 | 0.00854072 | None |
| B/B | 0, 1 | 120 | valid_feasible | 97.588363257 | -4.04208e-37 | 0.0530525 | None |
| B/B | 0, 1 | 60 | valid_feasible | 97.588363245 | -4.04207e-37 | 0.0132631 | None |
| B/B | 1, 0 | 120 | valid_feasible | 97.580990162 | -4.04208e-37 | 0.0595123 | None |
| B/B | 1, 0 | 60 | valid_feasible | 97.581000502 | -4.04207e-37 | 0.014878 | None |
| B/B | 1, 1 | 120 | valid_feasible | 97.575150561 | -4.04208e-37 | 0.0875865 | None |
| B/B | 1, 1 | 60 | valid_feasible | 97.575150547 | -4.04207e-37 | 0.0218962 | None |

## Timestep comparison

| Pair | Covers | Absolute peak change, 120 to 60 s (mm) | Same validity/feasibility classification | 30 s needed |
|---|---|---:|---|---|
| type_a/type_a | [0.0, 0.0] | 0.007738 | True | False |
| type_a/type_a | [0.0, 1.0] | 0.000012 | True | False |
| type_a/type_a | [1.0, 0.0] | 0.000012 | True | False |
| type_a/type_a | [1.0, 1.0] | 0.000015 | True | False |
| type_a/type_b | [0.0, 0.0] | 0.007738 | True | False |
| type_a/type_b | [0.0, 1.0] | 0.000012 | True | False |
| type_a/type_b | [1.0, 0.0] | 0.000012 | True | False |
| type_a/type_b | [1.0, 1.0] | 0.003302 | True | False |
| type_b/type_a | [0.0, 0.0] | 0.007738 | True | False |
| type_b/type_a | [0.0, 1.0] | 0.000012 | True | False |
| type_b/type_a | [1.0, 0.0] | 0.010340 | True | False |
| type_b/type_a | [1.0, 1.0] | 0.003052 | True | False |
| type_b/type_b | [0.0, 0.0] | 0.007738 | True | False |
| type_b/type_b | [0.0, 1.0] | 0.000012 | True | False |
| type_b/type_b | [1.0, 0.0] | 0.010340 | True | False |
| type_b/type_b | [1.0, 1.0] | 0.000013 | True | False |
