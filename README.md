# Reproduction materials: from a panel of goodness-of-fit tests to one decision

Author: Ebrahim Khaled Ebrahim, Department of Statistics, Faculty of Science, Alexandria University
(ORCID 0009-0006-7839-8778).

These files reproduce every table of the paper "From a Panel of Goodness-of-Fit Tests to One
Decision: Redundancy, the Limits of Learned Weights, and Robustness in Logistic Regression" and
of its supplementary material.

## Contents

| Folder | What it holds |
|---|---|
| `pvalues/archive/` | The simulation archive. One file per role, model and sample size (`<role>_<model>_n<n>_pvalues.csv`). Each row is one simulated data set: its replicate number, its random seed, and the p-values of all tests computed on that same data set. Roles: `test` (power and size), `conformal` (the null reference set used for critical values), `train` (training sets used only to choose weights). |
| `pvalues/corruption/` | The corruption study: data from the correct model with k = 0, 1, 5, 10 corrupted covariate values, n = 500 and 1000. |
| `package_2.4.1_27d5ef2/` | The source files of the R package `ebrahim.gof` at the frozen commit used for every simulation (version 2.4.1). |
| `simulation/` | The data-generating models and the scripts that produced the archive (`run_stage1_features.R`), the corruption study (`run_corrupt.R`), the earlier held-out evaluation (`cand_nonneg.R`, `select_final_heldout.R`) and the timing (`exp_timing_frozen.R`). `check_reproduce_row.R` regenerates stored rows from their seeds. |
| `analysis/` | The Python scripts that turn the stored p-values into results (`result12_panel.py`, `result3_ceiling.py`, `result3b_transfer.py`, `result4_robust.py`, and `result5_rules_ci.py` for the other combination rules and the bootstrap intervals), their output files (`result*.csv`, `result*.txt`) next to them, and `make_tables.py`, which writes every LaTeX table into `analysis/tables/`. |
| `heldout_and_timing/` | The outputs behind Supplementary Table S2 and the computing times quoted in the Discussion. |

## How to reproduce

1. Tables from the stored p-values: from inside `analysis/`, run the five `result*.py` scripts, then
   `make_tables.py` (Python 3 with numpy, pandas and scipy). The first four take minutes; `result5_rules_ci.py`
   (1,000 bootstrap replicates) takes about ten minutes.
2. Regenerate stored data sets (optional): source `simulation/config.R`, then run
   `simulation/check_reproduce_row.R`, which rebuilds data sets from their seeds and compares the
   p-values with the stored ones.
3. The full archive (optional, many CPU hours): `simulation/run_stage1_features.R` and
   `simulation/run_corrupt.R`.

The R scripts use paths relative to the root of this folder, so set the working directory there
first. The Python scripts find their files from their own location and can be run from anywhere.
The simulations used R 4.4 on Windows.

## Licence

Code: MIT licence. Data (the p-value files and outputs): CC BY 4.0.
