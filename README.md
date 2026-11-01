# Reproduction materials: from a panel of goodness-of-fit tests to one decision

Author: Ebrahim Khaled Ebrahim, Department of Applied Statistics, Faculty of Business, Alexandria University
(ORCID 0009-0006-7839-8778).

These files reproduce every table and figure of the paper "From a Panel of Goodness-of-Fit Tests to One
Decision: Redundancy, the Limits of Learned Weights, and Robustness in Logistic Regression" and
of its supplementary material.

## Contents

| Folder | What it holds |
|---|---|
| `pvalues/archive/` | The simulation archive. One file per role, model and sample size (`<role>_<model>_n<n>_pvalues.csv`). Each row is one simulated data set: its replicate number, its random seed, and the p-values of all tests computed on that same data set. Roles: `test` (power and size), `conformal` (the null reference set used for critical values), `train` (training sets used only to choose weights). |
| `pvalues/corruption/` | The corruption study: data from the correct model with k = 0, 1, 5, 10 corrupted covariate values, n = 500 and 1000. |
| `pvalues/extra/` | The extra simulations of Section 3.2 and the end of Section 6 of the paper: other recording errors (`corrupt_flipy_*`, `corrupt_x2_*`), two other null designs (`design_<rare or six>_<reference or test>_n*`), and the parametric-bootstrap calibration (`boot_*` in the original design, `bootD_*` in the two other designs, `bootE_*` for power at n = 200; each holds the observed p-values and the bootstrap p-values of each combination rule). |
| `package_2.4.1_27d5ef2/` | The source files of the R package `ebrahim.gof` at the frozen commit used for every simulation (version 2.4.1). |
| `package_legoft_panels/` | The source of `ebrahim.gof` with `legoft()` extended to the panels and combination rules of the paper (`panel = "equal"`, `"full"`, `"classical"`, `"one.per.group"`, `"protected"`; `rule = "cauchy"`, `"hmp"`, `"fisher"`, `"minp"`), as committed for version 3.0.0 of the package. Install with `install.packages("package_legoft_panels/ebrahim.gof_3.0.0.tar.gz", repos = NULL, type = "source")`, or use version 3.0.0 or later from CRAN. |
| `simulation/` | The data-generating models and the scripts that produced the archive (`run_stage1_features.R`), the corruption study (`run_corrupt.R`), the extra simulations (`run_extra.R`, `run_extra2.R`), the earlier held-out evaluation (`cand_nonneg.R`, `select_final_heldout.R`) and the timing (`exp_timing_frozen.R`). `check_reproduce_row.R` and `check_reproduce_extra.R` regenerate stored rows from their seeds and compare them with the stored p-values. |
| `analysis/` | The Python scripts that turn the stored p-values into results (`result12_panel.py`, `result3_ceiling.py`, `result3b_transfer.py`, `result4_robust.py`, `result5_rules_ci.py` for the other combination rules and the bootstrap intervals, `result6_extra.py` and `result7_extra2.py` for the extra simulations, and `result8_legoft.py` for the shipped LEGofT weights), their output files (`result*.csv`, `result*.txt`) next to them, `make_tables.py`, which writes every LaTeX table into `analysis/tables/`, and `make_figures.py`, which draws every figure into `analysis/figures/`. |
| `heldout_and_timing/` | The outputs behind Supplementary Table S2 and the computing times quoted in the Discussion. |

## How to reproduce

1. Tables and figures from the stored p-values: from inside `analysis/`, run the eight `result*.py` scripts, then
   `make_tables.py` and `make_figures.py` (Python 3 with numpy, pandas, scipy and matplotlib). Most scripts take
   minutes; `result5_rules_ci.py` (1,000 bootstrap replicates) takes about ten minutes.
2. Regenerate stored data sets (optional): from the root of this folder, run
   `Rscript simulation/check_reproduce_row.R` and `Rscript simulation/check_reproduce_extra.R`. They rebuild
   data sets from their seeds and print the largest difference from the stored p-values (about 1e-16).
3. The full simulations (optional, many CPU hours): `simulation/run_stage1_features.R`,
   `simulation/run_corrupt.R`, `simulation/run_extra.R` and `simulation/run_extra2.R`. Every data set is written to disk as it finishes,
   and a stopped run resumes where it stopped.

The R scripts use paths relative to the root of this folder, so set the working directory there
first. The Python scripts find their files from their own location and can be run from anywhere.
The simulations used R 4.4 on Windows.

## Changes

- v1.0.3: adds `result8_legoft.py` (the shipped LEGofT weights on the archive) and the package source with the
  panels of the paper in `legoft()` (`package_legoft_panels/`).
- v1.0.2: adds the bootstrap calibration in the two other null designs and the paired comparison of
  bootstrap and fixed-reference power at n = 200 (`run_extra2.R`, `result7_extra2.py`, `pvalues/extra/bootD_*`,
  `bootE_*`), a per-alternative power table, and one rounding rule (half up) for every table and figure label.
- v1.0.1: adds the extra simulations (`pvalues/extra/`, `run_extra.R`, `result6_extra.py`), the figure script,
  and `check_reproduce_extra.R`. Corrects the R source paths of v1.0.0, which pointed to the folder layout of
  the original project rather than to this archive, so that the R scripts run from the root of this folder.
  Corrects the author's affiliation.
- v1.0.0: first release.

## Licence

Code: MIT licence. Data (the p-value files and outputs): CC BY 4.0.
