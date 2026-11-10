# SBHS: Sparse Bethe-Hessian Shrinkage

R code for the paper

> Yonghan Kwon and Inkyung Jung.
> *Community-Aware Variable Selection via Bethe-Hessian Regularization in High-Dimensional Regression.*

SBHS is a penalized regression estimator for correlated predictors. It
combines the MCP sparsity penalty with a quadratic penalty built from the
Bethe-Hessian matrix of a signed graph over the predictors:

```
minimise over b:  (1/2n) ||y - X b||^2 + sum_j MCP(|b_j|; lambda1, gamma) + (lambda2/2) b' H(r) b
H(r) = (r^2 - 1) I - r A + D
```

`A` is the signed adjacency matrix of the graph and `D` its degree matrix.
With `r = 1` the matrix is the graph Laplacian and the estimator is Sparse
Laplacian Shrinkage (SLS; Huang et al., 2011).

This repository contains the implementation and the scripts that reproduce
every table and figure of the paper that is based on computation.

## Quick start

```r
source("R/sbhs.R")            # run from the repository root

fit <- sbhs(X, y)             # X: n x p matrix (centred columns), y: centred response
which(fit$beta != 0)          # selected variables
fit$r; fit$lambda1; fit$lambda2
```

`sbhs()` does four things, each of which is also available as a function:

| Step | Function | Paper |
|---|---|---|
| Signed graph from thresholded sample correlations | `signed_correlation_graph()` | Section 6.1 |
| Degree-adaptive parameter `r = 1 + 1 / log(mean degree + 1)` | `adaptive_r()` | Section 5.1 |
| Upper bound of `lambda2` (objective bounded below) | `lambda2_upper_bound()`, `is_bounded()` | Proposition 5.1 |
| 5-fold CV over `(lambda1, lambda2)` and coordinate descent | `cv_sbhs()`, `sbhs_cd()` | Section 5.2, Algorithm 1 |

Use `sbhs(X, y, penalty = "laplacian")` for SLS on the same graph, folds and
tuning grid, and `sbhs(X, y, A = my_graph)` to supply your own signed
adjacency matrix. A complete example on simulated data:

```bash
Rscript scripts/00_quick_example.R
```

## Requirements

R (version 4.5.2 was used) with a C++ compiler for Rcpp, and the packages

| Package | Version used | Needed for |
|---|---|---|
| Rcpp | 1.1.0 | coordinate descent |
| RSpectra | 0.16.2 | eigenvalues in the simulations |
| mvtnorm | 1.3.3 | simulated data |
| ncvreg | 3.16.0 | MCP |
| glmnet | 4.1.10 | elastic net |
| ggplot2 | 4.0.1 | figures |
| fuzzyforest | 1.0.8 | mouse liver data |
| hdi | 0.1.10 | riboflavin data |

```r
install.packages(c("Rcpp", "RSpectra", "mvtnorm", "ncvreg", "glmnet",
                   "ggplot2", "fuzzyforest", "hdi"))
```

## Reproducing the paper

Run every script from the repository root. The scripts that fit models run in
parallel; set the number of workers with the environment variable
`SBHS_CORES` (default: all cores but one).

| Result in the paper | Script | Output in `results/` | Time (9 cores) |
|---|---|---|---|
| Simulation study, Section 6.1-6.2 | `scripts/01_simulation.R` | `simulation_metrics.csv`, `simulation_beta.rds` | about 19 min |
| Tables 1 and 2 | `scripts/02_tables_simulation.R` | `table1_selection.csv`, `table2_prediction.csv` | seconds |
| Figure 3 | `scripts/03_figure_correlation.R` | `figure3_correlation_structures.pdf` | seconds |
| Figure 4 | `scripts/04_figure_paired_comparison.R` | `figure4_paired_comparison.pdf`, `figure4_counts.csv` | seconds |
| Ablation, Section 6.2 and Table S1 | `scripts/05_ablation.R` | `ablation_metrics.csv`, `tableS1_ablation.csv` | about 3.5 min |
| Real data, Section 6.3 and Table 3 | `scripts/06_real_data.R` | `liver_cv.csv`, `riboflavin_cv.csv`, `table3_real_data.csv` | about 2 min |

Scripts 02, 04 and 05 read the output of script 01. The output files of all
scripts are included in `results/`, so tables and figures can be rebuilt
without rerunning the simulations.

```bash
Rscript scripts/01_simulation.R
```

```bash
Rscript scripts/02_tables_simulation.R
```

### Reproducibility

All random steps are seeded inside the code (data generation, fold
assignments of every method), so the results do not depend on the number of
cores. Replicate `b` of the simulation uses the base seed `1000 * b`; the
real-data analyses use the seeds 2024-2033 for the ten repetitions. Exact
agreement of the last digits can depend on the platform and on package
versions (linear algebra library, `ncvreg`, `glmnet`).

The two figure files are drawn exactly as Figures 3 and 4 of the paper. The
files used in the paper additionally had their fonts embedded with
Ghostscript, which does not change how they look.

## Data

No data files are stored in this repository.

- **Mouse liver** (Ghazalpour et al., 2006): `Liver_Expr` in the CRAN package
  `fuzzyforest`. 17 of its 66 rows are exact copies of other rows; they are
  removed, which leaves 49 mice.
- **Riboflavin** (Bühlmann, Kalisch and Meier, 2014): `riboflavin` in the CRAN
  package `hdi`.

## Files

```
R/sbhs.R              SBHS and SLS: graph, Bethe-Hessian, lambda2 bound, CV, coordinate descent
R/simulation.R        simulation design, MCP and elastic net wrappers, metrics
src/sbhs_cd.cpp       one coordinate-descent sweep (C++)
scripts/              one script per table or figure (see above)
results/              output of the scripts
LICENSE               MIT license
```

## Citation

If you use this code, please cite the paper above.

## License

The code is released under the MIT license (see `LICENSE`).
