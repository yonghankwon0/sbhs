################################################################################
#  Quick example: fit SBHS and SLS to one simulated data set (about 10 seconds)
#
#  Run from the repository root:   Rscript scripts/00_quick_example.R
################################################################################
source("R/sbhs.R")
source("R/simulation.R")

# One data set from the "dense" scenario: 120 predictors in three groups of 40.
# Groups 1 and 3 carry signal; Group 2 is noise that is spuriously correlated
# with the signal groups.
d <- generate_data(n = 120, scenario = "dense", seed = 3)
X <- d$X; y <- d$y

# SBHS: the graph is built from the sample correlations of X, r from the mean
# degree, and (lambda1, lambda2) by 5-fold cross-validation.
fit_bh <- sbhs(X, y, penalty = "bethe-hessian", seed = 3)

# SLS: the same graph, folds and tuning grid, but the Laplacian (r = 1).
fit_l  <- sbhs(X, y, penalty = "laplacian", seed = 3)

cat(sprintf("mean degree of the graph: %.1f,  adaptive r: %.3f\n",
            mean(rowSums(abs(fit_bh$A))), fit_bh$r))
cat(sprintf("lambda2 upper bound: %.4f,  selected lambda2: SBHS %.4f, SLS %.4f\n",
            fit_bh$lambda2_max, fit_bh$lambda2, fit_l$lambda2))

group <- rep(c("Group 1 (signal)", "Group 2 (noise)", "Group 3 (signal)"), each = 40)
cat("\nNumber of selected variables (out of 40 per group)\n")
print(rbind(SBHS = tapply(fit_bh$beta != 0, group, sum),
            SLS  = tapply(fit_l$beta  != 0, group, sum)))

# This is a single data set; the counts change from one data set to the next.
# Table 1 of the paper reports averages over 100 replicates per setting.

# To use your own graph, pass a signed adjacency matrix (entries -1, 0, 1):
#   fit <- sbhs(X, y, A = my_adjacency)
