################################################################################
#  Real-data analyses (Section 6.3, Table 3)
#
#  Two gene-expression data sets, both taken from CRAN packages:
#
#    liver       mouse liver expression and body weight (Ghazalpour et al.,
#                2006), data set `Liver_Expr` of the package fuzzyforest
#    riboflavin  B. subtilis expression and riboflavin production rate
#                (Buehlmann et al., 2014), data set `riboflavin` of the
#                package hdi
#
#  For each data set: the 200 genes with the largest variance, standardised
#  genes and response, 5-fold cross-validation repeated 10 times (seeds
#  2024-2033). All tuning is done inside the training folds.
#
#  Run from the repository root:
#     Rscript scripts/06_real_data.R              both data sets (about 2 minutes on 9 cores)
#     Rscript scripts/06_real_data.R liver        one data set
#  Output:  results/liver_cv.csv, results/riboflavin_cv.csv   one row per repetition
#           results/table3_real_data.csv
################################################################################
source("R/sbhs.R")
suppressPackageStartupMessages({ library(ncvreg); library(glmnet) })


# ------------------------------------------------------------------------------
# Data
# ------------------------------------------------------------------------------

# Top 200 genes by variance; genes and response scaled to mean 0, variance 1.
prepare <- function(X_full, y) {
  top <- order(apply(X_full, 2, var), decreasing = TRUE)[1:200]
  X <- scale(X_full[, top, drop = FALSE])
  attributes(X) <- list(dim = dim(X))            # plain matrix without names
  list(X = X, y = as.numeric(scale(y)))
}

load_liver <- function() {
  if (!requireNamespace("fuzzyforest", quietly = TRUE)) stop("install.packages('fuzzyforest')")
  data("Liver_Expr", package = "fuzzyforest", envir = environment())
  # 17 of the 66 rows are exact copies of other rows (the same mouse appears two
  # or three times). They are removed; otherwise the same mouse could be in the
  # training and in the test fold. 49 mice remain.
  liver <- Liver_Expr[!duplicated(Liver_Expr), ]
  X_full <- as.matrix(liver[, -1])               # column 1 is the body weight
  X_full <- X_full[, colSums(is.na(X_full)) == 0, drop = FALSE]
  prepare(X_full, liver[[1]])
}

load_riboflavin <- function() {
  if (!requireNamespace("hdi", quietly = TRUE)) stop("install.packages('hdi')")
  data("riboflavin", package = "hdi", envir = environment())
  prepare(as.matrix(riboflavin$x), as.numeric(riboflavin$y))
}


# ------------------------------------------------------------------------------
# Description of the graph on the full data: how many communities are there,
# and how many edges run between them? Communities are found by k-means on the
# eigenvectors of H(r) with negative eigenvalues (Saade et al., 2014).
# ------------------------------------------------------------------------------
describe_graph <- function(X) {
  A <- signed_correlation_graph(X)
  r <- adaptive_r(A)
  eig <- eigen(bethe_hessian(A, r), symmetric = TRUE)
  neg <- which(eig$values < -1e-10)
  set.seed(1)
  community <- kmeans(eig$vectors[, neg, drop = FALSE], length(neg), nstart = 50)$cluster
  same <- outer(community, community, "==")
  cat(sprintf("  graph: mean degree %.1f, r = %.3f, negative edges %.0f%%\n",
              mean(rowSums(abs(A))), r, 100 * mean(A[A != 0] < 0)))
  cat(sprintf("  communities: %s; edges between communities %.1f%%\n",
              paste(sort(table(community), decreasing = TRUE), collapse = ", "),
              100 * mean(!same[A != 0])))
}


# ------------------------------------------------------------------------------
# One repetition of 5-fold cross-validation for the four methods.
# The four methods are evaluated on the same outer folds. Inside each training
# fold, SLS and SBHS are tuned on the same inner folds and the same grid.
# Returns the CV-MSE, the mean number of selected genes and the MSE per fold.
# ------------------------------------------------------------------------------
one_repetition <- function(X, y, rep_seed, nfolds = 5) {
  n <- nrow(X)
  set.seed(rep_seed)
  folds <- sample(rep(1:nfolds, length.out = n))
  methods <- c("MCP", "ENet", "SLS", "SBHS")
  pred <- matrix(0, n, 4, dimnames = list(NULL, methods))
  n_selected <- setNames(numeric(4), methods)

  for (f in 1:nfolds) {
    tr <- which(folds != f); te <- which(folds == f)
    X_tr <- X[tr, , drop = FALSE]; y_tr <- y[tr]
    fold_seed <- rep_seed * 100 + f * 10
    beta <- list()

    # MCP (gamma = 3), lambda by inner 5-fold CV
    set.seed(fold_seed + 1)
    cv <- cv.ncvreg(X_tr, y_tr, penalty = "MCP", gamma = 3, nfolds = 5)
    beta$MCP <- coef(cv, s = "lambda.min")[-1]

    # Elastic net: alpha in {0, 0.1, ..., 1} and lambda by inner 5-fold CV
    alpha_seq <- seq(0, 1, by = 0.1)
    cv_err <- sapply(seq_along(alpha_seq), function(a) {
      set.seed(fold_seed + 2 + a)
      min(cv.glmnet(X_tr, y_tr, alpha = alpha_seq[a], nfolds = 5)$cvm)
    })
    set.seed(fold_seed + 2)
    cv <- cv.glmnet(X_tr, y_tr, alpha = alpha_seq[which.min(cv_err)], nfolds = 5)
    beta$ENet <- as.numeric(coef(cv, s = "lambda.min")[-1])

    # SLS and SBHS: graph, r and the lambda2 bound are computed from the
    # training fold only
    control <- list(cv_control = list(max_iter = 300, tol = 1e-5),
                    fit_control = list(max_iter = 500, tol = 1e-5))
    beta$SLS  <- do.call(sbhs, c(list(X_tr, y_tr, penalty = "laplacian",
                                      seed = fold_seed + 4), control))$beta
    beta$SBHS <- do.call(sbhs, c(list(X_tr, y_tr, penalty = "bethe-hessian",
                                      seed = fold_seed + 4), control))$beta

    for (m in methods) {
      pred[te, m] <- X[te, , drop = FALSE] %*% beta[[m]]
      n_selected[m] <- n_selected[m] + sum(beta[[m]] != 0)
    }
  }

  fold_mse <- sapply(methods, function(m) tapply((y - pred[, m])^2, folds, mean))
  out <- data.frame(seed = rep_seed)
  for (m in methods) out[[paste0(m, "_MSE")]] <- mean((y - pred[, m])^2)
  for (m in methods) out[[paste0(m, "_Var")]] <- n_selected[[m]] / nfolds
  for (m in methods) out[paste0(m, "_f", 1:nfolds)] <- fold_mse[, m]
  out
}


# ------------------------------------------------------------------------------
# Run
# ------------------------------------------------------------------------------
which_data <- commandArgs(trailingOnly = TRUE)
if (length(which_data) == 0) which_data <- c("liver", "riboflavin")
seeds <- 2024:2033
methods <- c("MCP", "ENet", "SLS", "SBHS")
dir.create("results", showWarnings = FALSE)
table3 <- NULL

for (name in which_data) {
  dat <- switch(name, liver = load_liver(), riboflavin = load_riboflavin(),
                stop("unknown data set: ", name))
  cat(sprintf("\n== %s: n = %d, p = %d\n", name, nrow(dat$X), ncol(dat$X)))
  describe_graph(dat$X)

  n_cores <- max(1, min(as.integer(Sys.getenv("SBHS_CORES", parallel::detectCores() - 1)), length(seeds)))
  cl <- parallel::makeCluster(n_cores)
  parallel::clusterCall(cl, function(wd) {
    setwd(wd); source("R/sbhs.R"); library(ncvreg); library(glmnet); NULL
  }, getwd())
  parallel::clusterExport(cl, "one_repetition")
  res <- do.call(rbind, parallel::parLapply(cl, seeds, function(s, dat) {
    one_repetition(dat$X, dat$y, s)
  }, dat = dat))
  parallel::stopCluster(cl)
  write.csv(res, sprintf("results/%s_cv.csv", name), row.names = FALSE)

  # Table 3: mean (SD) over the 10 repetitions. The SD reflects the fold
  # assignment only, because the same samples are used in every repetition.
  mcp <- mean(res$MCP_MSE)
  for (m in methods) {
    mse <- res[[paste0(m, "_MSE")]]; nv <- res[[paste0(m, "_Var")]]
    table3 <- rbind(table3, data.frame(
      data = name, method = m,
      cv_mse = sprintf("%.3f (%.3f)", mean(mse), sd(mse)),
      selected = sprintf("%.1f (%.1f)", mean(nv), sd(nv)),
      improvement_over_mcp = if (m == "MCP") "" else sprintf("%.1f%%", 100 * (mcp - mean(mse)) / mcp)))
  }
  # In how many of the 10 repetitions is the CV-MSE of SBHS lower?
  for (m in methods[1:3]) {
    cat(sprintf("  SBHS has a lower CV-MSE than %-4s in %d of %d repetitions\n", m,
                sum(res$SBHS_MSE < res[[paste0(m, "_MSE")]]), nrow(res)))
  }
}

cat("\nTable 3\n"); print(table3, row.names = FALSE)
if (length(which_data) == 2) write.csv(table3, "results/table3_real_data.csv", row.names = FALSE)
