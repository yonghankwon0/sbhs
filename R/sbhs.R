################################################################################
#  Sparse Bethe-Hessian Shrinkage (SBHS): core functions
#
#  SBHS minimises
#
#     (1/2n) ||y - X b||^2 + sum_j MCP(|b_j|; lambda1, gamma) + (lambda2/2) b' H(r) b,
#
#  where H(r) = (r^2 - 1) I - r A + D is the Bethe-Hessian matrix of a signed
#  graph over the predictors (A = signed adjacency matrix, D = degree matrix).
#  With r = 1 the matrix is the graph Laplacian L = D - A and the estimator is
#  Sparse Laplacian Shrinkage (SLS; Huang et al., 2011).
#
#  Typical use (see scripts/00_quick_example.R):
#
#     source("R/sbhs.R")
#     fit <- sbhs(X, y)          # graph, r, lambda1 and lambda2 chosen from the data
#     which(fit$beta != 0)       # selected variables
#
#  Functions in this file
#     signed_correlation_graph()  signed graph from thresholded sample correlations
#     adaptive_r()                degree-adaptive Bethe-Hessian parameter
#     bethe_hessian()             H(r); r = 1 gives the Laplacian
#     lambda2_upper_bound()       lambda2^max of Proposition 5.1
#     is_bounded()                is the objective bounded below for a given lambda2?
#     sbhs_cd()                   coordinate descent for given (lambda1, lambda2)
#     cv_sbhs()                   K-fold cross-validation over the (lambda1, lambda2) grid
#     sbhs()                      everything above in one call
################################################################################

suppressPackageStartupMessages(library(Rcpp))

# Compile the C++ coordinate-descent sweep once per R session. `SBHS_HOME` can
# be set beforehand to the repository root if this file is sourced from
# another working directory.
if (!exists("SBHS_HOME")) SBHS_HOME <- "."
if (!exists("cd_sweep")) Rcpp::sourceCpp(file.path(SBHS_HOME, "src", "sbhs_cd.cpp"))


# ------------------------------------------------------------------------------
# Signed graph from sample correlations (Section 6.1, "Graph construction")
#
# Two predictors are joined by an edge when their sample correlation is larger
# in absolute value than a threshold tau; the edge carries the sign of the
# correlation. tau is the two-sided critical value of a correlation test at
# level alpha (Fisher's z-transformation), capped at `cap` so that small
# samples do not produce an empty graph.
#
#   X      n x p matrix of predictors
#   alpha  test level (0.001 in the paper)
#   cap    upper limit for the threshold (0.5 in the paper)
# Returns the p x p signed adjacency matrix with entries in {-1, 0, 1}.
# ------------------------------------------------------------------------------
signed_correlation_graph <- function(X, alpha = 0.001, cap = 0.5) {
  n <- nrow(X)
  tau <- min(tanh(qnorm(1 - alpha / 2) / sqrt(n - 3)), cap)
  R <- cor(X)
  A <- (abs(R) > tau) * sign(R)
  diag(A) <- 0
  A
}


# ------------------------------------------------------------------------------
# Degree-adaptive Bethe-Hessian parameter (Section 5.1):
#     r_hat = 1 + 1 / log(mean degree + 1)
# Dense graphs give r_hat close to 1 (close to the Laplacian); sparse graphs
# give a larger r_hat. For a graph without edges the formula is undefined and
# 1.05 is returned.
# ------------------------------------------------------------------------------
adaptive_r <- function(A) {
  d_bar <- sum(abs(A)) / nrow(A)
  if (d_bar == 0) return(1.05)
  1 + 1 / log(d_bar + 1)
}


# ------------------------------------------------------------------------------
# Bethe-Hessian matrix H(r) = (r^2 - 1) I - r A + D  (Definition 2.1).
# bethe_hessian(A, 1) is the (signed) graph Laplacian D - A.
# ------------------------------------------------------------------------------
bethe_hessian <- function(A, r) {
  p <- nrow(A)
  d <- rowSums(abs(A))
  (r^2 - 1) * diag(p) - r * A + diag(d, p)
}


# ------------------------------------------------------------------------------
# Upper bound lambda2^max of Proposition 5.1:
#     min over negative eigenvalues mu_i of H of  (||X u_i||^2 / n) / |mu_i|.
# For every lambda2 above this value the objective is unbounded below, so it is
# used as the top of the lambda2 grid. Returns Inf when H has no negative
# eigenvalue (for example the Laplacian).
#
#   n_eig  if not NULL and the package RSpectra is installed, only the n_eig
#          smallest eigenvalues are computed (faster for large p). The
#          simulations of the paper used n_eig = 10; the real-data analyses used
#          the full eigendecomposition (n_eig = NULL).
# ------------------------------------------------------------------------------
lambda2_upper_bound <- function(X, H, n_eig = NULL) {
  n <- nrow(X); p <- ncol(X)
  if (!is.null(n_eig) && requireNamespace("RSpectra", quietly = TRUE)) {
    eig <- RSpectra::eigs_sym(H, k = min(n_eig, p - 1), which = "SA")
  } else {
    eig <- eigen(H, symmetric = TRUE)
  }
  neg <- which(eig$values < -1e-10)
  if (length(neg) == 0) return(Inf)
  min(sapply(neg, function(i) {
    sum((X %*% eig$vectors[, i])^2) / n / abs(eig$values[i])
  }))
}


# ------------------------------------------------------------------------------
# Is the SBHS objective bounded below for this lambda2?  (Proposition 5.1)
#
# The objective is bounded below exactly when G(lambda2) = X'X/n + lambda2 H is
# positive semidefinite (up to the boundary case). X is scaled column-wise in
# the same way as in sbhs_cd(), because that is the scale on which the penalty
# acts. lambda2 = 0 is always admissible.
# ------------------------------------------------------------------------------
is_bounded <- function(X, H, lambda2, tol = 1e-8) {
  if (lambda2 == 0) return(TRUE)
  n <- nrow(X)
  cn <- sqrt(colSums(X^2)); cn[cn < 1e-10] <- 1
  G0 <- crossprod(sweep(X, 2, cn / sqrt(n), "/")) / n
  min(eigen(G0 + lambda2 * H, symmetric = TRUE, only.values = TRUE)$values) >= -tol
}


# ------------------------------------------------------------------------------
# Coordinate descent for fixed tuning parameters (Algorithm 1).
#
#   X, y       data (y and the columns of X are assumed to be centred)
#   H          penalty matrix, e.g. bethe_hessian(A, r)
#   lambda1    MCP (sparsity) parameter
#   lambda2    weight of the quadratic graph penalty
#   gamma      MCP shape parameter (3 throughout the paper)
#   max_iter   maximum number of passes over the coordinates
#   tol        stop when no coefficient changes by more than tol in a pass
#
# The columns of X are rescaled internally to ||x_j||^2 = n and the returned
# coefficients are transformed back to the scale of the input X.
# Returns list(beta, iterations, converged).
# ------------------------------------------------------------------------------
sbhs_cd <- function(X, y, H, lambda1, lambda2, gamma = 3, max_iter = 1000, tol = 1e-6) {
  X <- as.matrix(X); y <- as.numeric(y); H <- as.matrix(H)
  n <- nrow(X); p <- ncol(X)

  col_norms <- sqrt(colSums(X^2))
  col_norms[col_norms < 1e-10] <- 1            # guard against constant columns
  X_scaled <- sweep(X, 2, col_norms / sqrt(n), "/")

  H_diag <- diag(H)
  beta <- rep(0, p)                            # start from the null model
  residual <- y

  for (iter in seq_len(max_iter)) {
    sweep_res <- cd_sweep(X_scaled, H, beta, residual, H_diag,
                          lambda1, lambda2, gamma, as.integer(n), as.integer(p))
    beta <- as.numeric(sweep_res$beta)
    residual <- as.numeric(sweep_res$residual)
    if (sweep_res$max_change < tol) break
  }

  list(beta = beta / (col_norms / sqrt(n)),    # back to the original scale of X
       iterations = iter,
       converged = (iter < max_iter))
}


# ------------------------------------------------------------------------------
# K-fold cross-validation over the (lambda1, lambda2) grid of Section 5.2.
#
#   lambda1 grid: 15 values, log-spaced from lambda1^max = max_j |x_j'y| / n
#                 down to 0.01 * lambda1^max
#   lambda2 grid: 0 and 7 values, log-spaced from `lambda2_max` down to
#                 0.01 * lambda2_max; values at which the objective is not
#                 bounded below are discarded (is_bounded())
#
#   H            penalty matrix (the same matrix is used in every fold)
#   lambda2_max  top of the lambda2 grid, normally lambda2_upper_bound(X, H).
#                For the Laplacian (no negative eigenvalue) pass the bound of
#                the Bethe-Hessian so that both methods are tuned on the same
#                grid, as in the paper.
#   seed         seed for the fold assignment
#   max_iter, tol  settings of sbhs_cd() inside the folds
#
# Returns the selected pair, the full grid and the mean CV error at each grid
# point.
# ------------------------------------------------------------------------------
cv_sbhs <- function(X, y, H, lambda2_max, nfolds = 5, gamma = 3, seed = NULL,
                    max_iter = 500, tol = 1e-5) {
  n <- nrow(X)
  lambda1_max <- max(abs(crossprod(X, y))) / n
  lambda1_seq <- lambda1_max * 10^seq(0, -2, length.out = 15)

  top <- if (is.finite(lambda2_max)) lambda2_max else lambda1_max
  lambda2_seq <- c(0, top * 10^seq(0, -2, length.out = 7))
  lambda2_seq <- lambda2_seq[sapply(lambda2_seq, function(l2) is_bounded(X, H, l2))]

  if (!is.null(seed)) set.seed(seed)
  folds <- sample(rep(1:nfolds, length.out = n))
  grid <- expand.grid(l1 = lambda1_seq, l2 = lambda2_seq)

  cv_err <- matrix(NA_real_, nrow(grid), nfolds)
  for (i in seq_len(nrow(grid))) {
    for (f in 1:nfolds) {
      tr <- which(folds != f); te <- which(folds == f)
      fit <- sbhs_cd(X[tr, , drop = FALSE], y[tr], H, grid$l1[i], grid$l2[i], gamma, max_iter, tol)
      if (!anyNA(fit$beta)) {
        cv_err[i, f] <- mean((y[te] - X[te, , drop = FALSE] %*% fit$beta)^2)
      }
    }
  }
  mean_err <- rowMeans(cv_err, na.rm = TRUE)
  best <- which.min(mean_err)

  list(lambda1 = grid$l1[best], lambda2 = grid$l2[best],
       grid = grid, cv_error = mean_err,
       lambda1_seq = lambda1_seq, lambda2_seq = lambda2_seq)
}


# ------------------------------------------------------------------------------
# SBHS (or SLS) in one call: graph -> r -> penalty matrix -> CV -> final fit.
#
#   penalty  "bethe-hessian" (SBHS) or "laplacian" (SLS). SLS is tuned on the
#            grid of SBHS, i.e. the lambda2 grid is built from the upper bound
#            of the Bethe-Hessian, so that the two methods differ only in r.
#   A        optional signed adjacency matrix. By default the graph is built
#            from the sample correlations of X.
#   n_eig    passed to lambda2_upper_bound()
#   cv_control, fit_control  lists with max_iter and tol for the coordinate
#            descent inside CV and for the final fit
#
# Returns the coefficients and everything that was chosen along the way.
# ------------------------------------------------------------------------------
sbhs <- function(X, y, penalty = c("bethe-hessian", "laplacian"), A = NULL,
                 gamma = 3, nfolds = 5, seed = NULL, n_eig = NULL,
                 cv_control = list(max_iter = 500, tol = 1e-5),
                 fit_control = list(max_iter = 1000, tol = 1e-6)) {
  penalty <- match.arg(penalty)
  if (is.null(A)) A <- signed_correlation_graph(X)

  r_hat <- adaptive_r(A)
  H_bh <- bethe_hessian(A, r_hat)
  lambda2_max <- lambda2_upper_bound(X, H_bh, n_eig = n_eig)
  H <- if (penalty == "bethe-hessian") H_bh else bethe_hessian(A, 1)

  cv <- cv_sbhs(X, y, H, lambda2_max, nfolds = nfolds, gamma = gamma, seed = seed,
                max_iter = cv_control$max_iter, tol = cv_control$tol)
  fit <- sbhs_cd(X, y, H, cv$lambda1, cv$lambda2, gamma,
                 fit_control$max_iter, fit_control$tol)

  list(beta = fit$beta, lambda1 = cv$lambda1, lambda2 = cv$lambda2,
       lambda2_max = lambda2_max, r = if (penalty == "bethe-hessian") r_hat else 1,
       A = A, H = H, cv = cv, converged = fit$converged, penalty = penalty)
}
