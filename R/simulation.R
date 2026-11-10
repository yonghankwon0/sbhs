################################################################################
#  Simulation design of Section 6.1 and the competing methods
#
#  Functions in this file
#     simulation_sigma()   population correlation matrix of a scenario
#     generate_data()      one data set (X, y) from a scenario
#     fit_mcp()            MCP, lambda by 5-fold CV (ncvreg)
#     fit_enet()           elastic net, alpha and lambda by 5-fold CV (glmnet)
#     evaluate()           selection, prediction and estimation metrics
#     run_replicate()      one replicate: data, four methods, metrics
#
#  Requires R/sbhs.R to be sourced first.
################################################################################

suppressPackageStartupMessages({
  library(mvtnorm)   # rmvnorm()
  library(ncvreg)    # MCP
  library(glmnet)    # elastic net
})

P          <- 120        # number of predictors
GROUP_SIZE <- 40         # three groups: 1-40, 41-80, 81-120
GROUP2     <- 41:80      # the noise group (no signal in any scenario)
SCENARIOS  <- c("dense", "sparse", "mixed")


# ------------------------------------------------------------------------------
# True coefficient vector of a scenario.
#   dense : all of Group 1 (1.5) and Group 3 (2.0) are signals      -> 80 signals
#   sparse: the first half of Group 1 and of Group 3 are signals    -> 40 signals
#   mixed : as dense, but the second half of each signal group has
#           the opposite sign                                       -> 80 signals
# Group 2 never carries signal.
# ------------------------------------------------------------------------------
true_beta <- function(scenario) {
  switch(scenario,
    dense  = c(rep(1.5, 40), rep(0, 40), rep(2.0, 40)),
    sparse = c(rep(1.5, 20), rep(0, 20), rep(0, 40), rep(2.0, 20), rep(0, 20)),
    mixed  = c(rep(1.5, 20), rep(-1.5, 20), rep(0, 40), rep(2.0, 20), rep(-2.0, 20)),
    stop("scenario must be one of: ", paste(SCENARIOS, collapse = ", ")))
}


# ------------------------------------------------------------------------------
# Population correlation matrix.
#
#   * correlation 0.6 within each group, 0 between groups;
#   * "mixed" scenario: correlation -0.4 between the two halves of Group 1 and
#     between the two halves of Group 3 (the halves have opposite signs);
#   * spurious correlations: every variable of Group 1 is paired with one
#     variable of the noise Group 2 (correlation 0.6), and so is every variable
#     of Group 3. The pairing is random and is determined by `seed`.
#
# The matrix built this way is not positive definite, so eigenvalues below
# 0.01 are raised to 0.01 and the result is rescaled to unit diagonal. The
# correlations actually used therefore differ slightly from the nominal values
# (Figure 3 shows the matrices that are used).
#
# NOTE: this function calls set.seed(seed). generate_data() relies on the state
# of the random number generator that it leaves behind.
# ------------------------------------------------------------------------------
simulation_sigma <- function(scenario, seed) {
  Sigma <- matrix(0, P, P)
  for (g in 1:3) {
    idx <- ((g - 1) * GROUP_SIZE + 1):(g * GROUP_SIZE)
    Sigma[idx, idx] <- 0.6
  }
  diag(Sigma) <- 1

  if (scenario == "mixed") {
    Sigma[1:20, 21:40] <- Sigma[21:40, 1:20] <- -0.4
    Sigma[81:100, 101:120] <- Sigma[101:120, 81:100] <- -0.4
  }

  # spurious pairs between the signal groups and the noise group
  set.seed(seed)
  g1 <- sample(1:40, 40);   g2a <- sample(41:80, 40)   # Group 1 <-> Group 2
  for (i in 1:40) Sigma[g1[i], g2a[i]] <- Sigma[g2a[i], g1[i]] <- 0.6
  g3 <- sample(81:120, 40); g2b <- sample(41:80, 40)   # Group 3 <-> Group 2
  for (i in 1:40) Sigma[g3[i], g2b[i]] <- Sigma[g2b[i], g3[i]] <- 0.6

  # nearest valid correlation matrix: clip small eigenvalues, rescale
  Sigma <- (Sigma + t(Sigma)) / 2
  eig <- eigen(Sigma, symmetric = TRUE)
  if (min(eig$values) < 0.01) {
    eig$values <- pmax(eig$values, 0.01)
    Sigma <- eig$vectors %*% diag(eig$values) %*% t(eig$vectors)
    Sigma <- (Sigma + t(Sigma)) / 2
  }
  D_inv <- diag(1 / sqrt(diag(Sigma)))
  Sigma <- D_inv %*% Sigma %*% D_inv
  diag(Sigma) <- 1
  Sigma
}


# ------------------------------------------------------------------------------
# One data set: X ~ N(0, Sigma), y = X beta + noise, with the noise variance
# chosen so that the signal-to-noise ratio beta' Sigma beta / sigma^2 equals
# `snr` (1 in the paper).
#
# `seed` fixes, in this order, the spurious pairs, X and the noise.
# ------------------------------------------------------------------------------
generate_data <- function(n, scenario, seed, snr = 1) {
  Sigma <- simulation_sigma(scenario, seed)      # sets the seed
  X <- mvtnorm::rmvnorm(n, sigma = Sigma)
  beta <- true_beta(scenario)
  sigma <- sqrt(as.numeric(t(beta) %*% Sigma %*% beta) / snr)
  y <- as.vector(X %*% beta + rnorm(n, 0, sigma))
  list(X = X, y = y, beta = beta, Sigma = Sigma, sigma = sigma)
}


# ------------------------------------------------------------------------------
# Competing methods without graph information
# ------------------------------------------------------------------------------

# MCP with gamma = 3; lambda minimises the 5-fold CV error.
fit_mcp <- function(X, y, seed, gamma = 3) {
  set.seed(seed)
  cv <- cv.ncvreg(X, y, family = "gaussian", penalty = "MCP", gamma = gamma,
                  nfolds = 5, max.iter = 50000)
  as.vector(coef(cv, s = "lambda.min"))[-1]      # drop the intercept
}

# Elastic net; alpha in {0, 0.1, ..., 1} and lambda minimise the 5-fold CV
# error. The same folds are used for every alpha.
fit_enet <- function(X, y, seed, alpha_seq = seq(0, 1, by = 0.1)) {
  set.seed(seed)
  foldid <- sample(rep(1:5, length.out = nrow(X)))
  cvs <- lapply(alpha_seq, function(a) cv.glmnet(X, y, alpha = a, foldid = foldid))
  best <- which.min(sapply(cvs, function(f) min(f$cvm)))
  as.numeric(coef(cvs[[best]], s = "lambda.min")[-1])
}


# ------------------------------------------------------------------------------
# Metrics of one fitted coefficient vector.
#   TPR, PPV, F1  variable selection (stored to four decimals)
#   FP_G2         number of selected variables in the noise Group 2
#   test_mse      prediction error on the test set
#   est_mse       mean squared error of the coefficients
# ------------------------------------------------------------------------------
evaluate <- function(beta_hat, beta, X_test, y_test) {
  truth <- which(beta != 0); sel <- which(beta_hat != 0)
  TP <- length(intersect(sel, truth))
  FP <- length(setdiff(sel, truth))
  FN <- length(setdiff(truth, sel))
  TPR <- if (TP + FN > 0) TP / (TP + FN) else 0
  PPV <- if (TP + FP > 0) TP / (TP + FP) else 0
  F1  <- if (PPV + TPR > 0) 2 * PPV * TPR / (PPV + TPR) else 0
  data.frame(TP = TP, FP = FP, TPR = round(TPR, 4), PPV = round(PPV, 4), F1 = round(F1, 4),
             FP_G2 = sum(beta_hat[GROUP2] != 0),
             test_mse = mean((y_test - X_test %*% beta_hat)^2),
             est_mse = mean((beta_hat - beta)^2))
}


# ------------------------------------------------------------------------------
# Training and test data of replicate b.
#
# The test set has the same coefficients and noise variance as the training
# set, but its spurious pairs are drawn afresh (Section 6.1): spurious
# correlations of the training sample are not expected to carry over to new
# data.
# ------------------------------------------------------------------------------
replicate_data <- function(scenario, n, b) {
  seed <- 1000 * b
  train <- generate_data(n, scenario, seed = seed + 10)
  test  <- generate_data(n, scenario, seed = seed + 1010)
  set.seed(seed + 2000)
  test$y <- as.vector(test$X %*% train$beta + rnorm(n, 0, train$sigma))
  list(train = train, test = test, seed = seed)
}


# ------------------------------------------------------------------------------
# One replicate of the comparison: MCP, ENet, SLS, SBHS.
#
# SLS and SBHS use the same folds (seed + 400) and the same (lambda1, lambda2)
# grid, so they differ only in the penalty matrix.
# Returns list(metrics = one row per method, beta = 4 x p coefficient matrix).
# ------------------------------------------------------------------------------
run_replicate <- function(scenario, n, b) {
  d <- replicate_data(scenario, n, b)
  X <- d$train$X; y <- d$train$y; seed <- d$seed

  sls  <- sbhs(X, y, penalty = "laplacian",     seed = seed + 400, n_eig = 10)
  bh   <- sbhs(X, y, penalty = "bethe-hessian", seed = seed + 400, n_eig = 10)
  beta_hat <- rbind(MCP  = fit_mcp(X, y, seed = seed + 100),
                    ENet = fit_enet(X, y, seed = seed + 200),
                    SLS  = sls$beta,
                    SBHS = bh$beta)

  metrics <- do.call(rbind, lapply(rownames(beta_hat), function(m) {
    cbind(scenario = scenario, n = n, rep = b, method = m,
          evaluate(beta_hat[m, ], d$train$beta, d$test$X, d$test$y))
  }))
  metrics$sigma2  <- d$train$sigma^2             # noise variance of this replicate
  metrics$r       <- c(NA, NA, 1, bh$r)
  metrics$lambda1 <- c(NA, NA, sls$lambda1, bh$lambda1)
  metrics$lambda2 <- c(NA, NA, sls$lambda2, bh$lambda2)
  list(metrics = metrics, beta = beta_hat)
}


# ------------------------------------------------------------------------------
# Helper: a PSOCK cluster whose workers have the SBHS functions loaded.
# The number of workers can be set with the environment variable SBHS_CORES.
# ------------------------------------------------------------------------------
sbhs_cluster <- function(n_tasks) {
  n_cores <- as.integer(Sys.getenv("SBHS_CORES", parallel::detectCores() - 1))
  n_cores <- max(1, min(n_cores, n_tasks))
  cl <- parallel::makeCluster(n_cores)
  parallel::clusterCall(cl, function(wd) {
    setwd(wd); source("R/sbhs.R"); source("R/simulation.R"); NULL
  }, getwd())
  message(sprintf("running %d tasks on %d cores", n_tasks, n_cores))
  cl
}
