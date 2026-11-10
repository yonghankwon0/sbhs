################################################################################
#  Ablation experiment (Section 6.2, Table S1 of the supplementary material)
#
#  The Bethe-Hessian penalty decomposes as (Proposition 3.1)
#
#     b' H(r) b = b' L b  +  (r^2 - 1) ||b||^2  -  (r - 1) b' A b,
#                 Laplacian    ridge term           sign-coherence term
#
#  Is the reduction of false positives in the noise group due to the ridge
#  term? Three penalty matrices are compared on the data and folds of the main
#  simulation (r = adaptive r of each replicate):
#
#     SLS            L                    from results/simulation_metrics.csv
#     SLS + ridge    L + (r^2 - 1) I      fitted here
#     SBHS           H(r)                 from results/simulation_metrics.csv
#
#  All three use the same folds and the same (lambda1, lambda2) grid, so
#  "SLS + ridge" and SBHS differ only in the sign-coherence term.
#
#  Run from the repository root:   Rscript scripts/05_ablation.R   (about 3.5 minutes on 9 cores)
#  Input:   results/simulation_metrics.csv   (scripts/01_simulation.R)
#  Output:  results/ablation_metrics.csv     SLS + ridge, one row per replicate
#           results/tableS1_ablation.csv     Table S1: F1 and FP_G2, mean (SD)
################################################################################
source("R/sbhs.R")
source("R/simulation.R")

run_ablation <- function(scenario, n, b) {
  d <- replicate_data(scenario, n, b)
  X <- d$train$X; y <- d$train$y

  A <- signed_correlation_graph(X)
  r <- adaptive_r(A)
  L <- bethe_hessian(A, 1)
  bound_H <- lambda2_upper_bound(X, bethe_hessian(A, r), n_eig = 10)

  ridge <- L + (r^2 - 1) * diag(P)
  cv <- cv_sbhs(X, y, ridge, bound_H, seed = d$seed + 400)   # same folds and grid as SLS and SBHS
  fit <- sbhs_cd(X, y, ridge, cv$lambda1, cv$lambda2)
  cbind(scenario = scenario, n = n, rep = b, method = "SLS + ridge",
        evaluate(fit$beta, d$train$beta, d$test$X, d$test$y),
        r = r, lambda1 = cv$lambda1, lambda2 = cv$lambda2)
}

tasks <- expand.grid(rep = 1:100, n = c(60, 120, 180), scenario = SCENARIOS,
                     stringsAsFactors = FALSE)
cl <- sbhs_cluster(nrow(tasks))
parallel::clusterExport(cl, "run_ablation")
t0 <- Sys.time()
res <- parallel::parLapply(cl, seq_len(nrow(tasks)), function(i, tasks) {
  run_ablation(tasks$scenario[i], tasks$n[i], tasks$rep[i])
}, tasks = tasks)
parallel::stopCluster(cl)
message(sprintf("elapsed: %.1f minutes", as.numeric(difftime(Sys.time(), t0, units = "mins"))))

ablation <- do.call(rbind, res)
write.csv(ablation, "results/ablation_metrics.csv", row.names = FALSE)

# ---------------------------------------------------------------- Table S1
main <- read.csv("results/simulation_metrics.csv")
cols <- c("scenario", "n", "rep", "method", "F1", "FP_G2")
all <- rbind(main[main$method %in% c("SLS", "SBHS"), cols], ablation[, cols])

f3 <- function(x) sub("^0\\.", ".", sprintf("%.3f", x))
arms <- c("SLS", "SLS + ridge", "SBHS")
table_s1 <- NULL
for (sc in SCENARIOS) for (n in c(60, 120, 180)) {
  row <- data.frame(scenario = sc, n = n)
  for (arm in arms) {
    x <- all$F1[all$scenario == sc & all$n == n & all$method == arm]
    row[[paste("F1:", arm)]] <- sprintf("%s (%s)", f3(mean(x)), f3(sd(x)))
  }
  for (arm in arms) {
    x <- all$FP_G2[all$scenario == sc & all$n == n & all$method == arm]
    row[[paste("FP_G2:", arm)]] <- sprintf("%.1f (%.1f)", mean(x), sd(x))
  }
  table_s1 <- rbind(table_s1, row)
}
write.csv(table_s1, "results/tableS1_ablation.csv", row.names = FALSE)
print(table_s1, row.names = FALSE)
message("saved results/ablation_metrics.csv and results/tableS1_ablation.csv")
