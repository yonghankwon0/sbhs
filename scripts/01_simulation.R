################################################################################
#  Simulation study (Section 6.1-6.2): input of Tables 1-2 and Figure 4
#
#  3 scenarios x 3 sample sizes x 100 replicates = 900 replicates, four methods
#  each (MCP, ENet, SLS, SBHS). About 19 minutes on 9 cores.
#
#  Run from the repository root:   Rscript scripts/01_simulation.R
#  Output:  results/simulation_metrics.csv   one row per replicate and method
#           results/simulation_beta.rds      estimated coefficients
#
#  Every replicate sets its own seeds (see replicate_data() and
#  run_replicate() in R/simulation.R), so the results do not depend on the
#  number of cores or on the order in which the replicates are run.
################################################################################
source("R/sbhs.R")
source("R/simulation.R")

tasks <- expand.grid(rep = 1:100, n = c(60, 120, 180), scenario = SCENARIOS,
                     stringsAsFactors = FALSE)

cl <- sbhs_cluster(nrow(tasks))
t0 <- Sys.time()
res <- parallel::parLapply(cl, seq_len(nrow(tasks)), function(i, tasks) {
  run_replicate(tasks$scenario[i], tasks$n[i], tasks$rep[i])
}, tasks = tasks)
parallel::stopCluster(cl)
message(sprintf("elapsed: %.1f minutes", as.numeric(difftime(Sys.time(), t0, units = "mins"))))

metrics <- do.call(rbind, lapply(res, `[[`, "metrics"))
# The two mean squared errors are stored to four decimals; the tables of the
# paper were computed from the stored values.
metrics$test_mse <- round(metrics$test_mse, 4)
metrics$est_mse  <- round(metrics$est_mse, 4)

dir.create("results", showWarnings = FALSE)
write.csv(metrics, "results/simulation_metrics.csv", row.names = FALSE)
saveRDS(list(tasks = tasks, beta = lapply(res, `[[`, "beta")), "results/simulation_beta.rds")
message("saved results/simulation_metrics.csv and results/simulation_beta.rds")
