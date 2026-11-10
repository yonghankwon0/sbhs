################################################################################
#  Tables 1 and 2: simulation results, mean (SD) over 100 replicates
#
#  Run from the repository root:   Rscript scripts/02_tables_simulation.R
#  Input:   results/simulation_metrics.csv   (scripts/01_simulation.R)
#  Output:  results/table1_selection.csv     F1, PPV, TPR, FP_G2
#           results/table2_prediction.csv    Test MSE / sigma^2, Est. MSE
#
#  Test MSE / sigma^2 is the test error divided by the noise variance of the
#  replicate, so a value of 1 is the error of the true coefficient vector.
#  FP_G2 is the number of selected variables in the noise Group 2.
################################################################################
m <- read.csv("results/simulation_metrics.csv")
m$test_mse_rel <- m$test_mse / m$sigma2

f3 <- function(x) sub("^0\\.", ".", sprintf("%.3f", x))          # .123 instead of 0.123
cell <- list(
  F1    = function(x) sprintf("%s (%s)", f3(mean(x)), f3(sd(x))),
  PPV   = function(x) sprintf("%s (%s)", f3(mean(x)), f3(sd(x))),
  TPR   = function(x) sprintf("%s (%s)", f3(mean(x)), f3(sd(x))),
  FP_G2 = function(x) sprintf("%.1f (%.1f)", mean(x), sd(x)),
  test_mse_rel = function(x) sprintf("%.3f (%.3f)", mean(x), sd(x)),
  est_mse      = function(x) sprintf("%.2f (%.2f)", mean(x), sd(x)))

make_table <- function(cols) {
  out <- NULL
  for (sc in c("dense", "sparse", "mixed")) for (n in c(60, 120, 180)) {
    for (method in c("MCP", "ENet", "SLS", "SBHS")) {
      d <- m[m$scenario == sc & m$n == n & m$method == method, ]
      stopifnot(nrow(d) == 100)
      row <- data.frame(scenario = sc, n = n, method = method)
      for (k in cols) row[[k]] <- cell[[k]](d[[k]])
      out <- rbind(out, row)
    }
  }
  out
}

table1 <- make_table(c("F1", "PPV", "TPR", "FP_G2"))
table2 <- make_table(c("test_mse_rel", "est_mse"))
write.csv(table1, "results/table1_selection.csv", row.names = FALSE)
write.csv(table2, "results/table2_prediction.csv", row.names = FALSE)

cat("Table 1: variable selection\n");            print(table1, row.names = FALSE)
cat("\nTable 2: prediction and estimation\n");   print(table2, row.names = FALSE)
