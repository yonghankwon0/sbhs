################################################################################
#  Figure 4: replicate-wise differences between SBHS and SLS
#
#  Run from the repository root:   Rscript scripts/04_figure_paired_comparison.R
#  Input:   results/simulation_metrics.csv   (scripts/01_simulation.R)
#  Output:  results/figure4_paired_comparison.pdf
#           results/figure4_counts.csv       wins, ties and losses of SBHS
#
#  Each dot is one replicate (SBHS and SLS were fitted to the same data with
#  the same folds). The horizontal bar is the mean difference. The numbers
#  above each column are "replicates in which SBHS is better : worse"; the
#  remaining replicates are ties. Better means a higher F1 (top row) or fewer
#  selected variables in the noise Group 2 (bottom row).
################################################################################
library(ggplot2)

m <- read.csv("results/simulation_metrics.csv")
key <- c("scenario", "n", "rep")
w <- merge(m[m$method == "SBHS", c(key, "F1", "FP_G2")],
           m[m$method == "SLS",  c(key, "F1", "FP_G2")], by = key, suffixes = c(".sbhs", ".sls"))

d <- rbind(data.frame(w[key], metric = "F1", diff = w$F1.sbhs - w$F1.sls),
           data.frame(w[key], metric = "FP", diff = w$FP_G2.sbhs - w$FP_G2.sls))
d$better <- ifelse(d$metric == "F1", d$diff > 0, d$diff < 0)     # SBHS better
d$worse  <- ifelse(d$metric == "F1", d$diff < 0, d$diff > 0)     # SBHS worse
d$who <- ifelse(d$diff == 0, "Tie", ifelse(d$better, "SBHS better", "SLS better"))
d$scenario <- factor(d$scenario, levels = c("dense", "sparse", "mixed"),
                     labels = c("Dense", "Sparse", "Mixed signs"))
d$metric <- factor(d$metric, levels = c("F1", "FP"),
                   labels = c("F1~~(SBHS-SLS)", "FP[G2]~~(SBHS-SLS)"))
# Row order used for the figure in the paper. The horizontal jitter of the dots is
# random and is drawn in row order, so this order makes the dots land in the same places.
d <- d[order(d$n, d$scenario, d$metric, d$rep), ]

counts <- aggregate(cbind(better, worse) ~ scenario + n + metric, d, sum)
counts$tie <- 100 - counts$better - counts$worse
top <- tapply(d$diff, d$metric, max)
counts$y <- top[as.character(counts$metric)] + ifelse(grepl("^F1", counts$metric), 0.18, 7)   # clear of the highest dot
counts$label <- paste0(counts$better, ":", counts$worse)
write.csv(counts[, c("scenario", "n", "metric", "better", "tie", "worse")],
          "results/figure4_counts.csv", row.names = FALSE)

set.seed(1)   # horizontal jitter of the dots only
fig <- ggplot(d, aes(factor(n), diff)) +
  geom_hline(yintercept = 0, colour = "#8f8d87", linewidth = 0.4) +
  geom_jitter(aes(colour = who), width = 0.13, height = 0, size = 1.2, alpha = 0.8) +
  stat_summary(fun = mean, geom = "crossbar", width = 0.5, linewidth = 0.3, colour = "#1f1e1c") +
  geom_text(data = counts, aes(factor(n), y, label = label), size = 2.5, colour = "#55544f") +
  facet_grid(metric ~ scenario, scales = "free_y", switch = "y",
             labeller = labeller(metric = label_parsed)) +
  scale_colour_manual(values = c("SBHS better" = "#2a78d6", "SLS better" = "#eb6834", "Tie" = "#8f8d87"),
                      breaks = c("SBHS better", "SLS better", "Tie"), name = NULL) +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.14))) +
  labs(x = "Sample size n", y = NULL) +
  theme_minimal(base_size = 9) +
  theme(panel.grid.minor = element_blank(), panel.grid.major.x = element_blank(),
        panel.grid.major.y = element_line(colour = "#e8e7e3", linewidth = 0.3),
        strip.placement = "outside", strip.text = element_text(face = "bold", size = 9),
        axis.text = element_text(colour = "#55544f"), legend.position = "top")

ggsave("results/figure4_paired_comparison.pdf", fig, width = 6.3, height = 4.2)
message("saved results/figure4_paired_comparison.pdf and results/figure4_counts.csv")
