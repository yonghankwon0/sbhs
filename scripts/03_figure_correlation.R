################################################################################
#  Figure 3: population correlation matrices of the three scenarios
#
#  Run from the repository root:   Rscript scripts/03_figure_correlation.R
#  Output:  results/figure3_correlation_structures.pdf
#
#  The matrices are those of the training data of the first replicate
#  (seed 1000 * 1 + 10, see replicate_data() in R/simulation.R).
################################################################################
source("R/sbhs.R")
source("R/simulation.R")
library(ggplot2)

label <- c(dense = "Dense", sparse = "Sparse", mixed = "Mixed signs")
cells <- do.call(rbind, lapply(SCENARIOS, function(sc) {
  S <- simulation_sigma(sc, seed = 1010)
  data.frame(i = rep(1:P, times = P), j = rep(1:P, each = P),
             value = as.vector(S), scenario = label[[sc]])
}))
cells$scenario <- factor(cells$scenario, levels = label)

# borders between the three groups
borders <- data.frame(x = c(40.5, 80.5, 0.5, 0.5), xend = c(40.5, 80.5, 120.5, 120.5),
                      y = c(0.5, 0.5, 40.5, 80.5), yend = c(120.5, 120.5, 40.5, 80.5))
mid <- c(20.5, 60.5, 100.5)

fig <- ggplot(cells, aes(j, i, fill = value)) +
  geom_raster() +
  geom_segment(data = borders, aes(x = x, xend = xend, y = y, yend = yend),
               inherit.aes = FALSE, colour = "#55544f", linewidth = 0.3) +
  scale_x_continuous(breaks = mid, labels = c("G1", "G2", "G3"), expand = c(0, 0)) +
  scale_y_reverse(breaks = mid, labels = c("G1", "G2", "G3"), expand = c(0, 0)) +
  scale_fill_gradient2(low = "#2166AC", mid = "#F7F7F5", high = "#B2182B", midpoint = 0,
                       limits = c(-1, 1), breaks = c(-1, -0.5, 0, 0.5, 1), name = "Correlation") +
  coord_fixed() +
  facet_wrap(~ scenario, nrow = 1) +
  labs(x = NULL, y = NULL) +
  theme_minimal(base_size = 9) +
  theme(panel.grid = element_blank(), axis.ticks = element_blank(),
        axis.text = element_text(colour = "#55544f"),
        strip.text = element_text(face = "bold", size = 9),
        legend.title = element_text(size = 8),
        legend.key.height = unit(16, "pt"), legend.key.width = unit(7, "pt"))

dir.create("results", showWarnings = FALSE)
ggsave("results/figure3_correlation_structures.pdf", fig, width = 6.5, height = 2.45)
message("saved results/figure3_correlation_structures.pdf")
