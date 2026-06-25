library(coda)
library(tidyverse)

# Load it
chain <- readRDS("results/baseline_wave2.rds")

# Basic checks
dim(chain)               # rows = draws, cols = 6 parameters
colnames(chain)          # Beta1 Beta2 Beta3 Beta4 rho_y rho_o
effectiveSize(chain)     # all should be >= 200

# Posterior summary table
summary(chain)

# Trace plots (no density)
par(mfrow = c(3, 2), mar = c(3, 3, 2, 1))
for (p in colnames(chain)) {
  traceplot(chain[, p, drop = FALSE],
            main = p, col = "#154273", lwd = 0.4)
  abline(h = mean(chain[, p]), col = "#d73027", lty = 2)
}

# Posterior density plots
as.data.frame(chain) |>
  pivot_longer(everything(), names_to = "parameter", values_to = "value") |>
  ggplot(aes(x = value)) +
  geom_histogram(aes(y = after_stat(density)), bins = 50,
                 fill = "#154273", alpha = 0.6) +
  geom_density(colour = "black", linewidth = 0.8) +
  facet_wrap(~parameter, scales = "free") +
  theme_bw()

# Posterior means and 95% CIs
dom_eigen <- 4.074
as.data.frame(chain) |>
  pivot_longer(everything(), names_to = "parameter", values_to = "value") |>
  group_by(parameter) |>
  summarise(
    mean  = round(mean(value), 4),
    ci_lo = round(quantile(value, 0.025), 4),
    ci_hi = round(quantile(value, 0.975), 4)
  ) |>
  mutate(Rt = ifelse(grepl("Beta", parameter),
                     round(mean * dom_eigen, 3), NA)) |>
  print()
