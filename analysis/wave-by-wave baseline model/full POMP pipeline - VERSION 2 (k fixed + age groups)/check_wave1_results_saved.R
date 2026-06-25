library(coda)
library(tidyverse)

# Load the saved posterior chain
chain <- readRDS("results/baseline_wave1.rds")

# Basic checks
dim(chain)            # rows = posterior draws, cols = 5 parameters
colnames(chain)       # Beta1 Beta2 Beta3 rho_y rho_o
effectiveSize(chain)  # all should be >= 200

# Trace plots — check mixing and stationarity
plot(chain)

# Posterior summary table
summary(chain)

# Tidy posterior means and 95% CIs
post_df <- as.data.frame(chain)

dom_eigen <- 4.074  # R0 = beta * 4.074

post_df |>
  pivot_longer(everything(), names_to = "parameter", values_to = "value") |>
  group_by(parameter) |>
  summarise(
    mean  = mean(value),
    sd    = sd(value),
    ci_lo = quantile(value, 0.025),
    ci_hi = quantile(value, 0.975)
  ) |>
  mutate(
    Rt = case_when(
      grepl("Beta", parameter) ~ round(mean * dom_eigen, 3),
      TRUE ~ NA_real_
    )
  ) |>
  print()

# Density plots per parameter
post_df |>
  pivot_longer(everything(), names_to = "parameter", values_to = "value") |>
  ggplot(aes(x = value)) +
  geom_histogram(aes(y = after_stat(density)),
                 bins = 50, fill = "#154273", alpha = 0.6) +
  geom_density(colour = "#00a1d5", linewidth = 0.9) +
  facet_wrap(~parameter, scales = "free") +
  theme_bw()

