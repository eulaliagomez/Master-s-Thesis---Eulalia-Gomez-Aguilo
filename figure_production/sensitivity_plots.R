# =============================================================================
# sensitivity_plots.R
#
# Visualisation of sensitivity analysis results
# Run AFTER sensitivity_pipeline.R has completed.
#
# Produces 5 figures:
#   Fig 1: Posterior mean shift by parameter and degradation type
#   Fig 2: CI width ratio by parameter and degradation type
#   Fig 3: Coverage heatmap (all waves × parameters × scenarios)
#   Fig 4: Systematic vs stochastic comparison (isolate noise effect)
#   Fig 5: Theoretical prediction vs observed (rho shift under systematic)
# =============================================================================

library(tidyverse)
library(patchwork)

# Load results
metrics_df  <- readRDS("results/sensitivity_metrics.rds")
summary_df  <- readRDS("results/sensitivity_summary.rds")

# Colour palette (CBS-inspired)
col_sys  <- "#154273"   # navy   — systematic
col_sto  <- "#00a1d5"   # blue   — stochastic
col_agg  <- "#e17000"   # orange — aggregation
pal_deg  <- c("Systematic reduction"      = col_sys,
              "Stochastic under-reporting" = col_sto,
              "Temporal aggregation"       = col_agg)

wave_labels <- c(wave1 = "Wave 1\n(Ancestral)",
                 wave2 = "Wave 2\n(Double hump)",
                 wave3 = "Wave 3\n(Alpha)",
                 wave4 = "Wave 4\n(Delta)")

param_labels <- c(Beta1 = "beta[1]", Beta2 = "beta[2]",
                  Beta3 = "beta[3]", Beta4 = "beta[4]",
                  rho_y = "rho[y]",  rho_o = "rho[o]")

# ── Figure 1: Posterior mean shift ───────────────────────────────────────────
fig1 <- metrics_df |>
  filter(degradation != "Temporal aggregation") |>   # continuous f levels
  mutate(
    wave      = factor(wave, levels = names(wave_labels), labels = wave_labels),
    parameter = factor(parameter, levels = names(param_labels),
                       labels = param_labels)
  ) |>
  ggplot(aes(x = level_num, y = mean_shift_pct,
             colour = degradation, group = degradation)) +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "grey60",
             linewidth = 0.5) +
  geom_line(linewidth = 0.8) +
  geom_point(size = 2.2) +
  facet_grid(parameter ~ wave, scales = "free_y",
             labeller = labeller(parameter = label_parsed,
                                 wave      = label_value)) +
  scale_colour_manual(values = pal_deg, name = NULL) +
  scale_x_continuous(breaks = c(0.10, 0.30, 0.50),
                     labels = c("10%", "30%", "50%")) +
  labs(
    title    = "Posterior mean shift under data degradation",
    subtitle = "Relative shift (%) vs clean-data baseline",
    x        = "Degradation level (fraction removed)",
    y        = "Mean shift (% of baseline)"
  ) +
  theme_bw(base_size = 10) +
  theme(legend.position  = "bottom",
        strip.text.y     = element_text(angle = 0),
        panel.grid.minor = element_blank())

# Add temporal aggregation as separate panel annotation
fig1_agg <- metrics_df |>
  filter(degradation == "Temporal aggregation") |>
  mutate(
    wave      = factor(wave, levels = names(wave_labels), labels = wave_labels),
    parameter = factor(parameter, levels = names(param_labels),
                       labels = param_labels),
    block_lab = ifelse(level == "2", "2-week", "4-week")
  ) |>
  ggplot(aes(x = block_lab, y = mean_shift_pct, fill = block_lab)) +
  geom_col(width = 0.5) +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "grey60") +
  facet_grid(parameter ~ wave, scales = "free_y",
             labeller = labeller(parameter = label_parsed,
                                 wave      = label_value)) +
  scale_fill_manual(values = c("2-week" = "#e17000", "4-week" = "#b35500"),
                    name = NULL) +
  labs(
    title    = "Posterior mean shift — temporal aggregation",
    subtitle = "Relative shift (%) vs clean-data baseline",
    x        = "Aggregation block",
    y        = "Mean shift (% of baseline)"
  ) +
  theme_bw(base_size = 10) +
  theme(legend.position  = "bottom",
        strip.text.y     = element_text(angle = 0),
        panel.grid.minor = element_blank())

ggsave("results/fig1_mean_shift.pdf",    fig1,     width=14, height=10)
ggsave("results/fig1b_mean_shift_agg.pdf", fig1_agg, width=14, height=10)
cat("Fig 1 saved\n")

# ── Figure 2: CI width ratio ──────────────────────────────────────────────────
fig2_cont <- metrics_df |>
  filter(degradation != "Temporal aggregation") |>
  mutate(
    wave      = factor(wave, levels = names(wave_labels), labels = wave_labels),
    parameter = factor(parameter, levels = names(param_labels),
                       labels = param_labels)
  ) |>
  ggplot(aes(x = level_num, y = ci_width_ratio,
             colour = degradation, group = degradation)) +
  geom_hline(yintercept = 1, linetype = "dashed", colour = "grey60",
             linewidth = 0.5) +
  geom_line(linewidth = 0.8) +
  geom_point(size = 2.2) +
  facet_grid(parameter ~ wave, scales = "free_y",
             labeller = labeller(parameter = label_parsed,
                                 wave      = label_value)) +
  scale_colour_manual(values = pal_deg, name = NULL) +
  scale_x_continuous(breaks = c(0.10, 0.30, 0.50),
                     labels = c("10%", "30%", "50%")) +
  scale_y_continuous(limits = c(NA, NA)) +
  labs(
    title    = "95% credible interval width ratio under data degradation",
    subtitle = "Ratio > 1 = wider (more uncertain) than baseline",
    x        = "Degradation level (fraction removed)",
    y        = "CI width ratio"
  ) +
  theme_bw(base_size = 10) +
  theme(legend.position  = "bottom",
        strip.text.y     = element_text(angle = 0),
        panel.grid.minor = element_blank())

fig2_agg <- metrics_df |>
  filter(degradation == "Temporal aggregation") |>
  mutate(
    wave      = factor(wave, levels = names(wave_labels), labels = wave_labels),
    parameter = factor(parameter, levels = names(param_labels),
                       labels = param_labels),
    block_lab = ifelse(level == "2", "2-week", "4-week")
  ) |>
  ggplot(aes(x = block_lab, y = ci_width_ratio, fill = block_lab)) +
  geom_col(width = 0.5) +
  geom_hline(yintercept = 1, linetype = "dashed", colour = "grey60") +
  facet_grid(parameter ~ wave, scales = "free_y",
             labeller = labeller(parameter = label_parsed,
                                 wave      = label_value)) +
  scale_fill_manual(values = c("2-week" = "#e17000", "4-week" = "#b35500"),
                    name = NULL) +
  labs(
    title    = "CI width ratio — temporal aggregation",
    subtitle = "Ratio > 1 = wider than baseline",
    x        = "Aggregation block",
    y        = "CI width ratio"
  ) +
  theme_bw(base_size = 10) +
  theme(legend.position  = "bottom",
        strip.text.y     = element_text(angle = 0),
        panel.grid.minor = element_blank())

ggsave("results/fig2_ci_ratio.pdf",       fig2_cont, width=14, height=10)
ggsave("results/fig2b_ci_ratio_agg.pdf",  fig2_agg,  width=14, height=10)
cat("Fig 2 saved\n")

# ── Figure 3: Coverage heatmap ────────────────────────────────────────────────
# Each cell: does the clean-data 95% CI contain the degraded posterior mean?

fig3 <- metrics_df |>
  mutate(
    scenario  = paste0(as.character(degradation), "\n(", level, ")"),
    wave      = factor(wave, levels = names(wave_labels), labels = wave_labels),
    parameter = factor(parameter, levels = names(param_labels),
                       labels = param_labels)
  ) |>
  ggplot(aes(x = scenario, y = parameter, fill = coverage)) +
  geom_tile(colour = "white", linewidth = 0.5) +
  geom_text(aes(label = ifelse(coverage, "\u2713", "\u2717")),
            size = 3.5, colour = "white", fontface = "bold") +
  facet_wrap(~wave, nrow = 1) +
  scale_fill_manual(values = c("TRUE" = "#154273", "FALSE" = "#c0392b"),
                    labels = c("TRUE" = "Covered", "FALSE" = "Not covered"),
                    name = "Baseline CI\ncovers degraded mean") +
  scale_y_discrete(labels = parse(text = param_labels)) +
  labs(
    title    = "Coverage: does the clean-data 95% CI contain the degraded posterior mean?",
    subtitle = "Tick = yes (robust); Cross = no (baseline CI too narrow for degraded estimate)",
    x        = NULL,
    y        = NULL
  ) +
  theme_bw(base_size = 10) +
  theme(axis.text.x   = element_text(angle = 45, hjust = 1, size = 8),
        legend.position = "right",
        panel.grid    = element_blank())

ggsave("results/fig3_coverage.pdf", fig3, width = 16, height = 6)
cat("Fig 3 saved\n")

# ── Figure 4: Systematic vs stochastic (isolate noise effect) ─────────────────
# Side-by-side: at same f level, how much does binomial noise add?
# Focus on rho_y and rho_o since theory predicts beta is unaffected.

fig4 <- metrics_df |>
  filter(degradation %in% c("Systematic reduction",
                             "Stochastic under-reporting")) |>
  filter(parameter %in% c("rho_y", "rho_o",
                           "Beta1", "Beta2", "Beta3")) |>
  mutate(
    wave      = factor(wave, levels = names(wave_labels), labels = wave_labels),
    parameter = factor(parameter, levels = names(param_labels),
                       labels = param_labels)
  ) |>
  ggplot(aes(x = level_num, y = ci_width_ratio,
             colour = degradation, linetype = degradation, group = degradation)) +
  geom_hline(yintercept = 1, linetype = "dashed", colour = "grey70",
             linewidth = 0.4) +
  geom_line(linewidth = 0.9) +
  geom_point(size = 2.5) +
  facet_grid(parameter ~ wave, scales = "free_y",
             labeller = labeller(parameter = label_parsed,
                                 wave      = label_value)) +
  scale_colour_manual(values = c("Systematic reduction"      = col_sys,
                                 "Stochastic under-reporting" = col_sto),
                      name = NULL) +
  scale_linetype_manual(values = c("Systematic reduction"      = "solid",
                                   "Stochastic under-reporting" = "dashed"),
                        name = NULL) +
  scale_x_continuous(breaks = c(0.10, 0.30, 0.50),
                     labels = c("10%", "30%", "50%")) +
  labs(
    title    = "Systematic vs stochastic under-reporting: CI width ratio",
    subtitle = paste("Gap between lines = additional uncertainty from binomial noise.",
                     "Beta parameters should be unaffected by systematic reduction."),
    x        = "Degradation level (fraction removed)",
    y        = "CI width ratio"
  ) +
  theme_bw(base_size = 10) +
  theme(legend.position  = "bottom",
        strip.text.y     = element_text(angle = 0),
        panel.grid.minor = element_blank())

ggsave("results/fig4_sys_vs_sto.pdf", fig4, width = 14, height = 10)
cat("Fig 4 saved\n")

# ── Figure 5: Theoretical prediction vs observed ──────────────────────────────
# Under systematic reduction, theory predicts:
#   rho_a* ≈ rho_a / (1 - f)  → post_mean should scale by 1/(1-f)
#   beta_j* ≈ beta_j           → ratio should be 1.0
# Plot observed ratio (degraded / baseline) vs theoretical prediction.

baseline_means <- summary_df |>
  filter(degradation == "baseline") |>
  select(wave, parameter, base_mean = post_mean)

fig5_data <- summary_df |>
  filter(degradation == "systematic") |>
  left_join(baseline_means, by = c("wave", "parameter")) |>
  mutate(
    f             = as.numeric(level),
    observed_ratio = post_mean / base_mean,
    predicted_rho  = 1 / (1 - f),   # theory for rho
    predicted_beta = 1,              # theory for beta
    predicted      = ifelse(grepl("rho", parameter),
                            predicted_rho, predicted_beta),
    param_type     = ifelse(grepl("rho", parameter), "Detection rate",
                            "Transmission rate"),
    wave           = factor(wave, levels = names(wave_labels),
                            labels = wave_labels),
    parameter      = factor(parameter, levels = names(param_labels),
                            labels = param_labels)
  )

fig5 <- ggplot(fig5_data,
               aes(x = f, y = observed_ratio, colour = parameter,
                   group = parameter)) +
  # Theoretical prediction lines
  geom_line(aes(y = predicted, group = param_type, linetype = param_type),
            colour = "grey40", linewidth = 0.6, inherit.aes = FALSE,
            data = fig5_data |>
              distinct(wave, param_type, f, predicted)) +
  geom_hline(yintercept = 1, colour = "grey80", linewidth = 0.4) +
  geom_line(linewidth = 0.8) +
  geom_point(size = 2.2) +
  facet_wrap(~wave, nrow = 1) +
  scale_colour_brewer(palette = "Set1", name = "Parameter",
                      labels  = parse(text = levels(fig5_data$parameter))) +
  scale_linetype_manual(
    values = c("Detection rate" = "dashed", "Transmission rate" = "dotted"),
    name   = "Theoretical\nprediction"
  ) +
  scale_x_continuous(breaks = c(0.10, 0.30, 0.50),
                     labels = c("10%", "30%", "50%")) +
  labs(
    title    = "Theoretical prediction vs observed: systematic proportional reduction",
    subtitle = paste(
      "Dashed grey = theory for detection rates (1/(1-f)).",
      "Dotted grey = theory for transmission rates (no change).",
      "Coloured = observed PMCMC ratio."
    ),
    x = "Fraction removed (f)",
    y = "Ratio: degraded posterior mean / baseline posterior mean"
  ) +
  theme_bw(base_size = 10) +
  theme(legend.position  = "right",
        panel.grid.minor = element_blank())

ggsave("results/fig5_theory_vs_observed.pdf", fig5, width = 16, height = 5)
cat("Fig 5 saved\n")

cat("\nAll figures saved to results/\n")
