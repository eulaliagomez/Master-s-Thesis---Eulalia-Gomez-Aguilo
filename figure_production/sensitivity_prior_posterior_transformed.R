# =============================================================================
# sensitivity_prior_posterior_transformed.R
#
# Prior vs posterior plots on the INFERENCE scale:
#   log(beta)   for transmission parameters
#   logit(rho)  for detection rates
#
# On these scales:
#   - The prior is symmetric and Gaussian
#   - Width comparisons are honest (no geometric compression)
#   - A narrower posterior genuinely means more precise inference
#
# setwd("C:/Users/Usuario/Desktop/CBS - Internship/Analysis/DATA")
# source("sensitivity_prior_posterior_transformed.R")
# =============================================================================

library(tidyverse)
library(patchwork)

SENSITIVITY_FUNCTIONS_ONLY <- TRUE
source("sensitivity_pipeline.R")
rm(SENSITIVITY_FUNCTIONS_ONLY)

if (!dir.exists("results/prior_posterior"))
  dir.create("results/prior_posterior", recursive = TRUE)

set.seed(42)
N_PRIOR <- 50000

# ── Transformation functions ──────────────────────────────────────────────────
to_inference_scale <- function(x, param) {
  if (grepl("rho", param)) log(x / (1 - x))   # logit
  else                      log(x)              # log
}

from_inference_scale <- function(x, param) {
  if (grepl("rho", param)) exp(x) / (1 + exp(x))  # inverse logit
  else                      exp(x)                  # exp
}

axis_label <- function(param) {
  if (grepl("rho", param)) "logit(\u03c1)  —  inference scale"
  else                      "log(\u03b2)  —  inference scale"
}

param_title <- function(param) {
  switch(param,
    Beta1 = "log(\u03b21)",
    Beta2 = "log(\u03b22)",
    Beta3 = "log(\u03b23)",
    Beta4 = "log(\u03b24)",
    rho_y = "logit(\u03c1y)",
    rho_o = "logit(\u03c1o)",
    param
  )
}

# ── Colours ───────────────────────────────────────────────────────────────────
col_prior    <- "grey75"
col_baseline <- "#154273"

run_cols <- c(
  "systematic_f0.10" = "#c6dbef",
  "systematic_f0.30" = "#6baed6",
  "systematic_f0.50" = "#2171b5",
  "stochastic_f0.10" = "#74c476",
  "stochastic_f0.30" = "#31a354",
  "stochastic_f0.50" = "#006d2c",
  "aggregate_b2"     = "#e17000",
  "aggregate_b4"     = "#c0392b"
)

run_labels <- c(
  "systematic_f0.10" = "Systematic 10%",
  "systematic_f0.30" = "Systematic 30%",
  "systematic_f0.50" = "Systematic 50%",
  "stochastic_f0.10" = "Stochastic 10%",
  "stochastic_f0.30" = "Stochastic 30%",
  "stochastic_f0.50" = "Stochastic 50%",
  "aggregate_b2"     = "Aggregation 2-week",
  "aggregate_b4"     = "Aggregation 4-week"
)

# ── Helper: density on inference scale ────────────────────────────────────────
make_dens_transformed <- function(x_natural, param,
                                   from = NULL, to = NULL) {
  x <- to_inference_scale(x_natural, param)
  x <- x[is.finite(x) & !is.na(x)]
  if (length(x) < 10) return(NULL)
  d <- density(x, n = 512,
               from = if (is.null(from)) quantile(x, 0.001) else from,
               to   = if (is.null(to))   quantile(x, 0.999) else to)
  data.frame(x = d$x, y = d$y)
}

# ── Helper: one parameter panel ───────────────────────────────────────────────
make_panel <- function(param, prior_mu, prior_sd,
                        baseline_df, degraded_list,
                        wave_name) {

  # Prior on inference scale: simply N(prior_mu, prior_sd)
  prior_x <- seq(prior_mu - 4*prior_sd,
                 prior_mu + 4*prior_sd, length.out = 512)
  prior_y <- dnorm(prior_x, prior_mu, prior_sd)
  prior_df <- data.frame(x = prior_x, y = prior_y)

  # Compute x range from all posteriors on inference scale
  all_transformed <- c(
    to_inference_scale(baseline_df[[param]], param),
    unlist(lapply(degraded_list,
                  function(d) to_inference_scale(d[[param]], param)))
  )
  all_transformed <- all_transformed[is.finite(all_transformed)]
  xlo <- min(quantile(all_transformed, 0.001), prior_mu - 3.5*prior_sd)
  xhi <- max(quantile(all_transformed, 0.999), prior_mu + 3.5*prior_sd)

  # Baseline posterior density on inference scale
  base_d <- make_dens_transformed(baseline_df[[param]], param,
                                   from = xlo, to = xhi)

  # Degraded posterior densities
  deg_dfs <- lapply(names(degraded_list), function(rid) {
    d <- make_dens_transformed(degraded_list[[rid]][[param]], param,
                                from = xlo, to = xhi)
    if (is.null(d)) return(NULL)
    d$run <- rid
    d
  })
  deg_df <- bind_rows(deg_dfs[!sapply(deg_dfs, is.null)])

  # Panel
  p <- ggplot() +
    # Prior — filled Gaussian
    geom_area(data = prior_df,
              aes(x = x, y = y),
              fill = col_prior, alpha = 0.55, colour = "grey50",
              linewidth = 0.4) +
    # Degraded posteriors
    geom_line(data = deg_df,
              aes(x = x, y = y, colour = run, group = run),
              linewidth = 0.7, alpha = 0.90) +
    # Baseline posterior — on top
    geom_line(data = base_d,
              aes(x = x, y = y),
              colour = col_baseline, linewidth = 1.5) +
    # Prior mean reference line
    geom_vline(xintercept = prior_mu,
               linetype = "dotted", colour = "grey50",
               linewidth = 0.5) +
    scale_colour_manual(
      values = run_cols,
      labels = run_labels,
      name   = NULL,
      drop   = FALSE
    ) +
    scale_x_continuous(
      # Add natural-scale tick labels as secondary reference
      sec.axis = sec_axis(
        transform = function(x) from_inference_scale(x, param),
        name      = "Natural scale",
        labels    = function(x) round(x, 3)
      )
    ) +
    labs(
      title = param_title(param),
      x     = axis_label(param),
      y     = "Density"
    ) +
    theme_bw(base_size = 10) +
    theme(
      legend.position  = "none",
      plot.title       = element_text(size = 11, hjust = 0.5, face = "bold"),
      panel.grid.minor = element_blank(),
      axis.text.y      = element_blank(),
      axis.ticks.y     = element_blank(),
      axis.title.x.bottom = element_text(size = 8, colour = "grey40"),
      axis.title.x.top    = element_text(size = 8, colour = "grey40"),
      axis.text.x.top     = element_text(size = 7, colour = "grey50")
    )
  p
}

# ── Helper: shared legend ─────────────────────────────────────────────────────
make_shared_legend <- function() {
  df <- data.frame(
    x   = 1, y = 1,
    run = factor(names(run_labels), levels = names(run_labels))
  )
  p_leg <- ggplot(df, aes(x = x, y = y, colour = run)) +
    geom_line(linewidth = 1.2) +
    scale_colour_manual(
      values = run_cols,
      labels = run_labels,
      name   = "Degradation scenario"
    ) +
    # Add prior and baseline to legend manually
    annotate("segment", x=0, xend=0, y=0, yend=0,
             colour=col_baseline, linewidth=1.5) +
    theme_void() +
    theme(
      legend.position = "right",
      legend.title    = element_text(size = 9, face = "bold"),
      legend.text     = element_text(size = 8),
      legend.key.width = unit(1.2, "cm")
    )
  # Extract legend
  ggpubr::get_legend(p_leg)
}

# ── Main function: one wave ────────────────────────────────────────────────────
make_wave_transformed_plot <- function(wave_name) {

  cfg        <- wave_config[[wave_name]]
  params_est <- cfg$params_est
  prior_mu   <- cfg$prior_mu
  prior_sd   <- cfg$prior_sd

  rds_base <- file.path("results", paste0("baseline_", wave_name, ".rds"))
  if (!file.exists(rds_base)) {
    cat("  SKIP — baseline not found\n"); return(NULL)
  }

  baseline_df <- as.data.frame(readRDS(rds_base))

  # Load degraded posteriors
  degraded_list <- list()
  for (rid in names(run_labels)) {
    rds <- file.path("results", sprintf("%s_%s.rds", wave_name, rid))
    if (file.exists(rds)) {
      degraded_list[[rid]] <- as.data.frame(readRDS(rds)$chain)
    }
  }

  if (length(degraded_list) == 0) {
    cat("  SKIP — no degraded runs found\n"); return(NULL)
  }

  # Build one panel per parameter
  panels <- lapply(params_est, function(p) {
    make_panel(
      param        = p,
      prior_mu     = prior_mu[[p]],
      prior_sd     = prior_sd[[p]],
      baseline_df  = baseline_df,
      degraded_list = degraded_list,
      wave_name    = wave_name
    )
  })
  names(panels) <- params_est

  # Combine with patchwork
  n <- length(params_est)
  combined <- wrap_plots(panels, nrow = 1) +
    plot_annotation(
      title    = sprintf("%s \u2014 Prior vs posterior on inference scale",
                         switch(wave_name,
                           wave1 = "Wave 1 (Ancestral)",
                           wave2 = "Wave 2 (Ancestral, double hump)",
                           wave3 = "Wave 3 (Alpha)",
                           wave4 = "Wave 4 (Delta)")),
      subtitle = paste(
        "X-axis: log scale for \u03b2, logit scale for \u03c1  (scale on which PMCMC operates).",
        "Grey = prior  |  Navy = baseline posterior  |  Blue family = under-reporting",
        "|  Orange/red = temporal aggregation.",
        "Dotted vertical = prior mean.",
        "Width on this scale is honest: narrower = genuinely more precise."
      ),
      theme = theme(
        plot.title    = element_text(face = "bold", size = 12),
        plot.subtitle = element_text(colour = "grey40", size = 8,
                                     lineheight = 1.3)
      )
    )

  combined
}

# ── Generate plots ────────────────────────────────────────────────────────────
for (wn in c("wave1","wave2","wave3","wave4")) {

  cat("\nProcessing", wn, "...\n")

  p <- make_wave_transformed_plot(wn)
  if (is.null(p)) next

  n_params <- length(wave_config[[wn]]$params_est)
  w <- max(12, n_params * 3.0)

  ggsave(
    file.path("results/prior_posterior",
              paste0(wn, "_prior_posterior_transformed.pdf")),
    p, width = w, height = 5.5
  )
  ggsave(
    file.path("results/prior_posterior",
              paste0(wn, "_prior_posterior_transformed.png")),
    p, width = w, height = 5.5, dpi = 150
  )
  cat("  Saved:", wn, "\n")
}

cat("\n=== DONE ===\n")
cat("Figures saved to results/prior_posterior/\n\n")
cat("How to read:\n")
cat("  X-axis = log(beta) or logit(rho) — the scale PMCMC works on\n")
cat("  Grey filled = prior: symmetric Gaussian, width = prior SD\n")
cat("  Navy line = baseline posterior\n")
cat("  Blue lines (light to dark) = systematic 10%, 30%, 50%\n")
cat("  Green lines = stochastic 10%, 30%, 50%\n")
cat("  Orange = 2-week aggregation\n")
cat("  Red = 4-week aggregation\n")
cat("  Dotted vertical = prior mean\n\n")
cat("Key patterns:\n")
cat("  Baseline far from dotted line = data-dominated identification\n")
cat("  Red line near grey area = monthly aggregation collapsed to prior\n")
cat("  Blue/green lines near baseline = robust to under-reporting\n")
cat("  Width here is honest: wider = genuinely more uncertain\n")
