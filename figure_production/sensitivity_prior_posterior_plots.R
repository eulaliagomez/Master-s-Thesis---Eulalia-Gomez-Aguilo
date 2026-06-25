# =============================================================================
# sensitivity_prior_posterior_plots.R
#
# For each wave: one figure with one panel per estimated parameter showing
# (on the NATURAL/untransformed scale):
#   - The prior distribution (grey filled)
#   - The baseline posterior (blue filled)
#   - All degraded posteriors (coloured lines, one per scenario)
#
# This allows direct visual assessment of:
#   (a) How much the data moved the posterior away from the prior (informativeness)
#   (b) How much degradation moves the posterior further (sensitivity)
#   (c) Whether apparent robustness is genuine or just prior domination
#
# setwd("C:/Users/Usuario/Desktop/CBS - Internship/Analysis/DATA")
# source("sensitivity_prior_posterior_plots.R")
# =============================================================================

library(tidyverse)
library(coda)

SENSITIVITY_FUNCTIONS_ONLY <- TRUE
source("sensitivity_pipeline.R")
rm(SENSITIVITY_FUNCTIONS_ONLY)

if (!dir.exists("results/prior_posterior"))
  dir.create("results/prior_posterior", recursive = TRUE)

set.seed(42)
N_PRIOR <- 50000   # draws from prior for density estimation

# ── Colours ───────────────────────────────────────────────────────────────────
col_prior    <- "grey70"
col_baseline <- "#154273"
col_sys      <- "#2166ac"   # systematic — blue family
col_sto      <- "#00a1d5"   # stochastic — cyan
col_agg2     <- "#e17000"   # aggregation 2-week
col_agg4     <- "#c0392b"   # aggregation 4-week — red (most extreme)

# ── Helper: sample from prior on NATURAL scale ───────────────────────────────
# Prior is Normal on log(beta) or logit(rho)
# We sample on transformed scale and back-transform

sample_prior_natural <- function(param, mu, sd, n = N_PRIOR) {
  z <- rnorm(n, mu, sd)
  if (grepl("rho", param)) {
    # logit-normal: inverse logit
    return(exp(z) / (1 + exp(z)))
  } else {
    # log-normal: exp
    return(exp(z))
  }
}

# ── Helper: extract posterior on natural scale from chain ─────────────────────
get_natural <- function(chain_or_df, param) {
  if (inherits(chain_or_df, "mcmc") || inherits(chain_or_df, "mcmc.list")) {
    x <- as.data.frame(chain_or_df)[[param]]
  } else {
    x <- chain_or_df[[param]]
  }
  return(x)
}

# ── Helper: density df clipped to reasonable range ────────────────────────────
make_dens <- function(x, label, from = NULL, to = NULL) {
  x <- x[is.finite(x) & !is.na(x)]
  if (length(x) < 10) return(NULL)
  d <- density(x, n = 512,
               from = if (is.null(from)) min(x) * 0.5 else from,
               to   = if (is.null(to))   max(x) * 1.5 else to)
  data.frame(x = d$x, y = d$y, label = label)
}

# ── Main function: one wave ────────────────────────────────────────────────────
make_prior_posterior_plot <- function(wave_name) {
  
  cfg        <- wave_config[[wave_name]]
  params_est <- cfg$params_est
  prior_mu   <- cfg$prior_mu
  prior_sd   <- cfg$prior_sd
  
  # Check baseline exists
  rds_base <- file.path("results", paste0("baseline_", wave_name, ".rds"))
  if (!file.exists(rds_base)) {
    cat("  SKIP — baseline not found\n"); return(NULL)
  }
  
  baseline_chain <- readRDS(rds_base)
  baseline_df    <- as.data.frame(baseline_chain)
  
  # Load all degraded posteriors that exist
  run_ids <- c(
    "systematic_f0.10", "systematic_f0.30", "systematic_f0.50",
    "stochastic_f0.10", "stochastic_f0.30", "stochastic_f0.50",
    "aggregate_b2", "aggregate_b4"
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
  
  run_cols <- c(
    "systematic_f0.10" = "#c6dbef",
    "systematic_f0.30" = "#6baed6",
    "systematic_f0.50" = "#2171b5",
    "stochastic_f0.10" = "#c7e9b4",
    "stochastic_f0.30" = "#41b6c4",
    "stochastic_f0.50" = "#0c7bdc",
    "aggregate_b2"     = "#e17000",
    "aggregate_b4"     = "#c0392b"
  )
  
  degraded_data <- list()
  for (rid in run_ids) {
    rds <- file.path("results",
                     sprintf("%s_%s.rds", wave_name, rid))
    if (file.exists(rds)) {
      res <- readRDS(rds)
      degraded_data[[rid]] <- as.data.frame(res$chain)
    }
  }
  
  if (length(degraded_data) == 0) {
    cat("  SKIP — no degraded runs found\n"); return(NULL)
  }
  
  # ── Build density data for each parameter ──────────────────────────────────
  all_panels <- list()
  
  for (p in params_est) {
    
    mu <- prior_mu[[p]]
    sd <- prior_sd[[p]]
    
    # Prior samples on natural scale
    prior_samp <- sample_prior_natural(p, mu, sd)
    
    # Compute x range from POSTERIOR samples only (not prior)
    # This zooms into the posterior region so both prior shape and
    # posterior shape are visible — if we include the prior, the wide
    # prior dominates the range and everything looks flat
    post_vals <- c(
      get_natural(baseline_df, p),
      unlist(lapply(degraded_data, function(d) get_natural(d, p)))
    )
    post_vals <- post_vals[is.finite(post_vals)]
    # Add small padding beyond the posterior range to show prior tails
    post_range <- diff(range(post_vals))
    xlo <- max(0 + 1e-6,  quantile(post_vals, 0.001) - post_range * 0.30)
    xhi <- min(           quantile(post_vals, 0.999) + post_range * 0.30,
                          if (grepl("rho", p)) 1 - 1e-6 else Inf)
    
    dens_list <- list()
    
    # Prior
    d <- make_dens(prior_samp, "Prior", from = xlo, to = xhi)
    if (!is.null(d)) dens_list[["Prior"]] <- d
    
    # Baseline posterior
    d <- make_dens(get_natural(baseline_df, p), "Baseline", from = xlo, to = xhi)
    if (!is.null(d)) dens_list[["Baseline"]] <- d
    
    # Degraded posteriors
    for (rid in names(degraded_data)) {
      d <- make_dens(get_natural(degraded_data[[rid]], p),
                     run_labels[[rid]], from = xlo, to = xhi)
      if (!is.null(d)) dens_list[[rid]] <- d
    }
    
    dens_df <- bind_rows(dens_list)
    
    # Separate prior and baseline from degraded for layering
    dens_prior    <- dens_df |> filter(label == "Prior")
    dens_baseline <- dens_df |> filter(label == "Baseline")
    dens_degrade  <- dens_df |>
      filter(!label %in% c("Prior","Baseline")) |>
      mutate(label = factor(label, levels = run_labels))
    
    # Colour lookup for degraded lines
    col_lookup <- setNames(run_cols, run_labels)
    
    # Clean parameter label
    param_title <- switch(p,
                          Beta1 = expression(beta[1]),
                          Beta2 = expression(beta[2]),
                          Beta3 = expression(beta[3]),
                          Beta4 = expression(beta[4]),
                          rho_y = expression(rho[y]),
                          rho_o = expression(rho[o]),
                          p
    )
    
    # Panel
    panel <- ggplot() +
      # Prior — filled grey
      geom_area(data = dens_prior,
                aes(x = x, y = y),
                fill = col_prior, alpha = 0.50, colour = NA) +
      # Degraded posteriors — coloured lines (drawn before baseline)
      geom_line(data = dens_degrade,
                aes(x = x, y = y, colour = label),
                linewidth = 0.65, alpha = 0.85) +
      # Baseline posterior — thick navy line on top
      geom_line(data = dens_baseline,
                aes(x = x, y = y),
                colour = col_baseline, linewidth = 1.4) +
      scale_colour_manual(
        values = col_lookup,
        name   = NULL,
        drop   = FALSE
      ) +
      labs(x = NULL, y = "Density",
           title = param_title) +
      theme_bw(base_size = 10) +
      theme(
        legend.position  = "none",
        plot.title       = element_text(size = 11, hjust = 0.5),
        panel.grid.minor = element_blank(),
        axis.text.y      = element_blank(),
        axis.ticks.y     = element_blank()
      )
    
    all_panels[[p]] <- panel
  }
  
  # ── Assemble panels ────────────────────────────────────────────────────────
  n_params <- length(params_est)
  
  # Build shared legend data
  legend_df <- data.frame(
    x = 1, y = 1,
    label = factor(run_labels, levels = run_labels)
  )
  legend_plot <- ggplot(legend_df, aes(x = x, y = y, colour = label)) +
    geom_line(linewidth = 1) +
    scale_colour_manual(
      values = run_cols,
      name   = "Degradation scenario",
      labels = run_labels
    ) +
    # Add manual entries for prior and baseline
    guides(colour = guide_legend(
      override.aes = list(linewidth = 1.2)
    )) +
    theme_void() +
    theme(legend.position = "right",
          legend.title    = element_text(size = 9, face = "bold"),
          legend.text     = element_text(size = 8))
  
  # Full figure using patchwork
  library(patchwork)
  combined <- wrap_plots(all_panels, nrow = 1) +
    plot_annotation(
      title    = sprintf("%s: Prior vs posterior on natural scale",
                         toupper(wave_name)),
      subtitle = paste(
        "Grey filled = prior.",
        "Thick navy line = baseline posterior.",
        "Coloured lines = degraded posteriors",
        "(blue family = under-reporting; orange/red = aggregation)."
      ),
      theme = theme(
        plot.title    = element_text(face = "bold", size = 12),
        plot.subtitle = element_text(colour = "grey40", size = 9)
      )
    )
  
  combined
}

# ── Custom legend panel ────────────────────────────────────────────────────────
make_legend_panel <- function() {
  legend_data <- data.frame(
    x     = rep(1:3, 8),
    y     = rep(seq(8, 1, by = -1), each = 3),
    group = rep(c(
      "Systematic 10%","Systematic 30%","Systematic 50%",
      "Stochastic 10%","Stochastic 30%","Stochastic 50%",
      "Aggregation 2-week","Aggregation 4-week"
    ), each = 3)
  )
  cols <- c(
    "Systematic 10%" = "#c6dbef", "Systematic 30%" = "#6baed6",
    "Systematic 50%" = "#2171b5", "Stochastic 10%" = "#c7e9b4",
    "Stochastic 30%" = "#41b6c4", "Stochastic 50%" = "#0c7bdc",
    "Aggregation 2-week" = "#e17000", "Aggregation 4-week" = "#c0392b"
  )
  ggplot(legend_data, aes(x=x, y=y, colour=group, group=group)) +
    geom_line(linewidth=1.2) +
    scale_colour_manual(values=cols, name="Degradation scenario") +
    theme_void() +
    theme(legend.position="right",
          legend.title=element_text(size=9, face="bold"),
          legend.text=element_text(size=8))
}

# ── Generate plots ────────────────────────────────────────────────────────────
wave_titles <- c(
  wave1 = "Wave 1 (Ancestral)",
  wave2 = "Wave 2 (Ancestral, double hump)",
  wave3 = "Wave 3 (Alpha)",
  wave4 = "Wave 4 (Delta)"
)

for (wn in c("wave1","wave2","wave3","wave4")) {
  cat("\nProcessing", wn, "...\n")
  
  rds_base <- file.path("results", paste0("baseline_", wn, ".rds"))
  if (!file.exists(rds_base)) {
    cat("  SKIP — baseline not found\n"); next
  }
  
  p <- make_prior_posterior_plot(wn)
  if (is.null(p)) next
  
  n_params <- length(wave_config[[wn]]$params_est)
  w <- max(10, n_params * 2.8)
  
  ggsave(file.path("results/prior_posterior",
                   paste0(wn, "_prior_posterior.pdf")),
         p, width = w, height = 5)
  ggsave(file.path("results/prior_posterior",
                   paste0(wn, "_prior_posterior.png")),
         p, width = w, height = 5, dpi = 150)
  cat("  Saved:", wn, "\n")
}

cat("\n=== DONE ===\n")
cat("Figures saved to results/prior_posterior/\n")
cat("\nHow to read each figure:\n")
cat("  Grey filled area  = prior distribution (what we assumed before seeing data)\n")
cat("  Thick navy line   = baseline posterior (what data + prior together imply)\n")
cat("  Light blue lines  = systematic under-reporting posteriors (10%, 30%, 50%)\n")
cat("  Cyan lines        = stochastic under-reporting posteriors\n")
cat("  Orange line       = 2-week aggregation posterior\n")
cat("  Red line          = 4-week aggregation posterior\n")
cat("\nKey patterns to look for:\n")
cat("  Baseline far from prior  = data-dominated = genuine parameter identification\n")
cat("  Baseline near prior      = prior-dominated = robustness may be artificial\n")
cat("  Degraded lines near baseline = parameter is robust to that degradation\n")
cat("  Degraded lines near prior    = degradation destroyed data information\n")