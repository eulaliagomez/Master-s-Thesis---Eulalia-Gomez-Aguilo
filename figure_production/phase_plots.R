library(readr)

# import WAVE 1
wave1 <- read_csv("results/wave1_metrics.csv")
baseline1 <- read_csv("results/wave1_summary.csv")

# import WAVE 2
wave2 <- read_csv("results/wave2_metrics.csv")
baseline2 <- read_csv("results/wave2_summary.csv")

# import WAVE 3
wave3 <- read_csv("results/wave3_metrics.csv")
baseline3 <- read_csv("results/wave3_summary.csv")

# import WAVE 4
wave4 <- read_csv("results/wave4_metrics.csv")
baseline4 <- read_csv("results/wave4_summary.csv")

library(dplyr)
library(ggplot2)

# Combine baseline rows from all waves
ci_data <- bind_rows(
  baseline1 %>% filter(degradation == "baseline") %>% mutate(wave = "Wave 1"),
  baseline2 %>% filter(degradation == "baseline") %>% mutate(wave = "Wave 2"),
  baseline3 %>% filter(degradation == "baseline") %>% mutate(wave = "Wave 3"),
  baseline4 %>% filter(degradation == "baseline") %>% mutate(wave = "Wave 4")
)

# Plot
ggplot(ci_data, aes(x = wave, y = ci_width, fill = wave)) +
  geom_col() +
  facet_wrap(~ parameter) +
  scale_y_continuous(
    limits = c(0, max(ci_data$ci_width, na.rm = TRUE))
  ) +
  theme_minimal()


# =============================================================================
# wave_phase_plots.R
#
# For each wave: one plot showing the observed weekly case counts (young +
# old + total) with clearly labelled NPI phases, beta values, and Rt values.
# Saves to results/phase_plots/
#
# setwd("C:/Users/Usuario/Desktop/CBS - Internship/Analysis/DATA")
# source("wave_phase_plots.R")
# =============================================================================

library(tidyverse)

if (!dir.exists("results/phase_plots")) dir.create("results/phase_plots",
                                                   recursive = TRUE)

# ── Colours ───────────────────────────────────────────────────────────────────
# One colour per phase, consistent across all waves
phase_cols <- c(
  "0" = "#c6dbef",   # light blue  — phase 0
  "1" = "#9ecae1",   # medium blue — phase 1
  "2" = "#4292c6",   # blue        — phase 2
  "3" = "#08519c"    # dark blue   — phase 3
)
col_young <- "#2166ac"
col_old   <- "#d73027"
col_total <- "#1a7a1a"

# ── Wave definitions ──────────────────────────────────────────────────────────
wave_defs <- list(
  
  wave1 = list(
    label   = "Wave 1 — Ancestral strain",
    period  = "24 Feb \u2013 5 Jul 2020",
    file_y  = "covid19_wave1_young_NL.csv",
    file_o  = "covid19_wave1_old_NL.csv",
    phases  = data.frame(
      phase   = c(0, 1, 2),
      wk_from = c(1, 5, 12),
      wk_to   = c(4, 11, 19),
      beta    = c("beta[1]==0.520", "beta[2]==0.120", "beta[3]==0.180"),
      rt      = c("R[0]==2.12", "R[t]==0.49", "R[t]==0.73"),
      event   = c("Free transmission", "National lockdown\n23 Mar 2020",
                  "Phase-1 reopening\n11 May 2020")
    )
  ),
  
  wave2 = list(
    label   = "Wave 2 \u2014 Ancestral strain (double hump)",
    period  = "7 Sep 2020 \u2013 17 Jan 2021",
    file_y  = "covid19_wave2_young_NL.csv",
    file_o  = "covid19_wave2_old_NL.csv",
    phases  = data.frame(
      phase   = c(0, 1, 2, 3),
      wk_from = c(1, 7, 11, 14),
      wk_to   = c(6, 10, 13, 19),
      beta    = c("beta[1]==0.340", "beta[2]==0.200",
                  "beta[3]==0.430", "beta[4]==0.220"),
      rt      = c("R[t]==1.39", "R[t]==0.81",
                  "R[t]==1.75", "R[t]==0.90"),
      event   = c("Autumn resurgence", "Partial trough\nOct measures",
                  "Re-acceleration\nChristmas", "Hard lockdown\n14 Dec 2020")
    )
  ),
  
  wave3 = list(
    label   = "Wave 3 \u2014 Alpha variant (B.1.1.7)",
    period  = "1 Feb \u2013 20 Jun 2021",
    file_y  = "covid19_wave3_young_NL.csv",
    file_o  = "covid19_wave3_old_NL.csv",
    phases  = data.frame(
      phase   = c(0, 1, 2),
      wk_from = c(1, 5, 11),
      wk_to   = c(4, 10, 20),
      beta    = c("beta[1]==0.380", "beta[2]==0.360", "beta[3]==0.260"),
      rt      = c("R[0]==1.55", "R[0]==1.47", "R[t]==1.06"),
      event   = c("Lockdown\n(avondklok)", "Schools reopen\n1 Mar 2021",
                  "Easter lockdown 2\n6 Apr 2021")
    )
  ),
  
  wave4 = list(
    label   = "Wave 4 \u2014 Delta variant (B.1.617.2)",
    period  = "21 Jun \u2013 3 Oct 2021",
    file_y  = "covid19_wave4_young_NL.csv",
    file_o  = "covid19_wave4_old_NL.csv",
    phases  = data.frame(
      phase   = c(0, 1, 2),
      wk_from = c(1, 3, 8),
      wk_to   = c(2, 7, 15),
      beta    = c("beta[1]==1.426", "beta[2]==0.240", "beta[3]==0.450"),
      rt      = c("R[0]==5.81", "R[t]==0.98", "R[t]==1.83"),
      event   = c("All restrictions off\n26 Jun 2021",
                  "Emergency NPI\n9 Jul + holidays",
                  "Schools reopen\n6 Sep 2021")
    )
  )
)

# ── Helper: make one wave plot ────────────────────────────────────────────────
make_wave_plot <- function(wdef) {
  
  # Load data
  meas_y <- read_csv(wdef$file_y, show_col_types = FALSE) |>
    rename(young = reports) |> select(week, young)
  meas_o <- read_csv(wdef$file_o, show_col_types = FALSE) |>
    rename(old = reports) |> select(week, old)
  meas   <- left_join(meas_y, meas_o, by = "week") |>
    mutate(total = young + old)
  
  ph <- wdef$phases
  T  <- max(meas$week)
  
  # Background rectangles for each phase
  rects <- ph |>
    mutate(
      xmin = wk_from - 0.5,
      xmax = wk_to   + 0.5,
      fill = as.character(phase)
    )
  
  # Y position for beta labels: just above top of plot
  y_max <- max(meas$total) * 1.0
  
  # Mid-week of each phase for label placement
  ph <- ph |>
    mutate(
      wk_mid = (wk_from + wk_to) / 2,
      phase_chr = as.character(phase)
    )
  
  ggplot() +
    # Phase background shading
    geom_rect(data = rects,
              aes(xmin = xmin, xmax = xmax,
                  ymin = -Inf, ymax = Inf,
                  fill = fill),
              alpha = 0.18, inherit.aes = FALSE) +
    scale_fill_manual(values = phase_cols, guide = "none") +
    
    # Vertical phase boundaries
    geom_vline(data = ph[-1, ],
               aes(xintercept = wk_from - 0.5),
               linetype = "solid", colour = "grey50",
               linewidth = 0.5) +
    
    # Case counts
    geom_line(data = meas, aes(x = week, y = young),
              colour = col_young, linewidth = 0.9) +
    geom_point(data = meas, aes(x = week, y = young),
               colour = col_young, size = 1.8, shape = 16) +
    geom_line(data = meas, aes(x = week, y = old),
              colour = col_old, linewidth = 0.9) +
    geom_point(data = meas, aes(x = week, y = old),
               colour = col_old, size = 1.8, shape = 16) +
    geom_line(data = meas, aes(x = week, y = total),
              colour = col_total, linewidth = 1.2, linetype = "dashed") +
    
    # Beta and Rt labels inside each phase band
    geom_text(data = ph,
              aes(x = wk_mid, y = Inf,
                  label = beta),
              parse    = TRUE,
              vjust    = 2.2,
              size     = 3.8,
              fontface = "bold",
              colour   = "#154273",
              inherit.aes = FALSE) +
    geom_text(data = ph,
              aes(x = wk_mid, y = Inf,
                  label = rt),
              parse   = TRUE,
              vjust   = 4.0,
              size    = 3.4,
              colour  = "#154273",
              inherit.aes = FALSE) +
    
    # Event label at bottom of each phase
    geom_text(data = ph,
              aes(x = wk_mid, y = -Inf,
                  label = event),
              vjust   = -0.4,
              size    = 2.8,
              colour  = "grey30",
              lineheight = 0.9,
              inherit.aes = FALSE) +
    
    # Manual legend
    annotate("segment", x = -Inf, xend = -Inf, y = -Inf, yend = -Inf) +
    
    scale_x_continuous(breaks = seq(1, T, by = 1),
                       minor_breaks = NULL,
                       expand = expansion(mult = c(0.01, 0.01))) +
    scale_y_continuous(labels = scales::comma,
                       expand = expansion(mult = c(0.12, 0.15))) +
    
    labs(
      title    = wdef$label,
      subtitle = wdef$period,
      x        = "Week",
      y        = "Weekly reported cases"
    ) +
    
    # Manual colour legend for lines
    geom_line(data = data.frame(x=NA_real_, y=NA_real_, grp="Young (<60)"),
              aes(x=x, y=y, colour=grp), linewidth=0.9) +
    geom_line(data = data.frame(x=NA_real_, y=NA_real_, grp="Old (\u226560)"),
              aes(x=x, y=y, colour=grp), linewidth=0.9) +
    geom_line(data = data.frame(x=NA_real_, y=NA_real_, grp="Total"),
              aes(x=x, y=y, colour=grp), linewidth=1.2, linetype="dashed") +
    scale_colour_manual(
      values = c("Young (<60)" = col_young,
                 "Old (\u226560)"   = col_old,
                 "Total"       = col_total),
      name   = NULL
    ) +
    
    theme_bw(base_size = 12) +
    theme(
      legend.position  = "bottom",
      panel.grid.minor = element_blank(),
      panel.grid.major.x = element_blank(),
      plot.title       = element_text(face = "bold", size = 13),
      plot.subtitle    = element_text(colour = "grey40", size = 10),
      axis.text.x      = element_text(size = 9)
    )
}

# ── Generate and save all four plots ─────────────────────────────────────────
for (wn in names(wave_defs)) {
  cat("Making", wn, "...\n")
  p <- make_wave_plot(wave_defs[[wn]])
  
  # Height depends on number of phases (more phases = more bottom labels)
  h <- if (wn == "wave2") 7.0 else 6.5
  
  ggsave(file.path("results/phase_plots", paste0(wn, "_phases.pdf")),
         p, width = 12, height = h)
  ggsave(file.path("results/phase_plots", paste0(wn, "_phases.png")),
         p, width = 12, height = h, dpi = 150)
  cat("  Saved:", wn, "\n")
}

# ── Combined 4-panel figure ───────────────────────────────────────────────────
cat("Making combined 4-panel figure...\n")
library(patchwork)

p1 <- make_wave_plot(wave_defs[["wave1"]])
p2 <- make_wave_plot(wave_defs[["wave2"]])
p3 <- make_wave_plot(wave_defs[["wave3"]])
p4 <- make_wave_plot(wave_defs[["wave4"]])

combined <- (p1 + p2) / (p3 + p4) +
  plot_annotation(
    title    = "COVID-19 Netherlands: NPI phases and transmission parameters by wave",
    subtitle = "Shading = NPI phase | Blue bold = beta (MLE) | Blue normal = Rt/R0 | Grey = event label",
    theme    = theme(
      plot.title    = element_text(face = "bold", size = 14),
      plot.subtitle = element_text(colour = "grey40", size = 10)
    )
  )

ggsave("results/phase_plots/all_waves_phases.pdf",
       combined, width = 20, height = 14)
ggsave("results/phase_plots/all_waves_phases.png",
       combined, width = 20, height = 14, dpi = 150)
cat("Combined figure saved.\n")

cat("\nAll phase plots saved to results/phase_plots/\n")
cat("  wave1_phases.pdf/png\n")
cat("  wave2_phases.pdf/png\n")
cat("  wave3_phases.pdf/png\n")
cat("  wave4_phases.pdf/png\n")
cat("  all_waves_phases.pdf/png  (combined 4-panel)\n")

# including CI
setwd("C:/Users/Usuario/Desktop/CBS - Internship/Analysis/DATA")
source("wave_phase_ppc_plots.R")

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
    
    # Compute x range from all distributions combined
    all_vals <- c(prior_samp,
                  get_natural(baseline_df, p),
                  unlist(lapply(degraded_data,
                                function(d) get_natural(d, p))))
    all_vals <- all_vals[is.finite(all_vals)]
    xlo <- quantile(all_vals, 0.001)
    xhi <- quantile(all_vals, 0.999)
    
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
  
  # Extract legend as grob
  library(ggplotify)
  library(gridExtra)
  library(grid)
  
  # Combine panels in a row
  panel_row <- do.call(gridExtra::grid.arrange,
                       c(lapply(all_panels, function(p) ggplotify::as.grob(p)),
                         ncol = n_params))
  
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
 