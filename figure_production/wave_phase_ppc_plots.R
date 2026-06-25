# =============================================================================
# wave_phase_ppc_plots.R
#
# For each wave: observed total cases + baseline posterior predictive 95% CI,
# with NPI phase bands and beta/Rt labels clearly shown.
#
# Requires:
#   - results/baseline_wave{1-4}.rds  (saved from wave scripts)
#   - sensitivity_pipeline.R          (for build_pomp, wave_config)
#   - All wave CSV data files
#
# setwd("C:/Users/Usuario/Desktop/CBS - Internship/Analysis/DATA")
# source("wave_phase_ppc_plots.R")
# =============================================================================

library(tidyverse)
library(pomp)
library(patchwork)

SENSITIVITY_FUNCTIONS_ONLY <- TRUE
source("sensitivity_pipeline.R")
rm(SENSITIVITY_FUNCTIONS_ONLY)

if (!dir.exists("results/phase_plots")) dir.create("results/phase_plots",
                                                    recursive = TRUE)

set.seed(42)

# ── Colours ───────────────────────────────────────────────────────────────────
phase_fills <- c("0" = "#ddeef6",
                 "1" = "#b8d9ed",
                 "2" = "#7db8db",
                 "3" = "#3a7fbf")
col_total   <- "#154273"
col_obs     <- "black"

# ── Wave phase definitions ────────────────────────────────────────────────────
wave_phases <- list(

  wave1 = data.frame(
    phase   = c(0, 1, 2),
    wk_from = c(1, 5, 12),
    wk_to   = c(4, 11, 19),
    beta    = c("beta[1]==0.520", "beta[2]==0.120", "beta[3]==0.180"),
    rt      = c("R[0]==2.12",     "R[t]==0.49",     "R[t]==0.73"),
    event   = c("Free\ntransmission", "National\nlockdown", "Phase-1\nreopening"),
    stringsAsFactors = FALSE
  ),

  wave2 = data.frame(
    phase   = c(0, 1, 2, 3),
    wk_from = c(1, 7, 11, 14),
    wk_to   = c(6, 10, 13, 19),
    beta    = c("beta[1]==0.340", "beta[2]==0.200",
                "beta[3]==0.430", "beta[4]==0.220"),
    rt      = c("R[t]==1.39", "R[t]==0.81",
                "R[t]==1.75", "R[t]==0.90"),
    event   = c("Autumn\nresurgence", "Partial\ntrough",
                "Re-\nacceleration", "Hard\nlockdown"),
    stringsAsFactors = FALSE
  ),

  wave3 = data.frame(
    phase   = c(0, 1, 2),
    wk_from = c(1, 5, 11),
    wk_to   = c(4, 10, 20),
    beta    = c("beta[1]==0.380", "beta[2]==0.360", "beta[3]==0.260"),
    rt      = c("R[0]==1.55",     "R[0]==1.47",     "R[t]==1.06"),
    event   = c("Lockdown\n(avondklok)", "Schools\nreopen", "Easter\nlockdown 2"),
    stringsAsFactors = FALSE
  ),

  wave4 = data.frame(
    phase   = c(0, 1, 2),
    wk_from = c(1, 3, 8),
    wk_to   = c(2, 7, 15),
    beta    = c("beta[1]==1.426", "beta[2]==0.240", "beta[3]==0.450"),
    rt      = c("R[0]==5.81",     "R[t]==0.98",     "R[t]==1.83"),
    event   = c("Restrictions\nlifted", "Emergency\nNPI + holidays",
                "Schools\nreopen"),
    stringsAsFactors = FALSE
  )
)

wave_labels <- c(
  wave1 = "Wave 1 \u2014 Ancestral  |  24 Feb \u2013 5 Jul 2020",
  wave2 = "Wave 2 \u2014 Ancestral (double hump)  |  7 Sep 2020 \u2013 17 Jan 2021",
  wave3 = "Wave 3 \u2014 Alpha (B.1.1.7)  |  1 Feb \u2013 20 Jun 2021",
  wave4 = "Wave 4 \u2014 Delta (B.1.617.2)  |  21 Jun \u2013 3 Oct 2021"
)

# ── Helper: build PPC for one wave ────────────────────────────────────────────
make_ppc_total <- function(wave_name, n_draws = 500) {

  cfg   <- wave_config[[wave_name]]
  chain <- readRDS(file.path("results",
                             paste0("baseline_", wave_name, ".rds")))
  post_df <- as.data.frame(chain)

  # Load clean data
  meas_y <- read_csv(cfg$file_y, show_col_types = FALSE) |>
    rename(reports_young = reports) |> select(week, reports_young)
  meas_o <- read_csv(cfg$file_o, show_col_types = FALSE) |>
    rename(reports_old = reports) |> select(week, reports_old)
  meas   <- left_join(meas_y, meas_o, by = "week") |>
    mutate(total = reports_young + reports_old)

  # Build pomp object
  po <- build_pomp(cfg, meas)
  po <- pomp(po, dprior = make_prior(cfg), paramnames = names(coef(po)))

  # Simulate from posterior
  draw_idx <- sample(nrow(post_df), min(n_draws, nrow(post_df)))
  sims <- lapply(draw_idx, function(i) {
    th                 <- coef(po)
    th[cfg$params_est] <- as.numeric(post_df[i, cfg$params_est])
    simulate(po, params = th, nsim = 1, format = "data.frame") |>
      mutate(total = reports_young + reports_old)
  })

  # Summarise
  ppc_band <- bind_rows(sims) |>
    group_by(week) |>
    summarise(
      lo  = quantile(total, 0.025, na.rm = TRUE),
      hi  = quantile(total, 0.975, na.rm = TRUE),
      med = median(total, na.rm = TRUE),
      .groups = "drop"
    )

  list(ppc = ppc_band, obs = meas)
}

# ── Helper: make one wave plot ────────────────────────────────────────────────
make_phase_ppc_plot <- function(wave_name) {

  cat("  Building PPC for", wave_name, "...\n")
  dat <- make_ppc_total(wave_name)
  ppc <- dat$ppc
  obs <- dat$obs
  ph  <- wave_phases[[wave_name]]
  T   <- max(obs$week)

  # Phase rectangles
  rects <- ph |>
    mutate(xmin = wk_from - 0.5,
           xmax = wk_to   + 0.5,
           fill = as.character(phase))

  # Mid-week for labels
  ph <- ph |> mutate(wk_mid = (wk_from + wk_to) / 2,
                     phase_chr = as.character(phase))

  ggplot() +

    # Phase background
    geom_rect(data = rects,
              aes(xmin = xmin, xmax = xmax,
                  ymin = -Inf, ymax = Inf, fill = fill),
              alpha = 0.20, inherit.aes = FALSE) +
    scale_fill_manual(values = phase_fills, guide = "none") +

    # Phase boundaries
    geom_vline(data    = ph[-1, ],
               aes(xintercept = wk_from - 0.5),
               linetype = "solid", colour = "grey50", linewidth = 0.5) +

    # Posterior 95% CI band
    geom_ribbon(data = ppc,
                aes(x = week, ymin = lo, ymax = hi),
                fill = col_total, alpha = 0.20) +

    # Posterior median line
    geom_line(data = ppc,
              aes(x = week, y = med),
              colour = col_total, linewidth = 1.1) +

    # Observed total cases
    geom_point(data = obs,
               aes(x = week, y = total),
               colour = col_obs, size = 2.5, shape = 16) +
    geom_line(data = obs,
              aes(x = week, y = total),
              colour = col_obs, linewidth = 0.5, linetype = "dotted") +

    # Beta labels at top
    geom_text(data = ph,
              aes(x = wk_mid, y = Inf, label = beta),
              parse = TRUE, vjust = 2.1, size = 3.5,
              fontface = "bold", colour = "#154273",
              inherit.aes = FALSE) +

    # Rt labels just below beta
    geom_text(data = ph,
              aes(x = wk_mid, y = Inf, label = rt),
              parse = TRUE, vjust = 3.9, size = 3.0,
              colour = "#154273", inherit.aes = FALSE) +

    # Event labels at bottom
    geom_text(data = ph,
              aes(x = wk_mid, y = -Inf, label = event),
              vjust = -0.3, size = 2.7, colour = "grey30",
              lineheight = 0.85, inherit.aes = FALSE) +

    scale_x_continuous(breaks  = seq(1, T, by = 1),
                       minor_breaks = NULL,
                       expand  = expansion(mult = c(0.01, 0.01))) +
    scale_y_continuous(labels  = scales::comma,
                       expand  = expansion(mult = c(0.12, 0.16))) +

    labs(
      title    = wave_labels[[wave_name]],
      subtitle = paste("Shaded band = baseline 95% posterior predictive interval  |",
                       "Solid line = posterior median  |",
                       "Black dots = observed total cases"),
      x = "Week",
      y = "Weekly reported cases (total)"
    ) +

    theme_bw(base_size = 12) +
    theme(
      panel.grid.minor   = element_blank(),
      panel.grid.major.x = element_blank(),
      plot.title         = element_text(face = "bold", size = 12),
      plot.subtitle      = element_text(colour = "grey40", size = 9),
      axis.text.x        = element_text(size = 9)
    )
}

# ── Generate individual plots ─────────────────────────────────────────────────
plots <- list()
for (wn in names(wave_phases)) {
  cat("Processing", wn, "...\n")

  rds <- file.path("results", paste0("baseline_", wn, ".rds"))
  if (!file.exists(rds)) {
    cat("  SKIP — baseline not found:", rds, "\n\n")
    next
  }

  p <- make_phase_ppc_plot(wn)
  plots[[wn]] <- p

  h <- ifelse(wn == "wave2", 7.0, 6.5)
  ggsave(file.path("results/phase_plots",
                   paste0(wn, "_phases_ppc.pdf")), p, width = 12, height = h)
  ggsave(file.path("results/phase_plots",
                   paste0(wn, "_phases_ppc.png")), p, width = 12, height = h,
         dpi = 150)
  cat("  Saved:", wn, "\n\n")
}

# ── Combined 4-panel figure ───────────────────────────────────────────────────
if (length(plots) == 4) {
  cat("Making combined 4-panel figure...\n")
  combined <- (plots[["wave1"]] + plots[["wave2"]]) /
              (plots[["wave3"]] + plots[["wave4"]]) +
    plot_annotation(
      title    = "COVID-19 Netherlands: NPI phase structure and model fit (total cases)",
      subtitle = paste(
        "Shaded bands = baseline 95% posterior predictive interval.",
        "Solid lines = posterior median.",
        "Black dots = observed total cases.",
        "Beta and Rt values are MLE starting values."
      ),
      theme = theme(
        plot.title    = element_text(face = "bold", size = 14),
        plot.subtitle = element_text(colour = "grey40", size = 9)
      )
    )

  ggsave("results/phase_plots/all_waves_phases_ppc.pdf",
         combined, width = 20, height = 14)
  ggsave("results/phase_plots/all_waves_phases_ppc.png",
         combined, width = 20, height = 14, dpi = 150)
  cat("Combined figure saved.\n")
} else {
  cat("Not all waves available — combined figure skipped.\n")
}

cat("\nAll plots saved to results/phase_plots/\n")
