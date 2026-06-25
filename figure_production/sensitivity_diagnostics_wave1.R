# =============================================================================
# sensitivity_diagnostics_wave1.R
#
# Produces posterior diagnostics for every Wave 1 sensitivity run:
#   (1) Posterior marginal densities — with baseline overlaid in red
#   (2) Pairwise posterior scatter (correlations)
#   (3) Posterior predictive check — young / old / total panels
#
# Plots are saved as PNG and PDF in results/diagnostics/
# One subfolder per run: results/diagnostics/systematic_f0.10/ etc.
#
# USAGE:
#   setwd("C:/Users/Usuario/Desktop/CBS - Internship/Analysis/DATA")
#   source("sensitivity_diagnostics_wave1.R")
#
# REQUIRES: sensitivity_pipeline.R sourced first (for build_pomp, wave_config)
# =============================================================================

library(tidyverse)
library(pomp)
library(coda)
library(GGally)

# ── 0. Setup ──────────────────────────────────────────────────────────────────
SENSITIVITY_FUNCTIONS_ONLY <- TRUE
source("sensitivity_pipeline.R")
rm(SENSITIVITY_FUNCTIONS_ONLY)

if (!dir.exists("results/diagnostics")) dir.create("results/diagnostics")

cfg       <- wave_config[["wave1"]]
params_est <- cfg$params_est
dom_eigen  <- 4.074

series_cols <- c("Young (<60)" = "#2166ac",
                 "Old (\u226560)"   = "#d73027",
                 "Total"       = "#1a7a1a")

# ── 1. Load clean data and baseline ───────────────────────────────────────────
meas_y <- read_csv(cfg$file_y, show_col_types = FALSE) |>
  rename(reports_young = reports) |> select(week, reports_young)
meas_o <- read_csv(cfg$file_o, show_col_types = FALSE) |>
  rename(reports_old = reports) |> select(week, reports_old)
meas_clean <- left_join(meas_y, meas_o, by = "week")

baseline_chain <- readRDS("results/baseline_wave1.rds")
baseline_df    <- as.data.frame(baseline_chain)

cat("Baseline loaded. ESS:", round(effectiveSize(baseline_chain)), "\n\n")

# ── 2. Define all runs ────────────────────────────────────────────────────────
runs <- list(
  list(id="systematic_f0.10",  label="Systematic f=10%",
       type="systematic", level="0.1",  f=0.10, block=NA),
  list(id="systematic_f0.30",  label="Systematic f=30%",
       type="systematic", level="0.3",  f=0.30, block=NA),
  list(id="systematic_f0.50",  label="Systematic f=50%",
       type="systematic", level="0.5",  f=0.50, block=NA),
  list(id="stochastic_f0.10",  label="Stochastic f=10%",
       type="stochastic", level="0.1",  f=0.10, block=NA),
  list(id="stochastic_f0.30",  label="Stochastic f=30%",
       type="stochastic", level="0.3",  f=0.30, block=NA),
  list(id="stochastic_f0.50",  label="Stochastic f=50%",
       type="stochastic", level="0.5",  f=0.50, block=NA),
  list(id="aggregate_b2",      label="Aggregation 2-week",
       type="aggregate",  level="2",    f=NA,   block=2),
  list(id="aggregate_b4",      label="Aggregation 4-week",
       type="aggregate",  level="4",    f=NA,   block=4)
)

# ── 3. Helper: get degraded data for a run ────────────────────────────────────
get_degraded_data <- function(run) {
  csv_path <- file.path("results/degraded_data",
                        sprintf("wave1_%s.csv", run$id))
  if (file.exists(csv_path)) {
    df <- read_csv(csv_path, show_col_types = FALSE)
    # ensure column names match
    if (!"reports_young" %in% colnames(df)) {
      df <- df |> rename(reports_young = reports_young,
                         reports_old   = reports_old)
    }
    return(df)
  }
  # Recreate if CSV missing
  if (run$type == "systematic") {
    return(degrade_systematic(meas_clean, run$f))
  } else if (run$type == "stochastic") {
    deg_seed <- as.integer(paste0(1, round(run$f * 100)))
    return(degrade_stochastic(meas_clean, run$f, seed = deg_seed))
  } else {
    return(degrade_aggregate(meas_clean, run$block)$data)
  }
}

# ── 4. Helper: PPC for one run ────────────────────────────────────────────────
make_ppc <- function(po, post_df_run, meas_df, run_label, baseline_df_arg,
                     n_draws = 500) {

  set.seed(42)
  draw_idx <- sample(nrow(post_df_run), min(n_draws, nrow(post_df_run)))

  # Degraded-data PPC
  sims_deg <- lapply(draw_idx, function(i) {
    th             <- coef(po)
    th[params_est] <- as.numeric(post_df_run[i, params_est])
    simulate(po, params = th, nsim = 1, format = "data.frame") |>
      mutate(draw = i, reports_total = reports_young + reports_old,
             source = "Degraded posterior")
  })

  # Baseline PPC (same number of draws, on same data for comparison)
  po_base <- build_pomp(cfg, meas_clean)
  po_base <- pomp(po_base, dprior = make_prior(cfg),
                  paramnames = names(coef(po_base)))
  draw_idx_b <- sample(nrow(baseline_df_arg),
                       min(n_draws, nrow(baseline_df_arg)))
  sims_base <- lapply(draw_idx_b, function(i) {
    th             <- coef(po_base)
    th[params_est] <- as.numeric(baseline_df_arg[i, params_est])
    simulate(po_base, params = th, nsim = 1, format = "data.frame") |>
      mutate(draw = i, reports_total = reports_young + reports_old,
             source = "Baseline posterior")
  })

  all_sims <- bind_rows(c(sims_deg, sims_base))

  # Summarise bands
  bands <- all_sims |>
    pivot_longer(c(reports_young, reports_old, reports_total),
                 names_to = "series", values_to = "cases") |>
    mutate(series = recode(series,
      "reports_young" = "Young (<60)",
      "reports_old"   = "Old (\u226560)",
      "reports_total" = "Total"),
      series = factor(series,
                      levels = c("Young (<60)", "Old (\u226560)", "Total"))) |>
    group_by(week, series, source) |>
    summarise(lo  = quantile(cases, 0.025, na.rm = TRUE),
              hi  = quantile(cases, 0.975, na.rm = TRUE),
              med = median(cases, na.rm = TRUE),
              .groups = "drop")

  # Observed data (always the full clean data for reference)
  obs <- meas_clean |>
    mutate(reports_total = reports_young + reports_old) |>
    pivot_longer(c(reports_young, reports_old, reports_total),
                 names_to = "series", values_to = "cases") |>
    mutate(series = recode(series,
      "reports_young" = "Young (<60)",
      "reports_old"   = "Old (\u226560)",
      "reports_total" = "Total"),
      series = factor(series,
                      levels = c("Young (<60)", "Old (\u226560)", "Total")))

  ggplot() +
    geom_ribbon(data = bands |> filter(source == "Baseline posterior"),
                aes(x = week, ymin = lo, ymax = hi),
                fill = "#2166ac", alpha = 0.15) +
    geom_line(data = bands |> filter(source == "Baseline posterior"),
              aes(x = week, y = med),
              colour = "#2166ac", linewidth = 0.8, linetype = "dashed") +
    geom_ribbon(data = bands |> filter(source == "Degraded posterior"),
                aes(x = week, ymin = lo, ymax = hi),
                fill = "#d73027", alpha = 0.20) +
    geom_line(data = bands |> filter(source == "Degraded posterior"),
              aes(x = week, y = med),
              colour = "#d73027", linewidth = 1.0) +
    geom_point(data = obs,
               aes(x = week, y = cases),
               colour = "black", size = 2, shape = 16) +
    facet_wrap(~series, scales = "free_y", nrow = 3) +
    scale_y_continuous(labels = scales::comma) +
    labs(
      title    = sprintf("Wave 1 PPC: %s", run_label),
      subtitle = paste(
        "Red = degraded posterior (solid line + dark band).",
        "Blue = baseline posterior (dashed line + light band).",
        "Black dots = observed data."
      ),
      x = "Week", y = "Weekly reported cases"
    ) +
    theme_bw(base_size = 11)
}

# ── 5. Helper: posterior density plot with baseline overlay ───────────────────
make_density <- function(post_df_run, run_label) {

  deg_long  <- post_df_run |>
    select(all_of(params_est)) |>
    pivot_longer(everything(), names_to = "parameter", values_to = "value") |>
    mutate(source = "Degraded")

  base_long <- baseline_df |>
    select(all_of(params_est)) |>
    pivot_longer(everything(), names_to = "parameter", values_to = "value") |>
    mutate(source = "Baseline")

  bind_rows(deg_long, base_long) |>
    mutate(parameter = factor(parameter, levels = params_est)) |>
    ggplot(aes(x = value, fill = source, colour = source)) +
    geom_density(alpha = 0.35, linewidth = 0.8) +
    facet_wrap(~parameter, scales = "free", nrow = 2,
               labeller = as_labeller(c(
                 Beta1 = "beta[1]", Beta2 = "beta[2]", Beta3 = "beta[3]",
                 rho_y = "rho[y]",  rho_o = "rho[o]"
               ), label_parsed)) +
    scale_fill_manual(values   = c("Baseline" = "#2166ac",
                                   "Degraded" = "#d73027"),
                      name = NULL) +
    scale_colour_manual(values = c("Baseline" = "#2166ac",
                                   "Degraded" = "#d73027"),
                        name = NULL) +
    labs(
      title    = sprintf("Wave 1 posterior densities: %s", run_label),
      subtitle = "Red = degraded posterior  |  Blue = baseline posterior",
      x = "Parameter value", y = "Density"
    ) +
    theme_bw(base_size = 11) +
    theme(legend.position = "bottom",
          strip.text      = element_text(size = 9))
}

# ── 6. Helper: pairwise scatter ───────────────────────────────────────────────
make_pairs <- function(post_df_run, run_label) {
  tryCatch(
    ggpairs(
      post_df_run |> select(all_of(params_est)),
      lower = list(continuous = wrap("points", alpha = 0.04, size = 0.3,
                                     colour = "#d73027")),
      upper = list(continuous = wrap("cor", size = 3.0)),
      diag  = list(continuous = wrap("densityDiag", colour = "#d73027",
                                     fill = "#d73027", alpha = 0.4))
    ) +
      theme_bw(base_size = 9) +
      labs(title = sprintf("Wave 1 pairwise posterior: %s", run_label)),
    error = function(e) {
      ggplot() +
        annotate("text", x = 0.5, y = 0.5,
                 label = paste("Pairs plot failed:\n", e$message),
                 size = 4) +
        theme_void()
    }
  )
}

# ── 7. MAIN LOOP — produce diagnostics for every run ─────────────────────────
for (run in runs) {

  cat(strrep("-", 60), "\n")
  cat("Processing:", run$label, "\n")

  # Load posterior chain
  rds_path <- file.path("results",
                        sprintf("wave1_%s.rds", run$id))
  if (!file.exists(rds_path)) {
    cat("  SKIP — rds not found:", rds_path, "\n\n")
    next
  }

  res        <- readRDS(rds_path)
  post_df_run <- as.data.frame(res$chain)

  cat("  ESS:", paste(round(res$ess), collapse = " "), "\n")
  cat("  Draws:", nrow(post_df_run), "\n")

  # Create output directory for this run
  out_dir <- file.path("results/diagnostics", run$id)
  if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

  # Get degraded data and rebuild pomp object
  meas_deg  <- get_degraded_data(run)

  if (run$type == "aggregate") {
    agg       <- degrade_aggregate(meas_clean, run$block)
    po        <- build_pomp(cfg, agg$data, obs_times = agg$obs_times)
  } else {
    po        <- build_pomp(cfg, meas_deg)
  }
  po <- pomp(po, dprior = make_prior(cfg), paramnames = names(coef(po)))

  # ── Plot 1: Posterior densities ──────────────────────────────────────────
  cat("  Making density plot...\n")
  p_dens <- make_density(post_df_run, run$label)
  ggsave(file.path(out_dir, "posterior_densities.pdf"),
         p_dens, width = 12, height = 7)
  ggsave(file.path(out_dir, "posterior_densities.png"),
         p_dens, width = 12, height = 7, dpi = 150)

  # ── Plot 2: Pairwise scatter ─────────────────────────────────────────────
  cat("  Making pairs plot...\n")
  p_pairs <- make_pairs(post_df_run, run$label)
  ggsave(file.path(out_dir, "pairwise_scatter.pdf"),
         p_pairs, width = 10, height = 10)
  ggsave(file.path(out_dir, "pairwise_scatter.png"),
         p_pairs, width = 10, height = 10, dpi = 150)

  # ── Plot 3: PPC ──────────────────────────────────────────────────────────
  cat("  Making PPC...\n")
  p_ppc <- make_ppc(po, post_df_run, meas_deg, run$label, baseline_df)
  ggsave(file.path(out_dir, "ppc.pdf"),
         p_ppc, width = 10, height = 10)
  ggsave(file.path(out_dir, "ppc.png"),
         p_ppc, width = 10, height = 10, dpi = 150)

  cat("  Saved to:", out_dir, "\n\n")
}

# ── 8. Baseline diagnostics (for completeness) ────────────────────────────────
cat(strrep("-", 60), "\n")
cat("Processing: Baseline\n")

out_dir_base <- "results/diagnostics/baseline"
if (!dir.exists(out_dir_base)) dir.create(out_dir_base, recursive = TRUE)

po_base <- build_pomp(cfg, meas_clean)
po_base <- pomp(po_base, dprior = make_prior(cfg),
                paramnames = names(coef(po_base)))

# Baseline density (self-overlay for consistency)
p_base_dens <- baseline_df |>
  select(all_of(params_est)) |>
  pivot_longer(everything(), names_to = "parameter", values_to = "value") |>
  mutate(parameter = factor(parameter, levels = params_est)) |>
  ggplot(aes(x = value)) +
  geom_histogram(aes(y = after_stat(density)), bins = 50,
                 fill = "#2166ac", alpha = 0.55) +
  geom_density(colour = "black", linewidth = 0.8) +
  facet_wrap(~parameter, scales = "free", nrow = 2,
             labeller = as_labeller(c(
               Beta1 = "beta[1]", Beta2 = "beta[2]", Beta3 = "beta[3]",
               rho_y = "rho[y]",  rho_o = "rho[o]"
             ), label_parsed)) +
  labs(title    = "Wave 1 posterior densities: Baseline (clean data)",
       x = "Parameter value", y = "Density") +
  theme_bw(base_size = 11)

ggsave(file.path(out_dir_base, "posterior_densities.pdf"),
       p_base_dens, width = 12, height = 7)
ggsave(file.path(out_dir_base, "posterior_densities.png"),
       p_base_dens, width = 12, height = 7, dpi = 150)

# Baseline pairs
p_base_pairs <- make_pairs(baseline_df, "Baseline (clean data)")
ggsave(file.path(out_dir_base, "pairwise_scatter.pdf"),
       p_base_pairs, width = 10, height = 10)
ggsave(file.path(out_dir_base, "pairwise_scatter.png"),
       p_base_pairs, width = 10, height = 10, dpi = 150)

# Baseline PPC
set.seed(42)
draw_idx_b <- sample(nrow(baseline_df), 500)
sims_base_only <- lapply(draw_idx_b, function(i) {
  th             <- coef(po_base)
  th[params_est] <- as.numeric(baseline_df[i, params_est])
  simulate(po_base, params = th, nsim = 1, format = "data.frame") |>
    mutate(draw = i, reports_total = reports_young + reports_old)
})

pp_long_base <- bind_rows(sims_base_only) |>
  pivot_longer(c(reports_young, reports_old, reports_total),
               names_to = "series", values_to = "cases") |>
  mutate(series = recode(series,
    "reports_young" = "Young (<60)",
    "reports_old"   = "Old (\u226560)",
    "reports_total" = "Total"),
    series = factor(series,
                    levels = c("Young (<60)", "Old (\u226560)", "Total")))

pp_band_base <- pp_long_base |>
  group_by(week, series) |>
  summarise(lo  = quantile(cases, 0.025, na.rm = TRUE),
            hi  = quantile(cases, 0.975, na.rm = TRUE),
            med = median(cases, na.rm = TRUE),
            .groups = "drop")

obs_long_base <- meas_clean |>
  mutate(reports_total = reports_young + reports_old) |>
  pivot_longer(c(reports_young, reports_old, reports_total),
               names_to = "series", values_to = "cases") |>
  mutate(series = recode(series,
    "reports_young" = "Young (<60)",
    "reports_old"   = "Old (\u226560)",
    "reports_total" = "Total"),
    series = factor(series,
                    levels = c("Young (<60)", "Old (\u226560)", "Total")))

p_base_ppc <- ggplot() +
  geom_ribbon(data = pp_band_base,
              aes(x = week, ymin = lo, ymax = hi, fill = series),
              alpha = 0.25) +
  geom_line(data = pp_band_base,
            aes(x = week, y = med, colour = series),
            linewidth = 1) +
  geom_point(data = obs_long_base,
             aes(x = week, y = cases, colour = series),
             size = 2, shape = 16) +
  facet_wrap(~series, scales = "free_y", nrow = 3) +
  scale_colour_manual(values = series_cols, guide = "none") +
  scale_fill_manual(values   = series_cols, guide = "none") +
  scale_y_continuous(labels  = scales::comma) +
  labs(title    = "Wave 1 PPC: Baseline (clean data)",
       subtitle = "Shaded=95% interval; line=median; dots=observed",
       x = "Week", y = "Weekly reported cases") +
  theme_bw(base_size = 11)

ggsave(file.path(out_dir_base, "ppc.pdf"),
       p_base_ppc, width = 10, height = 10)
ggsave(file.path(out_dir_base, "ppc.png"),
       p_base_ppc, width = 10, height = 10, dpi = 150)

cat("  Baseline diagnostics saved to:", out_dir_base, "\n\n")

# ── 9. Summary ────────────────────────────────────────────────────────────────
cat(strrep("=", 60), "\n")
cat("ALL DIAGNOSTICS COMPLETE\n")
cat("Output structure:\n")
cat("  results/diagnostics/\n")
cat("    baseline/\n")
for (run in runs) cat(sprintf("    %s/\n", run$id))
cat("\n  Each folder contains:\n")
cat("    posterior_densities.pdf/.png\n")
cat("    pairwise_scatter.pdf/.png\n")
cat("    ppc.pdf/.png\n")
