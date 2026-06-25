# =============================================================================
# run_wave2_sensitivity.R
#
# Sensitivity analysis for Wave 2 only.
# Identical structure to run_wave1_sensitivity.R.
#
# USAGE:
#   1. Run covid19_age_wave2_v2.R and save baseline:
#      saveRDS(prod_chain, "results/baseline_wave2.rds")
#   2. source("run_wave2_sensitivity.R")
# =============================================================================

SENSITIVITY_FUNCTIONS_ONLY <- TRUE
source("sensitivity_pipeline.R")
rm(SENSITIVITY_FUNCTIONS_ONLY)

if (!dir.exists("results"))               dir.create("results")
if (!dir.exists("results/degraded_data")) dir.create("results/degraded_data")
if (!dir.exists("results/wave2"))         dir.create("results/wave2")

cfg        <- wave_config[["wave2"]]
wave_name  <- "wave2"
dom_eigen  <- 4.074
f_levels   <- c(0.10, 0.30, 0.50)
agg_blocks <- c(2, 4)

cat("\n", strrep("=", 70), "\n")
cat("SENSITIVITY ANALYSIS — WAVE 2 ONLY\n")
cat(cfg$label, "\n")
cat(strrep("=", 70), "\n")

# ── Load clean data ───────────────────────────────────────────────────────────
meas_y <- read_csv(cfg$file_y, show_col_types = FALSE) |>
  rename(reports_young = reports) |> select(week, reports_young)
meas_o <- read_csv(cfg$file_o, show_col_types = FALSE) |>
  rename(reports_old = reports) |> select(week, reports_old)
meas_clean <- left_join(meas_y, meas_o, by = "week")

cat("\nClean data loaded: T =", nrow(meas_clean), "weeks\n")

# ── Load baseline ─────────────────────────────────────────────────────────────
baseline_chain <- load_baseline(wave_name)
cat("Baseline posterior loaded.\n")
cat("ESS:", round(effectiveSize(baseline_chain)), "\n")

all_results  <- list()
summary_rows <- list()
summary_rows[[1]] <- posterior_summary(
  baseline_chain, wave_name, "baseline", "0", cfg$params_est
)

# ── 1. Systematic proportional reduction ──────────────────────────────────────
cat("\n", strrep("-", 60), "\n")
cat("TYPE 1: Systematic proportional reduction\n")
cat(strrep("-", 60), "\n")

for (f in f_levels) {
  run_id   <- sprintf("%s_systematic_f%.2f", wave_name, f)
  rds_path <- file.path("results", paste0(run_id, ".rds"))

  if (file.exists(rds_path)) {
    cat("[SKIP]", run_id, "\n"); res <- readRDS(rds_path)
  } else {
    cat(sprintf("[RUN] f = %.2f\n", f))
    meas_deg <- degrade_systematic(meas_clean, f)
    write_csv(meas_deg, file.path("results/degraded_data",
                                   sprintf("%s_systematic_f%.2f.csv", wave_name, f)))
    po  <- build_pomp(cfg, meas_deg)
    po  <- pomp(po, dprior = make_prior(cfg), paramnames = names(coef(po)))
    ll_check <- replicate(5, logLik(pfilter(po, Np = cfg$Np)))
    cat(sprintf("      pfilter LL: mean=%.2f, SD=%.3f\n",
                mean(ll_check, na.rm=TRUE), sd(ll_check, na.rm=TRUE)))
    res <- run_pmcmc(po, cfg, verbose = TRUE)
    saveRDS(res, rds_path)
  }
  all_results[[run_id]] <- res
  summary_rows[[length(summary_rows)+1]] <- posterior_summary(
    res$chain, wave_name, "systematic", f, cfg$params_est)
  cat("      ESS:", paste(round(res$ess), collapse=" "), "\n\n")
}

# ── 2. Stochastic under-reporting ─────────────────────────────────────────────
cat("\n", strrep("-", 60), "\n")
cat("TYPE 2: Stochastic under-reporting\n")
cat(strrep("-", 60), "\n")

for (f in f_levels) {
  run_id   <- sprintf("%s_stochastic_f%.2f", wave_name, f)
  rds_path <- file.path("results", paste0(run_id, ".rds"))

  if (file.exists(rds_path)) {
    cat("[SKIP]", run_id, "\n"); res <- readRDS(rds_path)
  } else {
    cat(sprintf("[RUN] f = %.2f\n", f))
    deg_seed <- as.integer(paste0(2, round(f * 100)))
    meas_deg <- degrade_stochastic(meas_clean, f, seed = deg_seed)
    write_csv(meas_deg, file.path("results/degraded_data",
                                   sprintf("%s_stochastic_f%.2f.csv", wave_name, f)))
    po  <- build_pomp(cfg, meas_deg)
    po  <- pomp(po, dprior = make_prior(cfg), paramnames = names(coef(po)))
    ll_check <- replicate(5, logLik(pfilter(po, Np = cfg$Np)))
    cat(sprintf("      pfilter LL: mean=%.2f, SD=%.3f\n",
                mean(ll_check, na.rm=TRUE), sd(ll_check, na.rm=TRUE)))
    res <- run_pmcmc(po, cfg, verbose = TRUE)
    saveRDS(res, rds_path)
  }
  all_results[[run_id]] <- res
  summary_rows[[length(summary_rows)+1]] <- posterior_summary(
    res$chain, wave_name, "stochastic", f, cfg$params_est)
  cat("      ESS:", paste(round(res$ess), collapse=" "), "\n\n")
}

# ── 3. Temporal aggregation ───────────────────────────────────────────────────
cat("\n", strrep("-", 60), "\n")
cat("TYPE 3: Temporal aggregation\n")
cat(strrep("-", 60), "\n")

for (block in agg_blocks) {
  run_id   <- sprintf("%s_aggregate_b%d", wave_name, block)
  rds_path <- file.path("results", paste0(run_id, ".rds"))

  if (file.exists(rds_path)) {
    cat("[SKIP]", run_id, "\n"); res <- readRDS(rds_path)
  } else {
    cat(sprintf("[RUN] %d-week blocks\n", block))
    agg <- degrade_aggregate(meas_clean, block)
    cat(sprintf("      Observation points: %d\n", nrow(agg$data)))
    write_csv(agg$data, file.path("results/degraded_data",
                                   sprintf("%s_aggregate_b%d.csv", wave_name, block)))
    po  <- build_pomp(cfg, agg$data, obs_times = agg$obs_times)
    po  <- pomp(po, dprior = make_prior(cfg), paramnames = names(coef(po)))
    ll_check <- replicate(5, logLik(pfilter(po, Np = cfg$Np)))
    cat(sprintf("      pfilter LL: mean=%.2f, SD=%.3f\n",
                mean(ll_check, na.rm=TRUE), sd(ll_check, na.rm=TRUE)))
    if (sd(ll_check, na.rm=TRUE) > 1.5)
      warning("pfilter SD high — Np may need increasing")
    res <- run_pmcmc(po, cfg, verbose = TRUE)
    saveRDS(res, rds_path)
  }
  all_results[[run_id]] <- res
  summary_rows[[length(summary_rows)+1]] <- posterior_summary(
    res$chain, wave_name, "aggregate", block, cfg$params_est)
  cat("      ESS:", paste(round(res$ess), collapse=" "), "\n\n")
}

# ── Compile and save ──────────────────────────────────────────────────────────
summary_df <- bind_rows(summary_rows)
saveRDS(summary_df, "results/wave2_summary.rds")
write_csv(summary_df, "results/wave2_summary.csv")

baseline_df <- summary_df |>
  filter(degradation == "baseline") |>
  select(parameter, base_mean=post_mean,
         base_ci_lo=ci_lo, base_ci_hi=ci_hi, base_width=ci_width)

metrics_df <- summary_df |>
  filter(degradation != "baseline") |>
  left_join(baseline_df, by = "parameter") |>
  mutate(
    mean_shift     = post_mean - base_mean,
    mean_shift_pct = 100 * (post_mean - base_mean) / abs(base_mean),
    ci_width_ratio = ci_width / base_width,
    coverage       = (post_mean >= base_ci_lo) & (post_mean <= base_ci_hi),
    degradation    = factor(degradation,
                            levels = c("systematic","stochastic","aggregate"),
                            labels = c("Systematic","Stochastic","Aggregation")),
    level_num      = suppressWarnings(as.numeric(level))
  )

saveRDS(metrics_df, "results/wave2_metrics.rds")
write_csv(metrics_df, "results/wave2_metrics.csv")

cat("\n", strrep("=", 70), "\n")
cat("WAVE 2 COMPLETE\n")
cat("Results saved to results/wave2_summary.rds and results/wave2_metrics.rds\n")
cat("Next: source('sensitivity_plots_wave2.R')\n")
cat(strrep("=", 70), "\n")
