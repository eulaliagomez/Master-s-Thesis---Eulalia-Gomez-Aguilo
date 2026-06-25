# =============================================================================
# run_wave1_sensitivity.R
#
# Standalone runner: sensitivity analysis for Wave 1 only.
# Sources sensitivity_pipeline.R for all shared functions and config,
# then executes only the Wave 1 degradation runs.
#
# USAGE:
#   1. Make sure results/baseline_wave1.rds exists (see Step 1 below)
#   2. setwd() to your project folder
#   3. source("run_wave1_sensitivity.R")
#
# Step 1 — save baseline if not already done:
#   source("covid19_age_wave1.R")          # runs PMCMC, creates prod_chain
#   dir.create("results", showWarnings=FALSE)
#   saveRDS(prod_chain, "results/baseline_wave1.rds")
# =============================================================================

# ── Load all shared functions and config from the pipeline ────────────────────
# sys.source() evaluates the pipeline file in this script's environment
# but we wrap the main loop in a guard so it doesn't execute automatically.
# We do this by temporarily defining a flag before sourcing.

SENSITIVITY_FUNCTIONS_ONLY <- TRUE          # pipeline checks this flag
source("sensitivity_pipeline.R")            # loads functions + wave_config
rm(SENSITIVITY_FUNCTIONS_ONLY)

# ── Directories ───────────────────────────────────────────────────────────────
if (!dir.exists("results"))                 dir.create("results")
if (!dir.exists("results/degraded_data"))   dir.create("results/degraded_data")

# ── Degradation levels ────────────────────────────────────────────────────────
f_levels   <- c(0.10, 0.30, 0.50)
agg_blocks <- c(2, 4)
dom_eigen  <- 4.074

# ── Wave 1 configuration ──────────────────────────────────────────────────────
wave_name <- "wave1"
cfg       <- wave_config[[wave_name]]

cat("\n", strrep("=", 70), "\n")
cat("SENSITIVITY ANALYSIS — WAVE 1 ONLY\n")
cat(cfg$label, "\n")
cat(strrep("=", 70), "\n")

# ── Load clean data ───────────────────────────────────────────────────────────
meas_y <- read_csv(cfg$file_y, show_col_types = FALSE) |>
  rename(reports_young = reports) |> select(week, reports_young)
meas_o <- read_csv(cfg$file_o, show_col_types = FALSE) |>
  rename(reports_old = reports) |> select(week, reports_old)
meas_clean <- left_join(meas_y, meas_o, by = "week")

cat("\nClean data loaded: T =", nrow(meas_clean), "weeks\n")
cat("Young: total =", sum(meas_clean$reports_young),
    "| peak =", max(meas_clean$reports_young),
    "at wk", which.max(meas_clean$reports_young), "\n")
cat("Old:   total =", sum(meas_clean$reports_old),
    "| peak =", max(meas_clean$reports_old),
    "at wk", which.max(meas_clean$reports_old), "\n")

# ── Load baseline posterior ───────────────────────────────────────────────────
baseline_chain <- load_baseline(wave_name)    # errors clearly if missing
cat("\nBaseline posterior loaded.\n")
cat("ESS:", round(effectiveSize(baseline_chain)), "\n")

# Collect results
all_results  <- list()
summary_rows <- list()

# Add baseline to summary
summary_rows[[1]] <- posterior_summary(
  baseline_chain, wave_name, "baseline", "0", cfg$params_est
)

# ── 1. Systematic proportional reduction ──────────────────────────────────────
cat("\n", strrep("-", 60), "\n")
cat("TYPE 1: Systematic proportional reduction\n")
cat(strrep("-", 60), "\n")
cat("Theory: rho shifts by 1/(1-f); beta unchanged\n\n")

for (f in f_levels) {

  run_id   <- sprintf("%s_systematic_f%.2f", wave_name, f)
  rds_path <- file.path("results", paste0(run_id, ".rds"))

  if (file.exists(rds_path)) {
    cat("[SKIP]", run_id, "(already exists)\n")
    res <- readRDS(rds_path)
  } else {
    cat(sprintf("[RUN] f = %.2f (remove %.0f%% of counts)\n", f, f*100))

    meas_deg <- degrade_systematic(meas_clean, f)

    # Sanity check: shape preserved
    ratio_check <- cor(meas_clean$reports_young,
                       meas_deg$reports_young, use = "complete.obs")
    cat(sprintf("      Shape check (Pearson r young): %.4f (should be 1.0)\n",
                ratio_check))

    write_csv(meas_deg,
              file.path("results/degraded_data",
                        sprintf("%s_systematic_f%.2f.csv", wave_name, f)))

    po  <- build_pomp(cfg, meas_deg)
    po  <- pomp(po, dprior = make_prior(cfg), paramnames = names(coef(po)))

    # pfilter diagnostic
    ll_check <- replicate(5, logLik(pfilter(po, Np = cfg$Np)))
    cat(sprintf("      pfilter LL: mean=%.2f, SD=%.3f\n",
                mean(ll_check, na.rm=TRUE), sd(ll_check, na.rm=TRUE)))
    if (sd(ll_check, na.rm=TRUE) > 1.0)
      warning("pfilter SD > 1.0 — consider increasing Np")

    res <- run_pmcmc(po, cfg, verbose = TRUE)
    saveRDS(res, rds_path)

    # Theoretical prediction check
    baseline_means <- colMeans(as.data.frame(baseline_chain))
    deg_means      <- colMeans(as.data.frame(res$chain))
    cat(sprintf("      Predicted rho_y shift: x%.3f | Observed: x%.3f\n",
                1/(1-f),
                deg_means["rho_y"] / baseline_means["rho_y"]))
    cat(sprintf("      Predicted beta1 shift: x1.000 | Observed: x%.3f\n",
                deg_means["Beta1"] / baseline_means["Beta1"]))
  }

  all_results[[run_id]] <- res
  summary_rows[[length(summary_rows)+1]] <- posterior_summary(
    res$chain, wave_name, "systematic", f, cfg$params_est
  )
  cat("      ESS:", paste(round(res$ess), collapse=" "), "\n\n")
}

# ── 2. Stochastic under-reporting ─────────────────────────────────────────────
cat("\n", strrep("-", 60), "\n")
cat("TYPE 2: Stochastic under-reporting\n")
cat(strrep("-", 60), "\n")
cat("Theory: same expected shift as systematic + additional noise\n\n")

for (f in f_levels) {

  run_id   <- sprintf("%s_stochastic_f%.2f", wave_name, f)
  rds_path <- file.path("results", paste0(run_id, ".rds"))

  if (file.exists(rds_path)) {
    cat("[SKIP]", run_id, "(already exists)\n")
    res <- readRDS(rds_path)
  } else {
    cat(sprintf("[RUN] f = %.2f\n", f))

    deg_seed <- as.integer(paste0(1, round(f * 100)))   # wave1 = prefix 1
    meas_deg <- degrade_stochastic(meas_clean, f, seed = deg_seed)

    # Noise added vs systematic
    meas_sys <- degrade_systematic(meas_clean, f)
    noise_sd  <- sd(meas_deg$reports_young - meas_sys$reports_young)
    cat(sprintf("      Extra noise SD (young): %.1f counts/week\n", noise_sd))

    write_csv(meas_deg,
              file.path("results/degraded_data",
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
    res$chain, wave_name, "stochastic", f, cfg$params_est
  )
  cat("      ESS:", paste(round(res$ess), collapse=" "), "\n\n")
}

# ── 3. Temporal aggregation ───────────────────────────────────────────────────
cat("\n", strrep("-", 60), "\n")
cat("TYPE 3: Temporal aggregation\n")
cat(strrep("-", 60), "\n")
cat("Theory: beta CI widens (shape information lost); rho mildly affected\n\n")

for (block in agg_blocks) {

  run_id   <- sprintf("%s_aggregate_b%d", wave_name, block)
  rds_path <- file.path("results", paste0(run_id, ".rds"))

  if (file.exists(rds_path)) {
    cat("[SKIP]", run_id, "(already exists)\n")
    res <- readRDS(rds_path)
  } else {
    cat(sprintf("[RUN] %d-week blocks\n", block))

    agg <- degrade_aggregate(meas_clean, block)
    cat(sprintf("      Observation points: %d (from %d weekly)\n",
                nrow(agg$data), nrow(meas_clean)))
    print(agg$data)

    write_csv(agg$data,
              file.path("results/degraded_data",
                        sprintf("%s_aggregate_b%d.csv", wave_name, block)))

    po  <- build_pomp(cfg, agg$data, obs_times = agg$obs_times)
    po  <- pomp(po, dprior = make_prior(cfg), paramnames = names(coef(po)))

    ll_check <- replicate(5, logLik(pfilter(po, Np = cfg$Np)))
    cat(sprintf("      pfilter LL: mean=%.2f, SD=%.3f\n",
                mean(ll_check, na.rm=TRUE), sd(ll_check, na.rm=TRUE)))
    if (sd(ll_check, na.rm=TRUE) > 1.5)
      warning("pfilter SD high for aggregated data — Np may need increasing")

    res <- run_pmcmc(po, cfg, verbose = TRUE)
    saveRDS(res, rds_path)
  }

  all_results[[run_id]] <- res
  summary_rows[[length(summary_rows)+1]] <- posterior_summary(
    res$chain, wave_name, "aggregate", block, cfg$params_est
  )
  cat("      ESS:", paste(round(res$ess), collapse=" "), "\n\n")
}

# ── Compile Wave 1 summary ────────────────────────────────────────────────────
summary_df <- bind_rows(summary_rows)
saveRDS(summary_df, "results/wave1_summary.rds")
write_csv(summary_df, "results/wave1_summary.csv")

# ── Sensitivity metrics ───────────────────────────────────────────────────────
baseline_df <- summary_df |>
  filter(degradation == "baseline") |>
  select(parameter, base_mean=post_mean,
         base_ci_lo=ci_lo, base_ci_hi=ci_hi, base_width=ci_width)

metrics_df <- summary_df |>
  filter(degradation != "baseline") |>
  left_join(baseline_df, by = "parameter") |>
  mutate(
    mean_shift      = post_mean - base_mean,
    mean_shift_pct  = 100 * (post_mean - base_mean) / abs(base_mean),
    ci_width_ratio  = ci_width / base_width,
    coverage        = (post_mean >= base_ci_lo) & (post_mean <= base_ci_hi),
    degradation     = factor(degradation,
                             levels = c("systematic","stochastic","aggregate"),
                             labels = c("Systematic","Stochastic","Aggregation")),
    level_num       = suppressWarnings(as.numeric(level))
  )

saveRDS(metrics_df, "results/wave1_metrics.rds")
write_csv(metrics_df, "results/wave1_metrics.csv")

# ── Console summary ───────────────────────────────────────────────────────────
cat("\n", strrep("=", 70), "\n")
cat("WAVE 1 COMPLETE\n")
cat(strrep("=", 70), "\n\n")

cat("--- Posterior summary ---\n")
summary_df |>
  mutate(Rt = ifelse(grepl("Beta", parameter),
                     round(post_mean * dom_eigen, 3), NA)) |>
  select(degradation, level, parameter, post_mean, ci_lo, ci_hi, Rt) |>
  print(n = 60)

cat("\n--- Key sensitivity metrics ---\n")
metrics_df |>
  select(degradation, level, parameter,
         mean_shift_pct, ci_width_ratio, coverage) |>
  arrange(degradation, level, parameter) |>
  print(n = 60)

cat("\n--- Theoretical prediction check (systematic reduction) ---\n")
cat("Expected: rho shifts by 1/(1-f); beta unchanged\n\n")
metrics_df |>
  filter(degradation == "Systematic") |>
  mutate(
    f = as.numeric(level),
    theoretical = ifelse(grepl("rho", parameter),
                         round(1/(1-f), 3), 1.000),
    observed    = round(1 + mean_shift_pct/100, 3)
  ) |>
  select(level, parameter, theoretical, observed, ci_width_ratio) |>
  print(n = 20)

cat("\nAll results saved to results/\n")
cat("  wave1_summary.rds / .csv\n")
cat("  wave1_metrics.rds / .csv\n")
cat("  degraded_data/wave1_*.csv\n")
cat("\nNext: source('sensitivity_plots.R') to generate figures\n")
