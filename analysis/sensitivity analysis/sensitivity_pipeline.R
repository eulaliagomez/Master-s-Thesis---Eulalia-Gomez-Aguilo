# =============================================================================
# sensitivity_pipeline.R
#
# Sensitivity analysis pipeline for age-structured SEIR COVID-19 models
# Eulàlia Gómez Aguiló — MSc EMOS Thesis — Leiden University / CBS
#
# PURPOSE:
#   For each of the 4 epidemic waves, apply 3 data degradation types at
#   specified levels, refit the identical PMCMC pipeline, and save posterior
#   summaries for comparison against the clean-data baseline.
#
# DEGRADATION TYPES:
#   1. Systematic proportional reduction  (deterministic scale reduction)
#   2. Stochastic under-reporting         (binomial thinning)
#   3. Temporal aggregation               (2-week and 4-week blocks)
#
# DESIGN:
#   - 4 waves × 3 types × (3 + 3 + 2) levels = 32 degradation runs
#   - Baseline posteriors loaded from .rds files saved by the wave scripts
#     (run those first and save with save_baseline() below)
#   - All PMCMC settings identical to baseline except Nmcmc is adaptive:
#     Phase 2 runs until ESS >= 200 for all params, up to max_iter
#   - Results saved as a structured list and as a flat comparison data frame
#
# WORKFLOW:
#   Step 1: Source baseline wave scripts and save posteriors
#           source("covid19_age_wave1.R"); save_baseline(prod_chain, "wave1")
#           ... repeat for waves 2-4 ...
#   Step 2: source("sensitivity_pipeline.R")
#   Step 3: source("sensitivity_plots.R")   [companion plotting script]
#
# OUTPUT FILES:
#   results/baseline_wave{1-4}.rds        — baseline posterior chains
#   results/degrade_{type}_w{wave}_f{level}.rds  — degraded posterior chains
#   results/sensitivity_summary.rds       — flat comparison data frame
#   results/sensitivity_summary.csv       — same, for inspection
#
# =============================================================================

library(tidyverse)
library(pomp)
library(coda)

set.seed(20240601)

# ── 0. DIRECTORY SETUP ────────────────────────────────────────────────────────
if (!dir.exists("results")) dir.create("results")
if (!dir.exists("results/degraded_data")) dir.create("results/degraded_data")

# ── 1. WAVE CONFIGURATION ─────────────────────────────────────────────────────
# All wave-specific settings in one place.
# Matches the final wave scripts exactly.

C_yy <- 9.4901; C_yo <- 0.7523; C_oy <- 2.8524; C_oo <- 3.2609
N_y  <- 12736000; N_o <- 3359000
dom_eigen <- 4.074

wave_config <- list(

  wave1 = list(
    label      = "Wave 1 (Ancestral)",
    file_y     = "covid19_wave1_young_NL.csv",
    file_o     = "covid19_wave1_old_NL.csv",
    params_est = c("Beta1", "Beta2", "Beta3", "rho_y", "rho_o"),
    theta_start = c(
      Beta1=0.520, Beta2=0.120, Beta3=0.180,
      mu_EI=1.37, mu_IR=2.41, eta_y=0.99, eta_o=0.99,
      rho_y=0.060, rho_o=0.300, k_y=10, k_o=10,
      N_y=N_y, N_o=N_o, C_yy=C_yy, C_yo=C_yo, C_oy=C_oy, C_oo=C_oo
    ),
    npi_df = data.frame(week=0:22,
                        npi_phase=c(rep(0,4), rep(1,6), rep(2,13))),
    prior_mu  = c(Beta1=-0.654, Beta2=-2.120, Beta3=-1.715,
                  rho_y=-2.752, rho_o=-0.847),
    prior_sd  = c(Beta1=0.25, Beta2=0.55, Beta3=0.55,
                  rho_y=0.70, rho_o=0.70),
    d = 5, Np = 3000,
    # Initial conditions at MLE rho
    I0y=2310, E0y=4063, I0o=279, E0o=491,
    n_phases = 3
  ),

  wave2 = list(
    label      = "Wave 2 (Ancestral, double hump)",
    file_y     = "covid19_wave2_young_NL.csv",
    file_o     = "covid19_wave2_old_NL.csv",
    params_est = c("Beta1", "Beta2", "Beta3", "Beta4", "rho_y", "rho_o"),
    theta_start = c(
      Beta1=0.340, Beta2=0.200, Beta3=0.430, Beta4=0.220,
      mu_EI=1.37, mu_IR=2.41, eta_y=0.90, eta_o=0.86,
      rho_y=0.160, rho_o=0.360, k_y=30, k_o=10,
      N_y=N_y, N_o=N_o, C_yy=C_yy, C_yo=C_yo, C_oy=C_oy, C_oo=C_oo
    ),
    npi_df = data.frame(week=0:21,
                        npi_phase=c(rep(0,6), rep(1,4), rep(2,4), rep(3,8))),
    prior_mu  = c(Beta1=-1.079, Beta2=-1.609, Beta3=-0.844, Beta4=-1.514,
                  rho_y=-1.658, rho_o=-0.575),
    prior_sd  = c(Beta1=0.20, Beta2=0.25, Beta3=0.20, Beta4=0.25,
                  rho_y=0.20, rho_o=0.35),
    d = 6, Np = 3000,
    I0y=28620, E0y=50347, I0o=1924, E0o=3384,
    n_phases = 4
  ),

  wave3 = list(
    label      = "Wave 3 (Alpha)",
    file_y     = "covid19_wave3_young_NL.csv",
    file_o     = "covid19_wave3_old_NL.csv",
    params_est = c("Beta1", "Beta2", "Beta3", "rho_y", "rho_o"),
    theta_start = c(
      Beta1=0.380, Beta2=0.360, Beta3=0.260,
      mu_EI=1.37, mu_IR=2.41, eta_y=0.72, eta_o=0.55,
      rho_y=0.240, rho_o=0.400, k_y=30, k_o=10,
      N_y=N_y, N_o=N_o, C_yy=C_yy, C_yo=C_yo, C_oy=C_oy, C_oo=C_oo
    ),
    npi_df = data.frame(week=0:22,
                        npi_phase=c(rep(0,4), rep(1,6), rep(2,13))),
    prior_mu  = c(Beta1=-0.968, Beta2=-1.022, Beta3=-1.347,
                  rho_y=-1.153, rho_o=-0.405),
    prior_sd  = c(Beta1=0.30, Beta2=0.35, Beta3=0.30,
                  rho_y=0.20, rho_o=0.35),
    d = 5, Np = 3000,
    I0y=32312, E0y=56840, I0o=6223, E0o=10947,
    n_phases = 3
  ),

  wave4 = list(
    label      = "Wave 4 (Delta)",
    file_y     = "covid19_wave4_young_NL.csv",
    file_o     = "covid19_wave4_old_NL.csv",
    params_est = c("Beta1", "Beta2", "Beta3", "rho_y", "rho_o"),
    theta_start = c(
      Beta1=1.426, Beta2=0.240, Beta3=0.450,
      mu_EI=1.37, mu_IR=2.41, eta_y=0.64, eta_o=0.25,
      rho_y=0.140, rho_o=0.320, k_y=30, k_o=10,
      N_y=N_y, N_o=N_o, C_yy=C_yy, C_yo=C_yo, C_oy=C_oy, C_oo=C_oo
    ),
    npi_df = data.frame(week=0:18,
                        npi_phase=c(rep(0,2), rep(1,5), rep(2,12))),
    prior_mu  = c(Beta1=0.355, Beta2=-1.427, Beta3=-0.799,
                  rho_y=-1.815, rho_o=-0.754),
    prior_sd  = c(Beta1=0.35, Beta2=0.40, Beta3=0.40,
                  rho_y=0.20, rho_o=0.35),
    d = 5, Np = 3000,
    I0y=12297, E0y=21632, I0o=296, E0o=520,
    n_phases = 3
  )
)

# ── 2. SHARED pomp COMPONENTS (Csnippets) ─────────────────────────────────────
# These are identical across all waves and all degradation runs.

seir_age_step <- Csnippet("
  double eff_beta;
  if      (npi_phase < 0.5) eff_beta = Beta1;
  else if (npi_phase < 1.5) eff_beta = Beta2;
  else if (npi_phase < 2.5) eff_beta = Beta3;
  else                       eff_beta = Beta4;

  double lam_y = eff_beta * dt * (C_yy*I_y/N_y + C_yo*I_o/N_o);
  double lam_o = eff_beta * dt * (C_oy*I_y/N_y + C_oo*I_o/N_o);

  double dN_SE_y = rbinom(S_y, 1.0 - exp(-lam_y));
  double dN_EI_y = rbinom(E_y, 1.0 - exp(-mu_EI*dt));
  double dN_IR_y = rbinom(I_y, 1.0 - exp(-mu_IR*dt));
  double dN_SE_o = rbinom(S_o, 1.0 - exp(-lam_o));
  double dN_EI_o = rbinom(E_o, 1.0 - exp(-mu_EI*dt));
  double dN_IR_o = rbinom(I_o, 1.0 - exp(-mu_IR*dt));

  S_y -= dN_SE_y; E_y += dN_SE_y-dN_EI_y; I_y += dN_EI_y-dN_IR_y;
  R_y += dN_IR_y; H_y += dN_IR_y;
  S_o -= dN_SE_o; E_o += dN_SE_o-dN_EI_o; I_o += dN_EI_o-dN_IR_o;
  R_o += dN_IR_o; H_o += dN_IR_o;
")

# Age-specific NegBin observation model
covid_dmeas <- Csnippet("
  lik = dnbinom_mu(reports_young, k_y, rho_y * H_y, 1)
      + dnbinom_mu(reports_old,   k_o, rho_o * H_o, 1);
  if (!give_log) lik = exp(lik);
")
covid_rmeas <- Csnippet("
  reports_young = rnbinom_mu(k_y, rho_y * H_y);
  reports_old   = rnbinom_mu(k_o, rho_o * H_o);
")

seir_rinit <- Csnippet("
  S_y = nearbyint(eta_y*N_y) - nearbyint(I0y) - nearbyint(E0y);
  E_y = nearbyint(E0y); I_y = nearbyint(I0y);
  R_y = nearbyint((1.0-eta_y)*N_y); H_y = 0.0;
  S_o = nearbyint(eta_o*N_o) - nearbyint(I0o) - nearbyint(E0o);
  E_o = nearbyint(E0o); I_o = nearbyint(I0o);
  R_o = nearbyint((1.0-eta_o)*N_o); H_o = 0.0;
")

# ── 3. HELPER: BUILD pomp OBJECT ──────────────────────────────────────────────
# Builds a pomp object from a config and a (possibly degraded) dataset.
# The times argument is overrideable for temporal aggregation.

build_pomp <- function(cfg, meas_df, obs_times = NULL) {

  # Default: weekly observation times from the data
  if (is.null(obs_times)) obs_times <- meas_df$week

  # Beta4 always included in paramnames and theta so the shared Csnippet
  # compiles on all waves. For 3-phase waves it is fixed at 0 and the
  # npi_phase covariate never reaches value 3, so the else branch is dead code.
  all_pnames <- c(
    "Beta1","Beta2","Beta3","Beta4",
    "mu_EI","mu_IR","eta_y","eta_o","rho_y","rho_o",
    "k_y","k_o","N_y","N_o","C_yy","C_yo","C_oy","C_oo",
    "I0y","E0y","I0o","E0o"
  )

  theta <- c(cfg$theta_start,
             I0y = cfg$I0y, E0y = cfg$E0y,
             I0o = cfg$I0o, E0o = cfg$E0o)
  # Ensure Beta4 is always present (= 0 for 3-phase waves)
  if (!"Beta4" %in% names(theta)) theta["Beta4"] <- 0

  npi_covar <- covariate_table(cfg$npi_df, times = "week")

  po <- meas_df |>
    filter(week %in% obs_times) |>
    pomp(
      times      = "week",
      t0         = 0,   # always start epidemic from week 0
      rprocess   = euler(seir_age_step, delta.t = 1/7),
      rinit      = seir_rinit,
      rmeasure   = covid_rmeas,
      dmeasure   = covid_dmeas,
      covar      = npi_covar,
      partrans   = parameter_trans(
        log   = intersect(c("Beta1","Beta2","Beta3","Beta4"),
                          cfg$params_est),
        logit = c("rho_y","rho_o")
      ),
      paramnames = all_pnames,
      statenames = c("S_y","E_y","I_y","R_y","H_y",
                     "S_o","E_o","I_o","R_o","H_o"),
      accumvars  = c("H_y","H_o"),
      obsnames   = c("reports_young","reports_old")
    )
  coef(po) <- theta
  po
}

# ── 4. HELPER: PRIOR Csnippet ─────────────────────────────────────────────────
make_prior <- function(cfg) {
  mu  <- cfg$prior_mu
  sd  <- cfg$prior_sd
  has_B4 <- "Beta4" %in% names(mu)

  lines <- c(
    "double lB1    = log(Beta1);",
    "double lB2    = log(Beta2);",
    "double lB3    = log(Beta3);",
    if (has_B4) "double lB4    = log(Beta4 + 1e-10);",
    "double lrho_y = log(rho_y / (1.0 - rho_y));",
    "double lrho_o = log(rho_o / (1.0 - rho_o));",
    "",
    sprintf("lik = dnorm(lB1,    %.4f, %.2f, 1)", mu["Beta1"], sd["Beta1"]),
    sprintf("    + dnorm(lB2,   %.4f, %.2f, 1)", mu["Beta2"], sd["Beta2"]),
    sprintf("    + dnorm(lB3,   %.4f, %.2f, 1)", mu["Beta3"], sd["Beta3"]),
    if (has_B4) sprintf("    + dnorm(lB4,   %.4f, %.2f, 1)", mu["Beta4"], sd["Beta4"]),
    sprintf("    + dnorm(lrho_y,%.4f, %.2f, 1)", mu["rho_y"], sd["rho_y"]),
    sprintf("    + dnorm(lrho_o,%.4f, %.2f, 1);", mu["rho_o"], sd["rho_o"]),
    "",
    "if (!give_log) lik = exp(lik);"
  )
  Csnippet(paste(lines, collapse = "\n"))
}

# ── 5. HELPER: RUN PMCMC ──────────────────────────────────────────────────────
# Phase 1: diagonal RW, 5000 iter (exploration)
# Phase 2: MVN from Phase-1 covariance, adaptive length
#   - starts at Nmcmc_min
#   - checks ESS every check_every iterations
#   - stops when ESS >= 200 for all params OR Nmcmc_max reached

run_pmcmc <- function(po, cfg,
                      Nmcmc_p1   = 5000,
                      Nmcmc_min  = 5000,
                      Nmcmc_max  = 15000,
                      check_every = 2000,
                      ess_target = 200,
                      verbose    = TRUE) {

  params_est <- cfg$params_est
  d          <- cfg$d
  Np         <- cfg$Np

  # ── Phase 1 ────────────────────────────────────────────────────────────────
  if (verbose) cat("  Phase 1 (", Nmcmc_p1, " iter)...\n", sep="")

  step_sizes <- setNames(rep(0.10, d), params_est)
  p1 <- pmcmc(po,
              Nmcmc    = Nmcmc_p1,
              Np       = Np,
              proposal = mvn_diag_rw(step_sizes^2))
  acc1 <- p1@accepts / Nmcmc_p1
  if (verbose) cat("  Phase 1 acceptance:", round(acc1, 3), "\n")

  # Estimate posterior covariance from second half of Phase 1
  chain1    <- as.matrix(as.data.frame(traces(p1)))
  warmup    <- floor(Nmcmc_p1 / 2)
  post_cov  <- cov(chain1[warmup:Nmcmc_p1, params_est, drop = FALSE])
  opt_scale <- (2.38^2) / d

  # ── Phase 2: adaptive length ───────────────────────────────────────────────
  if (verbose) cat("  Phase 2 (adaptive, ESS target =", ess_target, ")...\n")

  current_chain <- p1
  total_iter    <- 0

  # Phase 2: run in blocks, checking ESS after each block
  # We restart from scratch each block rather than chaining, to avoid
  # the trace concatenation indexing problem.
  # Collect all Phase 2 draws in a list and rbind at the end.
  p2_draws <- list()
  total_iter <- 0

  repeat {
    n_new <- min(check_every, Nmcmc_max - total_iter)
    if (n_new <= 0) break

    block_chain <- pmcmc(current_chain,
                         Nmcmc    = n_new,
                         Np       = Np,
                         proposal = mvn_rw(opt_scale * post_cov))
    total_iter    <- total_iter + n_new
    current_chain <- block_chain   # use last state as starting point

    # Extract only the new draws from this block
    block_mat  <- as.matrix(as.data.frame(traces(block_chain)))
    # traces() includes the starting value as row 1, so skip it
    new_rows   <- seq(2, nrow(block_mat))
    p2_draws[[length(p2_draws) + 1]] <- block_mat[new_rows,
                                                    params_est,
                                                    drop = FALSE]

    # Check ESS on all Phase 2 draws accumulated so far
    all_p2 <- do.call(rbind, p2_draws)
    if (nrow(all_p2) < 100) next

    prod_sub <- as.mcmc(all_p2)
    ess_vals <- effectiveSize(prod_sub)

    if (verbose) {
      cat(sprintf("    iter %d: ESS min=%.0f (target %d)\n",
                  total_iter, min(ess_vals), ess_target))
    }

    if (min(ess_vals) >= ess_target) break
    if (total_iter >= Nmcmc_max)     break
  }

  acc2 <- current_chain@accepts / n_new
  if (verbose) {
    cat("  Phase 2 acceptance:", round(acc2, 3),
        "| total iter:", total_iter, "\n")
  }

  # ── Extract production chain ───────────────────────────────────────────────
  all_p2     <- do.call(rbind, p2_draws)
  prod_chain <- as.mcmc(all_p2)
  ess_final  <- effectiveSize(prod_chain)

  if (any(ess_final < ess_target)) {
    warning(sprintf("ESS below target for: %s",
                    paste(names(ess_final)[ess_final < ess_target],
                          collapse = ", ")))
  }

  list(chain      = prod_chain,
       ess        = ess_final,
       acc_p1     = acc1,
       acc_p2     = acc2,
       total_iter = total_iter)
}

# ── 6. HELPER: POSTERIOR SUMMARY ─────────────────────────────────────────────
# Returns a tidy data frame with posterior mean, SD, and 95% CI for each param.

posterior_summary <- function(chain, wave, degradation, level, params_est) {
  df <- as.data.frame(chain)
  map_dfr(params_est, function(p) {
    x <- df[[p]]
    tibble(
      wave        = wave,
      degradation = degradation,
      level       = as.character(level),
      parameter   = p,
      post_mean   = mean(x),
      post_sd     = sd(x),
      ci_lo       = quantile(x, 0.025),
      ci_hi       = quantile(x, 0.975),
      ci_width    = quantile(x, 0.975) - quantile(x, 0.025)
    )
  })
}

# ── 7. DEGRADATION FUNCTIONS ──────────────────────────────────────────────────

# 7a. Systematic proportional reduction
# y* = floor((1 - f) * y)  — deterministic, shape exactly preserved
degrade_systematic <- function(meas_df, f) {
  meas_df |>
    mutate(
      reports_young = floor((1 - f) * reports_young),
      reports_old   = floor((1 - f) * reports_old)
    )
}

# 7b. Stochastic under-reporting
# y* ~ Binomial(y, 1 - f)  — random, shape distorted by binomial noise
degrade_stochastic <- function(meas_df, f, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  meas_df |>
    mutate(
      reports_young = rbinom(n(), reports_young, 1 - f),
      reports_old   = rbinom(n(), reports_old,   1 - f)
    )
}

# 7c. Temporal aggregation
# Aggregate weekly counts into blocks of `block` weeks.
# The observation times passed to build_pomp() become the LAST week
# of each block; the pomp model runs at full weekly resolution but
# the likelihood is only evaluated at block endpoints.
# Returns a list: aggregated data frame + observation times vector.
degrade_aggregate <- function(meas_df, block) {
  T   <- nrow(meas_df)
  # Assign each week to a block; only keep complete blocks
  meas_df <- meas_df |>
    mutate(block_id = ceiling(week / block)) |>
    group_by(block_id) |>
    filter(n() == block) |>           # drop incomplete trailing block
    summarise(
      week          = max(week),      # label = last week of block
      reports_young = sum(reports_young),
      reports_old   = sum(reports_old),
      .groups       = "drop"
    ) |>
    select(week, reports_young, reports_old)

  list(
    data      = meas_df,
    obs_times = meas_df$week
  )
}

# ── 8. HELPER: LOAD OR RUN BASELINE ──────────────────────────────────────────
# Loads baseline posterior from .rds if it exists; otherwise errors with
# instructions to run the wave script first.

load_baseline <- function(wave_name) {
  path <- file.path("results", paste0("baseline_", wave_name, ".rds"))
  if (!file.exists(path)) {
    stop(
      sprintf(
        "Baseline not found: %s\n",
        path
      ),
      sprintf(
        "Run the wave script first, then call:\n",
        "  saveRDS(prod_chain, '%s')\n", path
      )
    )
  }
  readRDS(path)
}

# Convenience: save baseline from within wave script
save_baseline <- function(prod_chain, wave_name) {
  path <- file.path("results", paste0("baseline_", wave_name, ".rds"))
  saveRDS(prod_chain, path)
  cat("Baseline saved:", path, "\n")
}

# ── 9. MAIN PIPELINE ──────────────────────────────────────────────────────────
# Runs all 32 degradation scenarios across all 4 waves.
# Results are collected into `all_results` (list) and `summary_df` (tibble).

# Degradation grid
f_levels   <- c(0.10, 0.30, 0.50)
agg_blocks <- c(2, 4)

# Guard: source with SENSITIVITY_FUNCTIONS_ONLY=TRUE to skip the main loop
if (!exists("SENSITIVITY_FUNCTIONS_ONLY") || !isTRUE(SENSITIVITY_FUNCTIONS_ONLY)) {


all_results <- list()
summary_rows <- list()

for (wave_name in names(wave_config)) {

  cfg <- wave_config[[wave_name]]
  cat("\n", strrep("=", 70), "\n")
  cat("WAVE:", cfg$label, "\n")
  cat(strrep("=", 70), "\n")

  # Load clean data
  meas_y <- read_csv(cfg$file_y, show_col_types = FALSE) |>
    rename(reports_young = reports) |> select(week, reports_young)
  meas_o <- read_csv(cfg$file_o, show_col_types = FALSE) |>
    rename(reports_old   = reports) |> select(week, reports_old)
  meas_clean <- left_join(meas_y, meas_o, by = "week")

  # Load baseline posterior
  baseline_chain <- tryCatch(
    load_baseline(wave_name),
    error = function(e) { warning(e$message); NULL }
  )

  if (!is.null(baseline_chain)) {
    row <- posterior_summary(baseline_chain, wave_name,
                             "baseline", "0", cfg$params_est)
    summary_rows[[length(summary_rows) + 1]] <- row
  }

  # ── 9a. Systematic proportional reduction ─────────────────────────────────
  for (f in f_levels) {

    run_id <- sprintf("%s_systematic_f%.2f", wave_name, f)
    rds_path <- file.path("results", paste0(run_id, ".rds"))

    if (file.exists(rds_path)) {
      cat("\n[SKIP — already exists]", run_id, "\n")
      res <- readRDS(rds_path)
    } else {
      cat(sprintf("\n[RUN] Systematic reduction f=%.2f | %s\n", f, wave_name))

      meas_deg <- degrade_systematic(meas_clean, f)
      # Save degraded data for reproducibility
      write_csv(meas_deg,
                file.path("results/degraded_data",
                          sprintf("%s_systematic_f%.2f.csv", wave_name, f)))

      po  <- build_pomp(cfg, meas_deg)
      po  <- pomp(po, dprior = make_prior(cfg),
                  paramnames = names(coef(po)))

      res <- run_pmcmc(po, cfg, verbose = TRUE)
      saveRDS(res, rds_path)
    }

    all_results[[run_id]] <- res
    row <- posterior_summary(res$chain, wave_name,
                             "systematic", f, cfg$params_est)
    summary_rows[[length(summary_rows) + 1]] <- row
    cat("  ESS:", paste(round(res$ess), collapse = " "), "\n")
  }

  # ── 9b. Stochastic under-reporting ────────────────────────────────────────
  for (f in f_levels) {

    run_id   <- sprintf("%s_stochastic_f%.2f", wave_name, f)
    rds_path <- file.path("results", paste0(run_id, ".rds"))

    if (file.exists(rds_path)) {
      cat("\n[SKIP — already exists]", run_id, "\n")
      res <- readRDS(rds_path)
    } else {
      cat(sprintf("\n[RUN] Stochastic under-reporting f=%.2f | %s\n",
                  f, wave_name))

      # Fixed seed per (wave × f) so degraded data is reproducible
      deg_seed <- as.integer(paste0(
        which(names(wave_config) == wave_name),
        round(f * 100)
      ))
      meas_deg <- degrade_stochastic(meas_clean, f, seed = deg_seed)
      write_csv(meas_deg,
                file.path("results/degraded_data",
                          sprintf("%s_stochastic_f%.2f.csv", wave_name, f)))

      po  <- build_pomp(cfg, meas_deg)
      po  <- pomp(po, dprior = make_prior(cfg),
                  paramnames = names(coef(po)))

      res <- run_pmcmc(po, cfg, verbose = TRUE)
      saveRDS(res, rds_path)
    }

    all_results[[run_id]] <- res
    row <- posterior_summary(res$chain, wave_name,
                             "stochastic", f, cfg$params_est)
    summary_rows[[length(summary_rows) + 1]] <- row
    cat("  ESS:", paste(round(res$ess), collapse = " "), "\n")
  }

  # ── 9c. Temporal aggregation ──────────────────────────────────────────────
  for (block in agg_blocks) {

    run_id   <- sprintf("%s_aggregate_b%d", wave_name, block)
    rds_path <- file.path("results", paste0(run_id, ".rds"))

    if (file.exists(rds_path)) {
      cat("\n[SKIP — already exists]", run_id, "\n")
      res <- readRDS(rds_path)
    } else {
      cat(sprintf("\n[RUN] Temporal aggregation block=%d | %s\n",
                  block, wave_name))

      agg      <- degrade_aggregate(meas_clean, block)
      write_csv(agg$data,
                file.path("results/degraded_data",
                          sprintf("%s_aggregate_b%d.csv", wave_name, block)))

      # build_pomp receives full weekly data for process model
      # but obs_times restricts likelihood evaluation to block endpoints
      # The aggregated COUNTS (not weekly) are what the measurement model sees
      po <- build_pomp(cfg, agg$data, obs_times = agg$obs_times)
      po <- pomp(po, dprior = make_prior(cfg),
                 paramnames = names(coef(po)))

      res <- run_pmcmc(po, cfg, verbose = TRUE)
      saveRDS(res, rds_path)
    }

    all_results[[run_id]] <- res
    row <- posterior_summary(res$chain, wave_name,
                             "aggregate", block, cfg$params_est)
    summary_rows[[length(summary_rows) + 1]] <- row
    cat("  ESS:", paste(round(res$ess), collapse = " "), "\n")
  }
}

# ── 10. COMPILE SUMMARY DATA FRAME ───────────────────────────────────────────
summary_df <- bind_rows(summary_rows)

# Add Rt / R0 derived columns for interpretability
summary_df <- summary_df |>
  mutate(
    Rt = case_when(
      parameter == "Beta1" ~ post_mean * dom_eigen,
      parameter == "Beta2" ~ post_mean * dom_eigen,
      parameter == "Beta3" ~ post_mean * dom_eigen,
      parameter == "Beta4" ~ post_mean * dom_eigen,
      TRUE ~ NA_real_
    )
  )

saveRDS(summary_df, "results/sensitivity_summary.rds")
write_csv(summary_df,  "results/sensitivity_summary.csv")

cat("\n", strrep("=", 70), "\n")
cat("PIPELINE COMPLETE\n")
cat("Results saved to results/sensitivity_summary.rds\n")
cat(sprintf("Total runs completed: %d\n", length(all_results)))
cat(strrep("=", 70), "\n")

# ── 11. COMPUTE SENSITIVITY METRICS ──────────────────────────────────────────
# For each (wave, degradation, level, parameter), compute:
#   (a) posterior mean shift vs baseline
#   (b) 95% CI width ratio vs baseline
#   (c) coverage: does baseline CI contain degraded posterior mean?

baseline_df <- summary_df |>
  filter(degradation == "baseline") |>
  select(wave, parameter,
         base_mean  = post_mean,
         base_ci_lo = ci_lo,
         base_ci_hi = ci_hi,
         base_width = ci_width)

metrics_df <- summary_df |>
  filter(degradation != "baseline") |>
  left_join(baseline_df, by = c("wave", "parameter")) |>
  mutate(
    # (a) absolute shift on natural scale
    mean_shift    = post_mean - base_mean,
    # relative shift (% of baseline mean)
    mean_shift_pct = 100 * (post_mean - base_mean) / abs(base_mean),
    # (b) CI width ratio (>1 = wider = more uncertain)
    ci_width_ratio = ci_width / base_width,
    # (c) coverage: is degraded posterior mean inside baseline CI?
    coverage       = (post_mean >= base_ci_lo) & (post_mean <= base_ci_hi),
    # clean factor for plotting
    degradation    = factor(degradation,
                            levels = c("systematic","stochastic","aggregate"),
                            labels = c("Systematic reduction",
                                       "Stochastic under-reporting",
                                       "Temporal aggregation")),
    level_num      = as.numeric(level)
  )

saveRDS(metrics_df, "results/sensitivity_metrics.rds")
write_csv(metrics_df, "results/sensitivity_metrics.csv")

cat("\nSensitivity metrics saved to results/sensitivity_metrics.rds\n")

# ── 12. QUICK CONSOLE SUMMARY ─────────────────────────────────────────────────
cat("\n--- Mean shift summary (|shift| > 10% of baseline) ---\n")
metrics_df |>
  filter(abs(mean_shift_pct) > 10) |>
  arrange(desc(abs(mean_shift_pct))) |>
  select(wave, degradation, level, parameter,
         mean_shift_pct, ci_width_ratio, coverage) |>
  print(n = 30)

cat("\n--- CI width ratio > 1.5 (substantially wider posterior) ---\n")
metrics_df |>
  filter(ci_width_ratio > 1.5) |>
  arrange(desc(ci_width_ratio)) |>
  select(wave, degradation, level, parameter,
         ci_width_ratio, mean_shift_pct) |>
  print(n = 30)

cat("\n--- Coverage failures (baseline CI does not contain degraded mean) ---\n")
metrics_df |>
  filter(!coverage) |>
  select(wave, degradation, level, parameter,
         base_ci_lo, base_ci_hi, post_mean, mean_shift_pct) |>
  print(n = 30)

} # end SENSITIVITY_FUNCTIONS_ONLY guard
