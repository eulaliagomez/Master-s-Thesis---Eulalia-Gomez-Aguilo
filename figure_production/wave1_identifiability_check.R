# =============================================================================
# wave1_identifiability_check.R
#
# Tests whether the Wave 1 posterior is genuinely driven by the data likelihood
# or is primarily determined by tight priors and good starting values.
#
# DESIGN:
#   Four runs, all for Wave 1:
#     (0) Baseline   — MLE starting values, original tight priors (already done)
#     (A) Off-Low    — βs halved, ρs tripled,    flat priors (SD=2 log/logit)
#     (B) Off-High   — βs doubled, ρs at 0.3x,  flat priors (SD=2 log/logit)
#     (C) Off-Mixed  — β1×3, β2×0.3, β3×0.3, ρy×5, ρo×0.5, flat priors
#
# INTERPRETATION:
#   If off chains converge to baseline posterior → likelihood-driven ✓
#   If off chains agree with each other but not baseline → prior anchoring ✗
#   If off chains don't converge → identifiability problem ✗
#   Parameter-by-parameter differences are also informative
#
# OUTPUT:
#   results/identifiability/wave1_off[A/B/C].rds
#   results/identifiability/wave1_identifiability_comparison.pdf/.png
#   results/identifiability/wave1_identifiability_traces.pdf
#
# USAGE:
#   setwd("C:/Users/Usuario/Desktop/CBS - Internship/Analysis/DATA")
#   source("wave1_identifiability_check.R")
# =============================================================================

library(tidyverse)
library(pomp)
library(coda)
library(patchwork)

# ── 0. Setup ──────────────────────────────────────────────────────────────────
SENSITIVITY_FUNCTIONS_ONLY <- TRUE
source("sensitivity_pipeline.R")
rm(SENSITIVITY_FUNCTIONS_ONLY)

dir.create("results/identifiability", showWarnings = FALSE, recursive = TRUE)

set.seed(2025)

cfg        <- wave_config[["wave1"]]
params_est <- cfg$params_est   # Beta1 Beta2 Beta3 rho_y rho_o
dom_eigen  <- 4.074

# Load clean data
meas_y <- read_csv(cfg$file_y, show_col_types = FALSE) |>
  rename(reports_young = reports) |> select(week, reports_young)
meas_o <- read_csv(cfg$file_o, show_col_types = FALSE) |>
  rename(reports_old = reports) |> select(week, reports_old)
meas   <- left_join(meas_y, meas_o, by = "week")

cat("Data loaded: T =", nrow(meas), "weeks\n")

# ── 1. Baseline (already run — just load) ─────────────────────────────────────
cat("\nLoading baseline...\n")
baseline_chain <- readRDS("results/baseline_wave1.rds")
baseline_df    <- as.data.frame(baseline_chain)
cat("Baseline ESS:", round(effectiveSize(baseline_chain)), "\n\n")

# ── 2. Define off-starting runs ───────────────────────────────────────────────
# MLE starting values (for reference)
mle <- c(Beta1 = 0.520, Beta2 = 0.120, Beta3 = 0.180,
         rho_y = 0.060, rho_o = 0.300)

cat("MLE starting values:\n")
for (p in names(mle)) cat(sprintf("  %s = %.4f\n", p, mle[p]))

off_starts <- list(

  # Off-A: all betas halved (too low transmission), rhos tripled (too high detection)
  offA = c(
    Beta1 = mle["Beta1"] * 0.50,    # 0.260
    Beta2 = mle["Beta2"] * 0.50,    # 0.060
    Beta3 = mle["Beta3"] * 0.50,    # 0.090
    rho_y = min(mle["rho_y"] * 3.0, 0.40),  # 0.180
    rho_o = min(mle["rho_o"] * 2.5, 0.80)   # 0.750
  ),

  # Off-B: all betas doubled (too high transmission), rhos reduced
  offB = c(
    Beta1 = mle["Beta1"] * 2.0,     # 1.040
    Beta2 = mle["Beta2"] * 2.0,     # 0.240
    Beta3 = mle["Beta3"] * 2.0,     # 0.360
    rho_y = mle["rho_y"] * 0.30,    # 0.018
    rho_o = mle["rho_o"] * 0.35     # 0.105
  ),

  # Off-C: pathological mix — wrong direction for each parameter
  offC = c(
    Beta1 = mle["Beta1"] * 3.0,     # 1.560  (way too high)
    Beta2 = mle["Beta2"] * 0.25,    # 0.030  (way too low)
    Beta3 = mle["Beta3"] * 0.25,    # 0.045  (way too low)
    rho_y = min(mle["rho_y"] * 5.0, 0.45),  # 0.300 (much too high)
    rho_o = mle["rho_o"] * 0.20     # 0.060  (much too low)
  )
)

cat("\nOff starting values:\n")
for (nm in names(off_starts)) {
  cat(sprintf("  %s:", nm))
  for (p in names(off_starts[[nm]])) {
    cat(sprintf("  %s=%.3f", p, off_starts[[nm]][p]))
  }
  cat("\n")
}

# ── 3. Flat (non-informative) prior ──────────────────────────────────────────
# Wide Normal on log(beta) and logit(rho): SD = 2.0
# This covers:
#   beta: [exp(-4), exp(4)] = [0.018, 54.6] with 95% probability
#   rho:  [expit(-4), expit(4)] = [0.018, 0.982] with 95% probability
# Centre at generic values NOT at the MLE (avoids replicating the tight prior):
#   log(beta): centred at 0 → beta centre = 1.0
#   logit(rho): centred at -1 → rho centre ≈ 0.27

flat_prior_mu <- c(
  Beta1 = 0.0,    # log scale: exp(0) = 1.0
  Beta2 = 0.0,
  Beta3 = 0.0,
  rho_y = -1.0,   # logit scale: expit(-1) ≈ 0.27
  rho_o = -1.0
)
flat_prior_sd <- c(
  Beta1 = 2.0,
  Beta2 = 2.0,
  Beta3 = 2.0,
  rho_y = 2.0,
  rho_o = 2.0
)

cat("\nFlat prior (non-informative):\n")
cat("  log(beta_j) ~ N(0, 2^2)     → beta 95% range: [0.018, 54.6]\n")
cat("  logit(rho_a) ~ N(-1, 2^2)   → rho  95% range: [0.018, 0.982]\n\n")

# Build flat prior Csnippet
make_flat_prior <- function() {
  prior_code <- paste(
    "double lp = 0;",
    sprintf("lp += dnorm(log(Beta1), %.3f, %.3f, 1);", flat_prior_mu["Beta1"], flat_prior_sd["Beta1"]),
    sprintf("lp += dnorm(log(Beta2), %.3f, %.3f, 1);", flat_prior_mu["Beta2"], flat_prior_sd["Beta2"]),
    sprintf("lp += dnorm(log(Beta3), %.3f, %.3f, 1);", flat_prior_mu["Beta3"], flat_prior_sd["Beta3"]),
    sprintf("lp += dnorm(log(rho_y/(1-rho_y)), %.3f, %.3f, 1);", flat_prior_mu["rho_y"], flat_prior_sd["rho_y"]),
    sprintf("lp += dnorm(log(rho_o/(1-rho_o)), %.3f, %.3f, 1);", flat_prior_mu["rho_o"], flat_prior_sd["rho_o"]),
    "lik = (give_log) ? lp : exp(lp);",
    sep = "\n"
  )
  Csnippet(prior_code)
}

# ── 4. Run PMCMC for each off-starting point ──────────────────────────────────
# Note: with flat priors and off starting values the chain may need more
# Phase 1 iterations to reach the likelihood peak before Phase 2 begins.
# We therefore increase Phase 1 to 8,000 iterations.

run_off <- function(start_name, start_vals) {

  rds_path <- file.path("results/identifiability",
                        sprintf("wave1_%s.rds", start_name))

  if (file.exists(rds_path)) {
    cat(sprintf("[SKIP] %s (already exists)\n", start_name))
    return(readRDS(rds_path))
  }

  cat(strrep("-", 60), "\n")
  cat(sprintf("[RUN] %s\n", start_name))
  cat("Starting values:\n")
  for (p in names(start_vals)) cat(sprintf("  %s = %.4f\n", p, start_vals[p]))
  cat("\n")

  # Build pomp with flat prior
  po <- build_pomp(cfg, meas)
  po <- pomp(po,
             dprior    = make_flat_prior(),
             paramnames = names(coef(po)))

  # Verify pfilter at starting values
  coef_start         <- coef(po)
  coef_start[params_est] <- start_vals
  po_start           <- pomp(po, params = coef_start)

  ll_check <- replicate(5, logLik(pfilter(po_start, Np = cfg$Np)))
  cat(sprintf("  pfilter LL at start: mean=%.2f, SD=%.3f\n",
              mean(ll_check, na.rm = TRUE), sd(ll_check, na.rm = TRUE)))

  if (all(is.nan(ll_check))) {
    cat("  ERROR: pfilter returns NaN at starting values.\n")
    cat("  Trying with Np=5000 for stability...\n")
    ll_check <- replicate(3, logLik(pfilter(po_start, Np = 5000)))
    cat(sprintf("  pfilter LL (Np=5000): mean=%.2f, SD=%.3f\n",
                mean(ll_check, na.rm = TRUE), sd(ll_check, na.rm = TRUE)))
  }

  # Phase 1: diagonal RW, more iterations to navigate from bad start
  cat("  Phase 1 (8,000 iter, diagonal proposal)...\n")
  chain1 <- pmcmc(po_start,
                  Nmcmc     = 8000,
                  Np        = cfg$Np,
                  proposal  = mvn.diag.rw(setNames(
                    rep(0.10, length(params_est)), params_est)))

  # Acceptance rate Phase 1
  acc1 <- chain1@accepts / 8000
  cat(sprintf("  Phase 1 acceptance: %.3f\n", acc1))

  # Posterior covariance from last 4000 iterations of Phase 1
  ph1_df  <- as.data.frame(chain1)
  ph1_sub <- tail(ph1_df[, params_est], 4000)
  vcov1   <- var(ph1_sub) * (2.38^2 / length(params_est))

  # Phase 2: MVN adaptive proposal
  cat("  Phase 2 (MVN adaptive, max 15,000 iter, stop ESS>=200)...\n")
  current_chain <- chain1
  current_po    <- po_start

  phase2_total <- 0
  ess_ok        <- FALSE
  BLOCK         <- 2000
  MAX_ITER      <- 15000

  while (!ess_ok && phase2_total < MAX_ITER) {
    current_po <- pomp(current_po, params = coef(current_chain))
    block <- pmcmc(current_po,
                   Nmcmc    = BLOCK,
                   Np       = cfg$Np,
                   proposal = mvn.rw(vcov1))
    phase2_total <- phase2_total + BLOCK
    combined     <- c(current_chain, block)
    # Keep only phase 2 draws (discard phase 1)
    post_draws   <- tail(as.data.frame(combined)[, params_est], phase2_total)
    ess_now      <- effectiveSize(as.mcmc(post_draws))
    cat(sprintf("    Phase 2 iter %d: ESS = %s\n",
                phase2_total,
                paste(round(ess_now), collapse = " ")))
    ess_ok       <- all(ess_now >= 200)
    current_chain <- combined
  }

  if (!ess_ok) cat("  WARNING: ESS < 200 for some parameters after max iterations\n")

  # Final chain: discard phase 1 burn-in, keep phase 2
  final_draws <- tail(as.data.frame(current_chain)[, params_est], phase2_total)
  final_chain <- as.mcmc(final_draws)

  ess_final <- round(effectiveSize(final_chain))
  cat(sprintf("  Final ESS: %s\n", paste(ess_final, collapse = " ")))

  # Save
  res <- list(
    chain      = final_chain,
    ess        = ess_final,
    start_name = start_name,
    start_vals = start_vals,
    prior      = "flat"
  )
  saveRDS(res, rds_path)
  cat(sprintf("  Saved: %s\n\n", rds_path))
  return(res)
}

# Run all three
results <- list()
for (nm in names(off_starts)) {
  results[[nm]] <- run_off(nm, off_starts[[nm]])
}

# ── 5. Comparison plots ───────────────────────────────────────────────────────
cat(strrep("=", 60), "\n")
cat("Generating comparison plots...\n\n")

# ── 5a. Posterior density comparison ──────────────────────────────────────────
make_density_comparison <- function() {

  # Collect all posterior draws
  all_df <- bind_rows(
    baseline_df[, params_est] |>
      mutate(run = "Baseline\n(tight prior, MLE start)"),
    as.data.frame(results$offA$chain) |>
      mutate(run = "Off-A: βs×0.5, ρs×3\n(flat prior)"),
    as.data.frame(results$offB$chain) |>
      mutate(run = "Off-B: βs×2, ρs×0.3\n(flat prior)"),
    as.data.frame(results$offC$chain) |>
      mutate(run = "Off-C: mixed extreme\n(flat prior)")
  ) |>
    mutate(run = factor(run, levels = c(
      "Baseline\n(tight prior, MLE start)",
      "Off-A: βs×0.5, ρs×3\n(flat prior)",
      "Off-B: βs×2, ρs×0.3\n(flat prior)",
      "Off-C: mixed extreme\n(flat prior)"
    )))

  run_cols <- c(
    "Baseline\n(tight prior, MLE start)"  = "#154273",
    "Off-A: βs×0.5, ρs×3\n(flat prior)"  = "#E05540",
    "Off-B: βs×2, ρs×0.3\n(flat prior)"  = "#0D9488",
    "Off-C: mixed extreme\n(flat prior)"  = "#7C3AED"
  )

  long_df <- all_df |>
    pivot_longer(all_of(params_est), names_to = "parameter", values_to = "value") |>
    mutate(parameter = factor(parameter, levels = params_est))

  # Also mark the MLE starting values as vertical lines
  mle_df <- data.frame(
    parameter = factor(names(mle), levels = params_est),
    value     = as.numeric(mle)
  )

  ggplot(long_df, aes(x = value, colour = run, fill = run)) +
    geom_density(alpha = 0.12, linewidth = 0.8) +
    geom_vline(data = mle_df, aes(xintercept = value),
               linetype = "dashed", colour = "grey40", linewidth = 0.5) +
    facet_wrap(~parameter, scales = "free", nrow = 2,
               labeller = as_labeller(c(
                 Beta1 = "beta[1]", Beta2 = "beta[2]", Beta3 = "beta[3]",
                 rho_y = "rho[y]", rho_o = "rho[o]"
               ), label_parsed)) +
    scale_colour_manual(values = run_cols, name = NULL) +
    scale_fill_manual(values   = run_cols, name = NULL) +
    labs(
      title    = "Wave 1 identifiability check: posterior comparison",
      subtitle = paste(
        "Dashed grey line = MLE starting value.",
        "If all four posteriors overlap substantially,",
        "the likelihood is genuinely identifying the parameters."
      ),
      x = "Parameter value (natural scale)", y = "Density"
    ) +
    theme_bw(base_size = 11) +
    theme(legend.position = "bottom",
          legend.text     = element_text(size = 8),
          strip.text      = element_text(size = 9))
}

p_dens <- make_density_comparison()
ggsave("results/identifiability/wave1_posterior_comparison.pdf",
       p_dens, width = 13, height = 7)
ggsave("results/identifiability/wave1_posterior_comparison.png",
       p_dens, width = 13, height = 7, dpi = 150)
cat("Density comparison saved.\n")

# ── 5b. Posterior summary table ───────────────────────────────────────────────
make_summary_table <- function() {
  all_res <- list(
    Baseline = baseline_df[, params_est],
    OffA     = as.data.frame(results$offA$chain),
    OffB     = as.data.frame(results$offB$chain),
    OffC     = as.data.frame(results$offC$chain)
  )

  rows <- list()
  for (nm in names(all_res)) {
    df <- all_res[[nm]]
    for (p in params_est) {
      rows[[length(rows)+1]] <- data.frame(
        run       = nm,
        parameter = p,
        mean      = mean(df[[p]]),
        sd        = sd(df[[p]]),
        ci_lo     = quantile(df[[p]], 0.025),
        ci_hi     = quantile(df[[p]], 0.975)
      )
    }
  }
  bind_rows(rows) |>
    mutate(across(where(is.numeric), ~round(., 4)))
}

summary_table <- make_summary_table()
write_csv(summary_table,
          "results/identifiability/wave1_identifiability_summary.csv")

# Print to console
cat("\nPosterior summary (mean [95% CI]):\n")
cat(sprintf("%-12s %-8s %12s %12s %12s %12s\n",
            "Run", "Param", "Mean", "SD", "CI_lo", "CI_hi"))
cat(strrep("-", 68), "\n")
for (i in seq_len(nrow(summary_table))) {
  r <- summary_table[i, ]
  cat(sprintf("%-12s %-8s %12.4f %12.4f %12.4f %12.4f\n",
              r$run, r$parameter, r$mean, r$sd, r$ci_lo, r$ci_hi))
}

# ── 5c. Trace plots ───────────────────────────────────────────────────────────
pdf("results/identifiability/wave1_identifiability_traces.pdf",
    width = 14, height = 10)

for (run_name in c("offA","offB","offC")) {
  chain_i <- results[[run_name]]$chain
  par(mfrow = c(length(params_est), 1),
      mar = c(2.5, 4, 2.5, 1),
      oma = c(1, 0, 3, 0))
  for (p in params_est) {
    x        <- as.numeric(chain_i[, p])
    run_mean <- cumsum(x) / seq_along(x)
    plot(x, type = "l", col = "#154273", lwd = 0.35,
         xlab = "", ylab = p, main = p, las = 1)
    lines(run_mean, col = "#e17000", lwd = 1.5)
    abline(h = mean(x), col = "#d73027", lwd = 1.0, lty = 2)
    # Add baseline posterior mean as horizontal green line
    abline(h = mean(baseline_df[[p]]), col = "#059669", lwd = 1.2, lty = 3)
    ess_p <- round(effectiveSize(chain_i[, p, drop = FALSE]))
    mtext(sprintf("ESS=%d", ess_p), side = 3, adj = 1,
          cex = 0.75, col = "grey40", line = 0.1)
  }
  start_nm <- switch(run_name,
    offA = "Off-A: betas x0.5, rhos x3",
    offB = "Off-B: betas x2, rhos x0.3",
    offC = "Off-C: mixed extreme"
  )
  mtext(sprintf("Wave 1 identifiability — %s (flat prior)\n"
                "Navy=trace  Orange=running mean  Red dashed=this posterior mean  "
                "Green dotted=baseline mean",
                start_nm),
        outer = TRUE, cex = 0.85, font = 2, col = "#154273")
}

dev.off()
cat("Trace plots saved.\n")

# ── 6. Diagnostic summary ─────────────────────────────────────────────────────
cat("\n", strrep("=", 60), "\n")
cat("IDENTIFIABILITY CHECK COMPLETE\n\n")
cat("Key question: do the Off posteriors overlap with the Baseline?\n\n")

for (p in params_est) {
  base_mean  <- mean(baseline_df[[p]])
  base_ci    <- quantile(baseline_df[[p]], c(0.025, 0.975))

  cat(sprintf("  %s  baseline: %.4f [%.4f, %.4f]\n", p, base_mean,
              base_ci[1], base_ci[2]))

  for (nm in c("offA","offB","offC")) {
    off_mean <- mean(as.data.frame(results[[nm]]$chain)[[p]])
    in_ci    <- off_mean >= base_ci[1] & off_mean <= base_ci[2]
    flag     <- if (in_ci) "  [within baseline CI]" else "  [OUTSIDE baseline CI]"
    cat(sprintf("      %-8s: %.4f%s\n", nm, off_mean, flag))
  }
  cat("\n")
}

cat("Outputs saved to results/identifiability/\n")
cat("  wave1_posterior_comparison.pdf/.png  — density overlays\n")
cat("  wave1_identifiability_traces.pdf     — trace plots\n")
cat("  wave1_identifiability_summary.csv    — numerical summary\n")
