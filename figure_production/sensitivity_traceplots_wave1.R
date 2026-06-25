# =============================================================================
# sensitivity_traceplots_wave1.R
#
# Trace plots only (no density) for all Wave 1 sensitivity runs + baseline.
# One PDF per run saved in results/diagnostics/<run_id>/traceplots.pdf
#
# setwd("C:/Users/Usuario/Desktop/CBS - Internship/Analysis/DATA")
# source("sensitivity_traceplots_wave1.R")
# =============================================================================

library(coda)

if (!dir.exists("results/diagnostics")) dir.create("results/diagnostics",
                                                    recursive = TRUE)

runs <- list(
  list(id = "baseline",         label = "Baseline (clean data)",  file = "results/baseline_wave1.rds",              chain_in_list = FALSE),
  list(id = "systematic_f0.10", label = "Systematic f=10%",       file = "results/wave1_systematic_f0.10.rds",      chain_in_list = TRUE),
  list(id = "systematic_f0.30", label = "Systematic f=30%",       file = "results/wave1_systematic_f0.30.rds",      chain_in_list = TRUE),
  list(id = "systematic_f0.50", label = "Systematic f=50%",       file = "results/wave1_systematic_f0.50.rds",      chain_in_list = TRUE),
  list(id = "stochastic_f0.10", label = "Stochastic f=10%",       file = "results/wave1_stochastic_f0.10.rds",      chain_in_list = TRUE),
  list(id = "stochastic_f0.30", label = "Stochastic f=30%",       file = "results/wave1_stochastic_f0.30.rds",      chain_in_list = TRUE),
  list(id = "stochastic_f0.50", label = "Stochastic f=50%",       file = "results/wave1_stochastic_f0.50.rds",      chain_in_list = TRUE),
  list(id = "aggregate_b2",     label = "Aggregation 2-week",     file = "results/wave1_aggregate_b2.rds",          chain_in_list = TRUE),
  list(id = "aggregate_b4",     label = "Aggregation 4-week",     file = "results/wave1_aggregate_b4.rds",          chain_in_list = TRUE)
)

for (run in runs) {

  cat("Processing:", run$label, "... ")

  if (!file.exists(run$file)) {
    cat("SKIP (file not found)\n")
    next
  }

  # Load chain
  obj   <- readRDS(run$file)
  chain <- if (run$chain_in_list) obj$chain else obj

  # ESS
  ess <- round(effectiveSize(chain))
  cat("ESS:", paste(ess, collapse = " "), "\n")

  # Output directory
  out_dir <- file.path("results/diagnostics", run$id)
  if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

  params  <- colnames(chain)
  n_params <- length(params)

  # Nice parameter labels for plot titles
  param_labels <- c(
    Beta1 = "beta[1]  (transmission phase 1)",
    Beta2 = "beta[2]  (transmission phase 2)",
    Beta3 = "beta[3]  (transmission phase 3)",
    rho_y = "rho[y]   (detection rate young)",
    rho_o = "rho[o]   (detection rate old)"
  )

  pdf(file.path(out_dir, "traceplots.pdf"), width = 12, height = 9)

  par(mfrow  = c(n_params, 1),
      mar    = c(2.5, 4, 2.5, 1),
      oma    = c(1, 0, 3, 0),
      cex.main = 0.95)

  for (p in params) {
    x <- as.numeric(chain[, p])
    n <- length(x)

    # Running mean for convergence check
    run_mean <- cumsum(x) / seq_along(x)

    # Trace
    plot(x,
         type = "l",
         col  = "#154273",
         lwd  = 0.35,
         xlab = "",
         ylab = p,
         main = param_labels[p],
         las  = 1)

    # Running mean overlaid in orange
    lines(run_mean, col = "#e17000", lwd = 1.5)

    # ESS annotation
    mtext(sprintf("ESS = %d", ess[p]),
          side = 3, adj = 1, cex = 0.75, col = "grey40", line = 0.1)

    # Horizontal line at posterior mean
    abline(h = mean(x), col = "#d73027", lwd = 1.0, lty = 2)
  }

  # Main title across all panels
  mtext(sprintf("Wave 1 trace plots: %s", run$label),
        outer = TRUE, cex = 1.1, font = 2, col = "#154273")

  dev.off()
  cat("  Saved:", file.path(out_dir, "traceplots.pdf"), "\n")
}

cat("\nDone. All trace plots saved in results/diagnostics/\n")
cat("Each PDF contains one panel per parameter with:\n")
cat("  Navy line   = MCMC trace\n")
cat("  Orange line = running mean (should stabilise)\n")
cat("  Red dashed  = posterior mean\n")
cat("  ESS         = effective sample size (top right)\n")
