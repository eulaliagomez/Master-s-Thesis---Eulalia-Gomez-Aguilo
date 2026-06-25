# =============================================================================
# sensitivity_traceplots_wave3.R
#
# Trace plots for all Wave 3 sensitivity runs + baseline.
# Saves to results/diagnostics/wave3_<run_id>/traceplots.pdf
#
# setwd("C:/Users/Usuario/Desktop/CBS - Internship/Analysis/DATA")
# source("sensitivity_traceplots_wave3.R")
# =============================================================================

library(coda)

if (!dir.exists("results/diagnostics")) dir.create("results/diagnostics",
                                                    recursive=TRUE)

runs <- list(
  list(id="baseline",         label="Baseline",          file="results/baseline_wave3.rds",         chain_in_list=FALSE),
  list(id="systematic_f0.10", label="Systematic f=10%",  file="results/wave3_systematic_f0.10.rds", chain_in_list=TRUE),
  list(id="systematic_f0.30", label="Systematic f=30%",  file="results/wave3_systematic_f0.30.rds", chain_in_list=TRUE),
  list(id="systematic_f0.50", label="Systematic f=50%",  file="results/wave3_systematic_f0.50.rds", chain_in_list=TRUE),
  list(id="stochastic_f0.10", label="Stochastic f=10%",  file="results/wave3_stochastic_f0.10.rds", chain_in_list=TRUE),
  list(id="stochastic_f0.30", label="Stochastic f=30%",  file="results/wave3_stochastic_f0.30.rds", chain_in_list=TRUE),
  list(id="stochastic_f0.50", label="Stochastic f=50%",  file="results/wave3_stochastic_f0.50.rds", chain_in_list=TRUE),
  list(id="aggregate_b2",     label="Aggregation 2-week",file="results/wave3_aggregate_b2.rds",     chain_in_list=TRUE),
  list(id="aggregate_b4",     label="Aggregation 4-week",file="results/wave3_aggregate_b4.rds",     chain_in_list=TRUE)
)

param_labels <- c(
  Beta1="beta[1]  (transmission phase 1)",
  Beta2="beta[2]  (transmission phase 2)",
  Beta3="beta[3]  (transmission phase 3)",
  rho_y="rho[y]   (detection rate young)",
  rho_o="rho[o]   (detection rate old)"
)

for (run in runs) {
  cat("Processing:", run$label, "... ")
  if (!file.exists(run$file)) { cat("SKIP\n"); next }
  obj   <- readRDS(run$file)
  chain <- if (run$chain_in_list) obj$chain else obj
  ess   <- round(effectiveSize(chain))
  cat("ESS:", paste(ess, collapse=" "), "\n")
  out_dir <- file.path("results/diagnostics", paste0("wave3_", run$id))
  if (!dir.exists(out_dir)) dir.create(out_dir, recursive=TRUE)
  params   <- colnames(chain)
  n_params <- length(params)
  pdf(file.path(out_dir,"traceplots.pdf"), width=12, height=9)
  par(mfrow=c(n_params,1), mar=c(2.5,4,2.5,1), oma=c(1,0,3,0), cex.main=0.90)
  for (p in params) {
    x        <- as.numeric(chain[,p])
    run_mean <- cumsum(x)/seq_along(x)
    plot(x, type="l", col="#154273", lwd=0.35, xlab="", ylab=p,
         main=param_labels[p], las=1)
    lines(run_mean, col="#e17000", lwd=1.5)
    abline(h=mean(x), col="#d73027", lwd=1.0, lty=2)
    mtext(sprintf("ESS = %d", ess[p]), side=3, adj=1, cex=0.75, col="grey40", line=0.1)
  }
  mtext(sprintf("Wave 3 trace plots: %s", run$label),
        outer=TRUE, cex=1.05, font=2, col="#154273")
  dev.off()
  cat("  Saved:", file.path(out_dir,"traceplots.pdf"), "\n")
}
cat("\nDone. All Wave 3 trace plots in results/diagnostics/wave3_*/\n")
