# =============================================================================
# sensitivity_diagnostics_wave3.R
#
# Posterior densities, pairwise scatter, and PPC for all Wave 3 runs.
# Saves to results/diagnostics/wave3_*/
#
# setwd("C:/Users/Usuario/Desktop/CBS - Internship/Analysis/DATA")
# source("sensitivity_diagnostics_wave3.R")
# =============================================================================

library(tidyverse)
library(pomp)
library(coda)
library(GGally)

SENSITIVITY_FUNCTIONS_ONLY <- TRUE
source("sensitivity_pipeline.R")
rm(SENSITIVITY_FUNCTIONS_ONLY)

if (!dir.exists("results/diagnostics")) dir.create("results/diagnostics",
                                                    recursive=TRUE)

cfg        <- wave_config[["wave3"]]
wave_name  <- "wave3"
params_est <- cfg$params_est
dom_eigen  <- 4.074

series_cols <- c("Young (<60)"="#2166ac",
                 "Old (\u226560)"  ="#d73027",
                 "Total"      ="#1a7a1a")

meas_y <- read_csv(cfg$file_y, show_col_types=FALSE) |>
  rename(reports_young=reports) |> select(week, reports_young)
meas_o <- read_csv(cfg$file_o, show_col_types=FALSE) |>
  rename(reports_old=reports) |> select(week, reports_old)
meas_clean <- left_join(meas_y, meas_o, by="week")

baseline_chain <- readRDS("results/baseline_wave3.rds")
baseline_df    <- as.data.frame(baseline_chain)
cat("Baseline loaded. ESS:", round(effectiveSize(baseline_chain)), "\n\n")

runs <- list(
  list(id="systematic_f0.10", label="Systematic f=10%",  type="systematic", f=0.10, block=NA),
  list(id="systematic_f0.30", label="Systematic f=30%",  type="systematic", f=0.30, block=NA),
  list(id="systematic_f0.50", label="Systematic f=50%",  type="systematic", f=0.50, block=NA),
  list(id="stochastic_f0.10", label="Stochastic f=10%",  type="stochastic", f=0.10, block=NA),
  list(id="stochastic_f0.30", label="Stochastic f=30%",  type="stochastic", f=0.30, block=NA),
  list(id="stochastic_f0.50", label="Stochastic f=50%",  type="stochastic", f=0.50, block=NA),
  list(id="aggregate_b2",     label="Aggregation 2-week",type="aggregate",  f=NA,   block=2),
  list(id="aggregate_b4",     label="Aggregation 4-week",type="aggregate",  f=NA,   block=4)
)

get_degraded_data <- function(run) {
  csv_path <- file.path("results/degraded_data",
                        sprintf("wave3_%s.csv", run$id))
  if (file.exists(csv_path)) return(read_csv(csv_path, show_col_types=FALSE))
  if (run$type == "systematic") return(degrade_systematic(meas_clean, run$f))
  if (run$type == "stochastic") {
    deg_seed <- as.integer(paste0(3, round(run$f*100)))
    return(degrade_stochastic(meas_clean, run$f, seed=deg_seed))
  }
  return(degrade_aggregate(meas_clean, run$block)$data)
}

make_density <- function(post_df_run, run_label) {
  bind_rows(
    post_df_run |> select(all_of(params_est)) |>
      pivot_longer(everything(), names_to="parameter", values_to="value") |>
      mutate(source="Degraded"),
    baseline_df |> select(all_of(params_est)) |>
      pivot_longer(everything(), names_to="parameter", values_to="value") |>
      mutate(source="Baseline")
  ) |>
    mutate(parameter=factor(parameter, levels=params_est)) |>
    ggplot(aes(x=value, fill=source, colour=source)) +
    geom_density(alpha=0.35, linewidth=0.8) +
    facet_wrap(~parameter, scales="free", nrow=2,
               labeller=as_labeller(c(
                 Beta1="beta[1]",Beta2="beta[2]",Beta3="beta[3]",
                 rho_y="rho[y]",rho_o="rho[o]"), label_parsed)) +
    scale_fill_manual(values=c("Baseline"="#2166ac","Degraded"="#d73027"), name=NULL) +
    scale_colour_manual(values=c("Baseline"="#2166ac","Degraded"="#d73027"), name=NULL) +
    labs(title=sprintf("Wave 3 posterior densities: %s", run_label),
         subtitle="Red = degraded | Blue = baseline",
         x="Parameter value", y="Density") +
    theme_bw(base_size=11) +
    theme(legend.position="bottom")
}

make_pairs <- function(post_df_run, run_label) {
  tryCatch(
    ggpairs(post_df_run |> select(all_of(params_est)),
            lower=list(continuous=wrap("points", alpha=0.04, size=0.3, colour="#d73027")),
            upper=list(continuous=wrap("cor", size=3.0)),
            diag =list(continuous=wrap("densityDiag", colour="#d73027",
                                       fill="#d73027", alpha=0.4))) +
      theme_bw(base_size=9) +
      labs(title=sprintf("Wave 3 pairwise posterior: %s", run_label)),
    error=function(e) ggplot() +
      annotate("text",x=0.5,y=0.5,label=paste("Failed:",e$message),size=4) +
      theme_void()
  )
}

make_ppc <- function(po, post_df_run, run_label, n_draws=500) {
  set.seed(42)
  draw_idx <- sample(nrow(post_df_run), min(n_draws, nrow(post_df_run)))
  sims_deg <- lapply(draw_idx, function(i) {
    th <- coef(po); th[params_est] <- as.numeric(post_df_run[i, params_est])
    simulate(po, params=th, nsim=1, format="data.frame") |>
      mutate(draw=i, reports_total=reports_young+reports_old, source="Degraded posterior")
  })
  po_base <- build_pomp(cfg, meas_clean)
  po_base <- pomp(po_base, dprior=make_prior(cfg), paramnames=names(coef(po_base)))
  draw_idx_b <- sample(nrow(baseline_df), min(n_draws, nrow(baseline_df)))
  sims_base <- lapply(draw_idx_b, function(i) {
    th <- coef(po_base); th[params_est] <- as.numeric(baseline_df[i, params_est])
    simulate(po_base, params=th, nsim=1, format="data.frame") |>
      mutate(draw=i, reports_total=reports_young+reports_old, source="Baseline posterior")
  })
  bands <- bind_rows(c(sims_deg, sims_base)) |>
    pivot_longer(c(reports_young,reports_old,reports_total),
                 names_to="series", values_to="cases") |>
    mutate(series=recode(series,"reports_young"="Young (<60)",
                         "reports_old"="Old (\u226560)","reports_total"="Total"),
           series=factor(series, levels=c("Young (<60)","Old (\u226560)","Total"))) |>
    group_by(week, series, source) |>
    summarise(lo=quantile(cases,0.025,na.rm=TRUE), hi=quantile(cases,0.975,na.rm=TRUE),
              med=median(cases,na.rm=TRUE), .groups="drop")
  obs <- meas_clean |>
    mutate(reports_total=reports_young+reports_old) |>
    pivot_longer(c(reports_young,reports_old,reports_total),
                 names_to="series", values_to="cases") |>
    mutate(series=recode(series,"reports_young"="Young (<60)",
                         "reports_old"="Old (\u226560)","reports_total"="Total"),
           series=factor(series, levels=c("Young (<60)","Old (\u226560)","Total")))
  ggplot() +
    geom_ribbon(data=bands |> filter(source=="Baseline posterior"),
                aes(x=week,ymin=lo,ymax=hi), fill="#2166ac", alpha=0.15) +
    geom_line(data=bands |> filter(source=="Baseline posterior"),
              aes(x=week,y=med), colour="#2166ac", linewidth=0.8, linetype="dashed") +
    geom_ribbon(data=bands |> filter(source=="Degraded posterior"),
                aes(x=week,ymin=lo,ymax=hi), fill="#d73027", alpha=0.20) +
    geom_line(data=bands |> filter(source=="Degraded posterior"),
              aes(x=week,y=med), colour="#d73027", linewidth=1.0) +
    geom_point(data=obs, aes(x=week,y=cases), colour="black", size=2, shape=16) +
    facet_wrap(~series, scales="free_y", nrow=3) +
    scale_y_continuous(labels=scales::comma) +
    labs(title=sprintf("Wave 3 PPC: %s", run_label),
         subtitle="Red=degraded | Blue dashed=baseline | Black dots=observed",
         x="Week", y="Weekly reported cases") +
    theme_bw(base_size=11)
}

for (run in runs) {
  cat(strrep("-",60),"\n")
  cat("Processing:", run$label, "\n")
  rds_path <- file.path("results", sprintf("wave3_%s.rds", run$id))
  if (!file.exists(rds_path)) { cat("  SKIP\n\n"); next }
  res         <- readRDS(rds_path)
  post_df_run <- as.data.frame(res$chain)
  cat("  ESS:", paste(round(res$ess), collapse=" "), "\n")
  out_dir <- file.path("results/diagnostics", paste0("wave3_", run$id))
  if (!dir.exists(out_dir)) dir.create(out_dir, recursive=TRUE)
  meas_deg <- get_degraded_data(run)
  if (run$type == "aggregate") {
    agg <- degrade_aggregate(meas_clean, run$block)
    po  <- build_pomp(cfg, agg$data, obs_times=agg$obs_times)
  } else {
    po  <- build_pomp(cfg, meas_deg)
  }
  po <- pomp(po, dprior=make_prior(cfg), paramnames=names(coef(po)))
  cat("  Density..."); p <- make_density(post_df_run, run$label)
  ggsave(file.path(out_dir,"posterior_densities.pdf"), p, width=12, height=7)
  ggsave(file.path(out_dir,"posterior_densities.png"), p, width=12, height=7, dpi=150)
  cat(" done\n")
  cat("  Pairs..."); p <- make_pairs(post_df_run, run$label)
  ggsave(file.path(out_dir,"pairwise_scatter.pdf"), p, width=10, height=10)
  ggsave(file.path(out_dir,"pairwise_scatter.png"), p, width=10, height=10, dpi=150)
  cat(" done\n")
  cat("  PPC..."); p <- make_ppc(po, post_df_run, run$label)
  ggsave(file.path(out_dir,"ppc.pdf"), p, width=10, height=10)
  ggsave(file.path(out_dir,"ppc.png"), p, width=10, height=10, dpi=150)
  cat(" done\n\n")
}

# Baseline
cat(strrep("-",60),"\n"); cat("Processing: Wave 3 Baseline\n")
out_dir_base <- "results/diagnostics/wave3_baseline"
if (!dir.exists(out_dir_base)) dir.create(out_dir_base, recursive=TRUE)
p <- make_density(baseline_df, "Baseline (clean data)")
# replace with clean single-colour version for baseline
p <- baseline_df |> select(all_of(params_est)) |>
  pivot_longer(everything(), names_to="parameter", values_to="value") |>
  mutate(parameter=factor(parameter, levels=params_est)) |>
  ggplot(aes(x=value)) +
  geom_histogram(aes(y=after_stat(density)), bins=50, fill="#2166ac", alpha=0.55) +
  geom_density(colour="black", linewidth=0.8) +
  facet_wrap(~parameter, scales="free", nrow=2,
             labeller=as_labeller(c(Beta1="beta[1]",Beta2="beta[2]",Beta3="beta[3]",
                                    rho_y="rho[y]",rho_o="rho[o]"), label_parsed)) +
  labs(title="Wave 3 posterior densities: Baseline", x="Parameter value", y="Density") +
  theme_bw(base_size=11)
ggsave(file.path(out_dir_base,"posterior_densities.pdf"), p, width=12, height=7)
ggsave(file.path(out_dir_base,"posterior_densities.png"), p, width=12, height=7, dpi=150)
p <- make_pairs(baseline_df, "Baseline")
ggsave(file.path(out_dir_base,"pairwise_scatter.pdf"), p, width=10, height=10)
ggsave(file.path(out_dir_base,"pairwise_scatter.png"), p, width=10, height=10, dpi=150)
cat("  Baseline saved.\n")
cat(strrep("=",60),"\nWAVE 3 DIAGNOSTICS COMPLETE\nOutput: results/diagnostics/wave3_*/\n")
