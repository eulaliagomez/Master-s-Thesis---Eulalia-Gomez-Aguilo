# =============================================================================
# covid19_age_wave1.R  — v3  (corrected parameters + total PPC)
#
# Wave 1: 24 Feb – 5 Jul 2020  |  Ancestral strain
# Two age groups: young (<60) and old (>=60)
#
# Key parameter decisions (all corrected from v1):
#
#   Beta1 = 0.520  → R0 = 2.12 (MLE)
#     In the contact-matrix model, R0 = Beta * dominant_eigenvalue(C/mu_IR)
#     = Beta * 4.074. Dutch Wave 1 R0 ≈ 1.9, not 4.56 (which was the
#     single-pop value absorbing the contact structure).
#
#   t_lock = 4  (not 3)
#     The lockdown was announced 23 March (week 3), but:
#     - incubation period = 5.1 days ≈ 1 week delay before effect visible
#     - compliance built up gradually over 1-2 weeks
#     Effective NPI impact on new infections: week 4-5.
#     Using t_lock=4 and t_open=11 gives correct peak timing (wk5-6).
#
#   Beta2 = 0.120  → Rt = 0.489 (MLE)
#     Slightly more effective suppression than previous version.
#     Consistent with RIVM Rt ≈ 0.5-0.6 at height of lockdown.
#
#   Beta3 = 0.180  → Rt = 0.733 (MLE)
#     Gradual reopening kept Rt below 1 but above lockdown level.
#     RIVM Rt ≈ 0.7-0.8 during May-June 2020.
#
#   t_open = 11  (not 12)
#     Phase 1 reopening 11 May = week 11 of wave. One week earlier than
#     previous version; consistent with policy date.
#
#   rho_y = 0.060  (MLE)
#     Wave 1 testing was hospital/clinical criteria only. Young,
#     mostly-asymptomatic cases were largely undetected.
#
#   rho_o = 0.300  (MLE)
#     Care home residents were tested comprehensively (institutional
#     setting). Dutch nursing home outbreaks had near-total ascertainment.
#     Prior: logit-N(logit(0.35), 0.60²). Source: Dutca et al. 2020;
#     van den Wijngaard et al. 2021 Lancet Reg Health Europe.
#
#   k = 10  FIXED  (Endo et al. 2020; not identifiable from T=19 obs)
#
# Contact matrix: Prem et al. 2021 NLD, all settings, 16x16 → 2x2 balanced
#   C_yy=9.4901, C_yo=0.7523, C_oy=2.8524, C_oo=3.2609
#   Source: https://doi.org/10.1371/journal.pcbi.1009098
# =============================================================================

library(tidyverse)
library(pomp)
library(coda)
library(GGally)

set.seed(123456)

# ── 1. DATA ──────────────────────────────────────────────────────────────────
meas_y <- read_csv("covid19_wave1_young_NL.csv") |>
  rename(reports_young = reports) |> select(week, reports_young)
meas_o <- read_csv("covid19_wave1_old_NL.csv") |>
  rename(reports_old = reports) |> select(week, reports_old)
meas_w1 <- left_join(meas_y, meas_o, by = "week")

cat("Wave 1 data:\n")
cat(sprintf("  T = %d weeks\n", nrow(meas_w1)))
cat(sprintf("  Young: peak=%d (wk%d), total=%d\n",
            max(meas_w1$reports_young), which.max(meas_w1$reports_young),
            sum(meas_w1$reports_young)))
cat(sprintf("  Old:   peak=%d (wk%d), total=%d\n",
            max(meas_w1$reports_old), which.max(meas_w1$reports_old),
            sum(meas_w1$reports_old)))

# ── 2. CONTACT MATRIX & POPULATION ───────────────────────────────────────────
C_yy <- 9.4901; C_yo <- 0.7523; C_oy <- 2.8524; C_oo <- 3.2609
N_y  <- 12736000; N_o <- 3359000

# ── 3. NPI COVARIATE ─────────────────────────────────────────────────────────
# Phase 0 (t < 4):  free transmission — Beta1
# Phase 1 (4≤t<11): lockdown — Beta2
# Phase 2 (t≥11):   gradual reopen — Beta3
#
# t_lock=4: lockdown announced wk3 but incubation lag means impact visible wk4
# t_open=11: Phase 1 reopening 11 May = wk11 from 24 Feb start

npi_df    <- data.frame(week      = 0:22,
                        npi_phase = c(rep(0, 4), rep(1, 6), rep(2, 13)))
npi_covar <- covariate_table(npi_df, times = "week")

# ── 4. PROCESS MODEL ─────────────────────────────────────────────────────────
seir_age_step <- Csnippet("
  double eff_beta;
  if      (npi_phase < 0.5) eff_beta = Beta1;
  else if (npi_phase < 1.5) eff_beta = Beta2;
  else                       eff_beta = Beta3;
 
  double lam_y = eff_beta * dt * (C_yy * I_y/N_y + C_yo * I_o/N_o);
  double lam_o = eff_beta * dt * (C_oy * I_y/N_y + C_oo * I_o/N_o);
 
  double dN_SE_y = rbinom(S_y, 1.0 - exp(-lam_y));
  double dN_EI_y = rbinom(E_y, 1.0 - exp(-mu_EI * dt));
  double dN_IR_y = rbinom(I_y, 1.0 - exp(-mu_IR * dt));
  double dN_SE_o = rbinom(S_o, 1.0 - exp(-lam_o));
  double dN_EI_o = rbinom(E_o, 1.0 - exp(-mu_EI * dt));
  double dN_IR_o = rbinom(I_o, 1.0 - exp(-mu_IR * dt));
 
  S_y -= dN_SE_y; E_y += dN_SE_y-dN_EI_y; I_y += dN_EI_y-dN_IR_y; R_y += dN_IR_y; H_y += dN_IR_y;
  S_o -= dN_SE_o; E_o += dN_SE_o-dN_EI_o; I_o += dN_EI_o-dN_IR_o; R_o += dN_IR_o; H_o += dN_IR_o;
")

# ── 5. OBSERVATION MODEL ──────────────────────────────────────────────────────
covid_dmeas <- Csnippet("
  lik = dnbinom_mu(reports_young, k, rho_y * H_y, 1)
      + dnbinom_mu(reports_old,   k, rho_o * H_o, 1);
  if (!give_log) lik = exp(lik);
")
covid_rmeas <- Csnippet("
  reports_young = rnbinom_mu(k, rho_y * H_y);
  reports_old   = rnbinom_mu(k, rho_o * H_o);
")

# ── 6. INITIAL CONDITIONS ─────────────────────────────────────────────────────
# I0_y = 334  / (0.060 * 2.41) = 2,310;  E0_y = 2,310 * 2.41/1.37 = 4,063
# I0_o = 202  / (0.300 * 2.41) =   279;  E0_o =   279 * 2.41/1.37 =   491

seir_rinit <- Csnippet("
  double I0_y = 2310.0;
  double E0_y = 4063.0;
  double I0_o =  279.0;
  double E0_o =  491.0;
  S_y = nearbyint(eta_y * N_y) - nearbyint(I0_y) - nearbyint(E0_y);
  E_y = nearbyint(E0_y); I_y = nearbyint(I0_y);
  R_y = nearbyint((1.0 - eta_y) * N_y); H_y = 0.0;
  S_o = nearbyint(eta_o * N_o) - nearbyint(I0_o) - nearbyint(E0_o);
  E_o = nearbyint(E0_o); I_o = nearbyint(I0_o);
  R_o = nearbyint((1.0 - eta_o) * N_o); H_o = 0.0;
")

# ── 7. BUILD pomp OBJECT ─────────────────────────────────────────────────────
params_est <- c("Beta1", "Beta2", "Beta3", "rho_y", "rho_o")

covid_w1 <- meas_w1 |>
  pomp(
    times      = "week",
    t0         = 0,
    rprocess   = euler(seir_age_step, delta.t = 1/7),
    rinit      = seir_rinit,
    rmeasure   = covid_rmeas,
    dmeasure   = covid_dmeas,
    covar      = npi_covar,
    partrans   = parameter_trans(
      log   = c("Beta1", "Beta2", "Beta3"),
      logit = c("rho_y", "rho_o")
    ),
    paramnames = c("Beta1","Beta2","Beta3","mu_EI","mu_IR",
                   "eta_y","eta_o","rho_y","rho_o","k",
                   "N_y","N_o","C_yy","C_yo","C_oy","C_oo"),
    statenames = c("S_y","E_y","I_y","R_y","H_y",
                   "S_o","E_o","I_o","R_o","H_o"),
    accumvars  = c("H_y","H_o"),
    obsnames   = c("reports_young","reports_old")
  )

# ── 8. STARTING PARAMETERS ────────────────────────────────────────────────────
theta_start <- c(
  Beta1 = 0.520,  # R0 = 2.12  (MLE; R0 = Beta * 4.074)
  Beta2 = 0.120,  # Rt = 0.489 (MLE; RIVM Rt ≈ 0.5 at lockdown peak)
  Beta3 = 0.180,  # Rt = 0.733 (MLE; gradual reopening)
  mu_EI = 1.37,   # FIXED: Lauer 2020
  mu_IR = 2.41,   # FIXED: He 2020
  eta_y = 0.99,   # FIXED: full susceptibility young
  eta_o = 0.99,   # FIXED: full susceptibility old (Wave 1 = first epidemic)
  rho_y = 0.060,  # hospital-era testing; young largely untested (MLE)
  rho_o = 0.300,  # care home/hospital; old near-fully ascertained (MLE)
  k     = 10,     # FIXED: Endo et al. 2020
  N_y   = N_y, N_o = N_o,
  C_yy  = C_yy, C_yo = C_yo, C_oy = C_oy, C_oo = C_oo
)
coef(covid_w1) <- theta_start

# ── 9. PRIOR ─────────────────────────────────────────────────────────────────
# log(0.460) = -0.7765; log(0.140) = -1.9661; log(0.170) = -1.7720
# logit(0.050) = -2.9444; logit(0.350) = -0.6190
# Wider priors than single-pop model (SD=0.50 for betas, 0.70 for rhos)
# to allow PMCMC to compensate for structural model misfit

covid_dprior <- Csnippet("
  double lB1   = log(Beta1);
  double lB2   = log(Beta2);
  double lB3   = log(Beta3);
  double lrho_y = log(rho_y / (1.0 - rho_y));
  double lrho_o = log(rho_o / (1.0 - rho_o));
 
  lik = dnorm(lB1,   -0.654, 0.50, 1)
      + dnorm(lB2,   -2.120, 0.55, 1)
      + dnorm(lB3,   -1.715, 0.55, 1)
      + dnorm(lrho_y,-2.752, 0.70, 1)
      + dnorm(lrho_o,-0.847, 0.70, 1);
 
  if (!give_log) lik = exp(lik);
")
covid_w1 <- pomp(covid_w1, dprior = covid_dprior,
                 paramnames = names(theta_start))

# ── 10. SIMULATION CHECK ──────────────────────────────────────────────────────
# What to look for:
#   - Peak timing: young wk5-6, old wk5
#   - Scale: young ~3000-4000, old ~3000-4000 at peak
#   - Gradual decline (not sharp cliff)
#   - Tail holding above zero through wk19

set.seed(2024)
sim_df <- covid_w1 |>
  simulate(nsim = 50, format = "data.frame", include.data = TRUE)

# Add total
sim_df <- sim_df |>
  mutate(reports_total = reports_young + reports_old)

p_sim <- sim_df |>
  pivot_longer(c(reports_young, reports_old, reports_total),
               names_to = "series", values_to = "cases") |>
  mutate(series = recode(series,
                         "reports_young" = "Young (<60)",
                         "reports_old"   = "Old (≥60)",
                         "reports_total" = "Total")) |>
  mutate(series = factor(series, levels = c("Young (<60)","Old (≥60)","Total"))) |>
  ggplot(aes(x = week, y = cases, group = .id,
             colour    = (.id == "data"),
             linewidth = (.id == "data"),
             alpha     = (.id == "data"))) +
  geom_line() +
  facet_wrap(~series, scales = "free_y", nrow = 3) +
  scale_colour_manual(values = c("TRUE"="black","FALSE"="#2166ac"), guide="none") +
  scale_linewidth_manual(values = c("TRUE"=1.0,"FALSE"=0.3), guide="none") +
  scale_alpha_manual(values = c("TRUE"=1.0,"FALSE"=0.4), guide="none") +
  scale_y_continuous(labels = scales::comma) +
  labs(title    = "Wave 1 (age-structured): Simulation check",
       subtitle = "Check: bracket data in scale and timing; all three panels",
       x = "Week", y = "Weekly reported cases") +
  theme_bw(base_size = 10)
print(p_sim)

# ── 11. PARTICLE FILTER DIAGNOSTIC ───────────────────────────────────────────
Np_use <- 3000
set.seed(999)
ll_check <- replicate(10, logLik(pfilter(covid_w1, Np = Np_use)))
ll_check  <- ll_check[is.finite(ll_check)]
cat(sprintf("\npfilter (Np=%d): mean=%.2f  SD=%.3f\n",
            Np_use, mean(ll_check), sd(ll_check)))
if (sd(ll_check) > 0.5) {
  cat("  CAUTION: SD > 0.5 — increase Np_use\n")
  Np_use <- Np_use * 2
  cat(sprintf("  Retrying with Np=%d\n", Np_use))
  ll2 <- replicate(10, logLik(pfilter(covid_w1, Np = Np_use)))
  cat(sprintf("  New SD=%.3f\n", sd(ll2[is.finite(ll2)])))
} else {
  cat("  OK — proceed to PMCMC\n")
}

# ── 12. PHASE 1 PMCMC ────────────────────────────────────────────────────────
pmcmc_p1 <- pmcmc(
  covid_w1,
  Nmcmc    = 5000,
  Np       = Np_use,
  proposal = mvn_diag_rw(c(Beta1=0.15, Beta2=0.15, Beta3=0.15,
                           rho_y=0.15, rho_o=0.15)^2)
)
acc_p1 <- pmcmc_p1@accepts / 5000
cat(sprintf("Phase 1 acceptance rate: %.3f\n", acc_p1))

# ── 13. PHASE 2 PMCMC ────────────────────────────────────────────────────────
chain_p1  <- as.matrix(as.data.frame(traces(pmcmc_p1)))
post_cov  <- cov(chain_p1[2001:5000, params_est])
opt_scale <- (2.38^2) / length(params_est)

pmcmc_p2 <- pmcmc(
  pmcmc_p1,
  Nmcmc    = 15000,
  Np       = Np_use,
  proposal = mvn_rw(opt_scale * post_cov)
)
acc_p2 <- pmcmc_p2@accepts / 15000
cat(sprintf("Phase 2 acceptance rate: %.3f  (target 0.20-0.40)\n", acc_p2))

# ── 14. PRODUCTION CHAIN ─────────────────────────────────────────────────────
chain_p2   <- as.matrix(as.data.frame(traces(pmcmc_p2)))
prod_chain <- as.mcmc(chain_p2[2001:nrow(chain_p2), params_est])

cat("\nESS:\n"); print(round(effectiveSize(prod_chain)))
plot(prod_chain, ask = FALSE)

# ── 15. POSTERIOR SUMMARIES ───────────────────────────────────────────────────
dom_eigen <- 4.074
post_df <- as.data.frame(prod_chain) |>
  mutate(
    R0          = Beta1 * dom_eigen,
    Rt_lockdown = Beta2 * dom_eigen,
    Rt_reopen   = Beta3 * dom_eigen
  )

cat("\nPosterior means:\n")
print(round(colMeans(post_df[, params_est]), 5))
cat("\n95% CIs:\n")
print(round(apply(post_df[, params_est], 2, quantile, c(0.025, 0.975)), 5))
cat("\nDerived R values:\n")
for (v in c("R0","Rt_lockdown","Rt_reopen"))
  cat(sprintf("  %s: %.2f [%.2f, %.2f]\n", v,
              mean(post_df[[v]]),
              quantile(post_df[[v]], 0.025),
              quantile(post_df[[v]], 0.975)))

# ── 16. POSTERIOR DENSITY PLOTS ──────────────────────────────────────────────
post_df |>
  select(all_of(params_est)) |>
  pivot_longer(everything(), names_to = "parameter", values_to = "value") |>
  ggplot(aes(x = value)) +
  geom_histogram(aes(y = after_stat(density)), bins=50,
                 fill="#2166ac", alpha=0.55) +
  geom_density(colour="black", linewidth=0.8) +
  facet_wrap(~parameter, scales="free", nrow=2) +
  labs(title="Wave 1 (age-structured): Posterior marginal distributions",
       x="Parameter value", y="Density") +
  theme_bw(base_size=11)

# ── 17. PAIRWISE SCATTER ──────────────────────────────────────────────────────
ggpairs(post_df |> select(all_of(params_est)),
        lower = list(continuous = wrap("points", alpha=0.04, size=0.4)),
        upper = list(continuous = wrap("cor", size=3)),
        diag  = list(continuous = wrap("densityDiag"))) +
  theme_bw(base_size=9) +
  labs(title="Wave 1 (age-structured): Pairwise posterior scatter")

# ── 18. POSTERIOR PREDICTIVE CHECK — ALL THREE PANELS ─────────────────────────
# Generates 95% predictive intervals for:
#   (a) young cases only
#   (b) old cases only
#   (c) total = young + old  ← NEW
# This allows checking both the age breakdown and the overall epidemic fit.

set.seed(42)
pp_sims <- lapply(sample(nrow(post_df), 500), function(i) {
  th            <- coef(covid_w1)
  th[params_est] <- as.numeric(post_df[i, params_est])
  simulate(covid_w1, params=th, nsim=1, format="data.frame") |>
    mutate(draw=i,
           reports_total = reports_young + reports_old)
})

pp_long <- bind_rows(pp_sims) |>
  pivot_longer(c(reports_young, reports_old, reports_total),
               names_to="series", values_to="cases") |>
  mutate(series = recode(series,
                         "reports_young" = "Young (<60)",
                         "reports_old"   = "Old (≥60)",
                         "reports_total" = "Total"))

pp_band <- pp_long |>
  group_by(week, series) |>
  summarise(lo  = quantile(cases, 0.025, na.rm=TRUE),
            hi  = quantile(cases, 0.975, na.rm=TRUE),
            med = median(cases,   na.rm=TRUE),
            .groups="drop") |>
  mutate(series = factor(series, levels=c("Young (<60)","Old (≥60)","Total")))

obs_long <- meas_w1 |>
  mutate(reports_total = reports_young + reports_old) |>
  pivot_longer(c(reports_young, reports_old, reports_total),
               names_to="series", values_to="cases") |>
  mutate(series = recode(series,
                         "reports_young" = "Young (<60)",
                         "reports_old"   = "Old (≥60)",
                         "reports_total" = "Total"),
         series = factor(series, levels=c("Young (<60)","Old (≥60)","Total")))

# Colour by series
series_cols <- c("Young (<60)"="#2166ac",
                 "Old (≥60)"  ="#d73027",
                 "Total"      ="#1a7a1a")

ggplot() +
  geom_ribbon(data=pp_band,
              aes(x=week, ymin=lo, ymax=hi, fill=series),
              alpha=0.25) +
  geom_line(data=pp_band,
            aes(x=week, y=med, colour=series),
            linewidth=1) +
  geom_point(data=obs_long,
             aes(x=week, y=cases, colour=series),
             size=2, shape=16) +
  facet_wrap(~series, scales="free_y", nrow=3) +
  scale_colour_manual(values=series_cols, guide="none") +
  scale_fill_manual(values=series_cols, guide="none") +
  scale_y_continuous(labels=scales::comma) +
  labs(title    = "Wave 1 (age-structured): Posterior predictive check",
       subtitle = "Shaded=95% interval; line=posterior median; dots=data",
       x="Week", y="Weekly reported cases") +
  theme_bw(base_size=11)

# ── Save baseline posterior for sensitivity analysis ──────────────────────────
dir.create("results", showWarnings = FALSE)
saveRDS(prod_chain, "results/baseline_wave1.rds")
cat("Baseline saved to results/baseline_wave1.rds\n")

source("run_wave1_sensitivity.R")

# ── Generate and save the plots ──────────────────────────
setwd("C:/Users/Usuario/Desktop/CBS - Internship/Analysis/DATA")
source("sensitivity_plots_wave1_fix2.R")
source("sensitivity_diagnostics_wave1.R")
source("sensitivity_traceplots_wave1.R")
