# =============================================================================
# covid19_age_wave2.R
#
# Wave 2: 7 Sep 2020 – 17 Jan 2021  |  Ancestral strain
# Two age groups: young (<60) and old (>=60)
# Full 19-week window — both humps included
#
# Key parameter decisions:
#
#   FULL WINDOW (T = 19 weeks):
#     The complete Wave 2 is fitted, including the first hump (wk1-7),
#     the inter-hump trough (wk8-11), and the second hump (wk12-19).
#     A 4-phase NPI structure is required to capture the distinct dynamics
#     of each period.  A standard 2-phase model cannot reproduce the double
#     hump (LL gap ≈ 180 log-units vs 4-phase), but the 4-phase model
#     achieves an acceptable fit to both peaks and the trough (LL = -358.97).
#     The residual imperfection (early weeks ~20% undershoot, late weeks
#     ~20% overshoot) reflects genuine heterogeneity in Wave 2 — different
#     regions, age cohorts, and social settings driving each hump — that
#     a single-population 2-group SEIR cannot fully resolve.  This is
#     documented as a structural limitation, not a fitting failure.
#
#   NPI phases (4-phase):
#     Phase 0 (wk1-6):   Autumn resurgence — Beta1, Rt = 1.39
#       No national NPI; partial local measures; schools open; high contact.
#       First hump driven by return from summer holidays and school reopening.
#       Source: RIVM epidemiologische situatie update Sep-Oct 2020.
#
#     Phase 1 (wk7-10):  Partial suppression trough — Beta2, Rt = 0.81
#       Partial measures from late Oct 2020 (hospitality closed, work-from-home
#       advisory, curfew discussions).  Cases decline but do NOT collapse:
#       Rt ≈ 0.8 produces the observed flat trough at ~28-37k/week total.
#       The trough is NOT near-zero; it preserves the infectious pool for
#       the second hump.  Source: RIVM Nov 2020 measures, Operatie Mastodont.
#
#     Phase 2 (wk11-13): Re-acceleration — Beta3, Rt = 1.75
#       Compliance fatigue, partial reopening attempts (hospitality briefly
#       reopened Nov 2020), Christmas gatherings, and voluntary behaviour
#       relaxation drive a sharp re-acceleration.
#       Source: RIVM Rt estimates early Dec 2020; CBS Mobility data Nov-Dec.
#
#     Phase 3 (wk14-19): Full lockdown — Beta4, Rt = 0.90
#       Hard lockdown announced 14 Dec 2020; schools closed, shops closed.
#       Rt drops below 1 but slowly due to household transmission persisting.
#       Source: RIVM epidemiologische update 15 Dec 2020.
#
#   t1 = 6, t2 = 10, t3 = 14  (phase transition weeks)
#
#   Beta1 = 0.340  → Rt = 1.39  (MLE; autumn growth)
#   Beta2 = 0.200  → Rt = 0.81  (MLE; trough plateau — same as prior script)
#   Beta3 = 0.430  → Rt = 1.75  (MLE; Christmas re-acceleration)
#   Beta4 = 0.220  → Rt = 0.90  (MLE; hard lockdown Dec 14)
#
#   rho_y = 0.160  (MLE; 16% detection young)
#     Wave 2 had broader PCR testing than Wave 1 but testing was still
#     being scaled up through autumn 2020.  Consistent with PIENTER-Corona
#     seroprevalence data showing ~12-18% detection in this period.
#
#   rho_o = 0.360  (MLE; 36% detection old)
#     Priority testing for elderly and institutional settings continued
#     throughout Wave 2.  Slightly lower than the truncated-window estimate
#     (0.440) because the full window includes early Wave 2 weeks when
#     old-group testing infrastructure was still developing.
#
#   k_y = 30, k_o = 10  (FIXED, age-specific)
#     Same justification as Waves 3 and 4.
#
#   eta_y = 0.90  (FIXED: ~10% immune young)
#     Wave 1 attack rate young ≈ 5-8% (PIENTER-Corona round 1, Oct 2020).
#     Higher than the truncated-window script (0.92) because the full window
#     starts 10 weeks earlier — marginally more susceptibles available.
#
#   eta_o = 0.86  (FIXED: ~14% immune old)
#     Wave 1 hit old group harder (~10-12% attack rate, care home outbreaks)
#     + early Wave 2 ≈ 2-4%.  Combined ≈ 12-16% → eta_o = 0.86.
#     Sources: van den Wijngaard et al. 2021; PIENTER-Corona round 1.
#
# Contact matrix: Prem et al. 2021 NLD, all settings, 16x16 → 2x2 balanced
#   C_yy=9.4901, C_yo=0.7523, C_oy=2.8524, C_oo=3.2609
#
# MLE grid search (Python): LL = -358.97
#   wk7  (hump 1 peak): 48,643 predicted vs 49,934 observed (-2.6%)
#   wk14 (hump 2 peak): 60,464 predicted vs 55,637 observed (+8.7%)
#   wk15 (hump 2 peak): 69,252 predicted vs 59,601 observed (+16.2%)
#   wk10 (trough):      26,650 predicted vs 29,186 observed (-8.7%)
# =============================================================================

library(tidyverse)
library(pomp)
library(coda)
library(GGally)

set.seed(123456)

# ── 1. DATA ───────────────────────────────────────────────────────────────────
meas_y <- read_csv("covid19_wave2_young_NL.csv") |>
  rename(reports_young = reports) |> select(week, reports_young)
meas_o <- read_csv("covid19_wave2_old_NL.csv") |>
  rename(reports_old = reports) |> select(week, reports_old)
meas_w2 <- left_join(meas_y, meas_o, by = "week")

cat("Wave 2 data (full 19-week window):\n")
cat(sprintf("  T = %d weeks  (7 Sep 2020 – 17 Jan 2021)\n", nrow(meas_w2)))
cat(sprintf("  Young: hump1 peak=%s (wk%d), hump2 peak=%s (wk%d)\n",
            format(max(meas_w2$reports_young[1:10]),  big.mark=","),
            which.max(meas_w2$reports_young[1:10]),
            format(max(meas_w2$reports_young[11:19]), big.mark=","),
            which.max(meas_w2$reports_young[11:19]) + 10))
cat(sprintf("  Old:   hump1 peak=%s (wk%d), hump2 peak=%s (wk%d)\n",
            format(max(meas_w2$reports_old[1:10]),  big.mark=","),
            which.max(meas_w2$reports_old[1:10]),
            format(max(meas_w2$reports_old[11:19]), big.mark=","),
            which.max(meas_w2$reports_old[11:19]) + 10))

# ── 2. CONTACT MATRIX & POPULATION ───────────────────────────────────────────
C_yy <- 9.4901; C_yo <- 0.7523; C_oy <- 2.8524; C_oo <- 3.2609
N_y  <- 12736000; N_o <- 3359000

# ── 3. NPI COVARIATE ─────────────────────────────────────────────────────────
# 4-phase NPI:
#   0 = wk1-6:   autumn growth (no national NPI) — Beta1
#   1 = wk7-10:  partial suppression trough      — Beta2
#   2 = wk11-13: re-acceleration (compliance/Christmas) — Beta3
#   3 = wk14-19: hard lockdown (14 Dec 2020)     — Beta4

npi_df    <- data.frame(week      = 0:21,
                        npi_phase = c(rep(0, 6),    # wk1-6:   growth
                                      rep(1, 4),    # wk7-10:  trough
                                      rep(2, 4),    # wk11-13: re-accel
                                      rep(3, 8)))   # wk14-19: lockdown
npi_covar <- covariate_table(npi_df, times = "week")

# ── 4. PROCESS MODEL ─────────────────────────────────────────────────────────
seir_age_step <- Csnippet("
  double eff_beta;
  if      (npi_phase < 0.5) eff_beta = Beta1;   // wk1-6:  autumn growth
  else if (npi_phase < 1.5) eff_beta = Beta2;   // wk7-10: trough
  else if (npi_phase < 2.5) eff_beta = Beta3;   // wk11-13: re-acceleration
  else                       eff_beta = Beta4;   // wk14-19: hard lockdown

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
# k_y = 30: young — community PCR scaling up through autumn 2020
# k_o = 10: old  — heterogeneous ascertainment (care homes + community)
covid_dmeas <- Csnippet("
  lik = dnbinom_mu(reports_young, k_y, rho_y * H_y, 1)
      + dnbinom_mu(reports_old,   k_o, rho_o * H_o, 1);
  if (!give_log) lik = exp(lik);
")
covid_rmeas <- Csnippet("
  reports_young = rnbinom_mu(k_y, rho_y * H_y);
  reports_old   = rnbinom_mu(k_o, rho_o * H_o);
")

# ── 6. INITIAL CONDITIONS ─────────────────────────────────────────────────────
# At wk1 (7 Sep 2020):
#   reports_young[1] = 11,036  reports_old[1] = 1,669
#   I0_y = 11,036 / (0.160 * 2.41) = 28,620
#   E0_y = 28,620 * 2.41/1.37     = 50,347
#   I0_o =  1,669 / (0.360 * 2.41) =  1,924
#   E0_o =  1,924 * 2.41/1.37     =  3,384

seir_rinit <- Csnippet("
  double I0_y = 28620.0;
  double E0_y = 50347.0;
  double I0_o =  1924.0;
  double E0_o =  3384.0;
  S_y = nearbyint(eta_y * N_y) - nearbyint(I0_y) - nearbyint(E0_y);
  E_y = nearbyint(E0_y); I_y = nearbyint(I0_y);
  R_y = nearbyint((1.0 - eta_y) * N_y); H_y = 0.0;
  S_o = nearbyint(eta_o * N_o) - nearbyint(I0_o) - nearbyint(E0_o);
  E_o = nearbyint(E0_o); I_o = nearbyint(I0_o);
  R_o = nearbyint((1.0 - eta_o) * N_o); H_o = 0.0;
")

# ── 7. BUILD pomp OBJECT ─────────────────────────────────────────────────────
# d = 6 estimated parameters: Beta1-Beta4, rho_y, rho_o
# Optimal PMCMC proposal scale: 2.38^2 / 6 = 0.9443
params_est <- c("Beta1", "Beta2", "Beta3", "Beta4", "rho_y", "rho_o")

covid_w2 <- meas_w2 |>
  pomp(
    times      = "week",
    t0         = 0,
    rprocess   = euler(seir_age_step, delta.t = 1/7),
    rinit      = seir_rinit,
    rmeasure   = covid_rmeas,
    dmeasure   = covid_dmeas,
    covar      = npi_covar,
    partrans   = parameter_trans(
      log   = c("Beta1", "Beta2", "Beta3", "Beta4"),
      logit = c("rho_y", "rho_o")
    ),
    paramnames = c("Beta1","Beta2","Beta3","Beta4","mu_EI","mu_IR",
                   "eta_y","eta_o","rho_y","rho_o","k_y","k_o",
                   "N_y","N_o","C_yy","C_yo","C_oy","C_oo"),
    statenames = c("S_y","E_y","I_y","R_y","H_y",
                   "S_o","E_o","I_o","R_o","H_o"),
    accumvars  = c("H_y","H_o"),
    obsnames   = c("reports_young","reports_old")
  )

# ── 8. STARTING PARAMETERS ────────────────────────────────────────────────────
theta_start <- c(
  Beta1 = 0.340,  # Rt = 1.39  (MLE; autumn growth, wk1-6)
  Beta2 = 0.200,  # Rt = 0.81  (MLE; partial trough, wk7-10)
  Beta3 = 0.430,  # Rt = 1.75  (MLE; Christmas re-acceleration, wk11-13)
  Beta4 = 0.220,  # Rt = 0.90  (MLE; hard lockdown 14 Dec, wk14-19)
  mu_EI = 1.37,   # FIXED: Lauer 2020
  mu_IR = 2.41,   # FIXED: He 2020
  eta_y = 0.90,   # FIXED: ~10% immune young (W1 attack rate)
  eta_o = 0.86,   # FIXED: ~14% immune old  (W1 hit old harder)
  rho_y = 0.160,  # 16% detection young (PCR scaling up Oct-Nov 2020)
  rho_o = 0.360,  # 36% detection old   (priority testing + care homes)
  k_y   = 30,     # FIXED: age-specific NegBin dispersion young
  k_o   = 10,     # FIXED: age-specific NegBin dispersion old
  N_y   = N_y, N_o = N_o,
  C_yy  = C_yy, C_yo = C_yo, C_oy = C_oy, C_oo = C_oo
)
coef(covid_w2) <- theta_start

# ── 9. PRIOR ─────────────────────────────────────────────────────────────────
# log(0.340) = -1.0788;  log(0.200) = -1.6094
# log(0.430) = -0.8440;  log(0.220) = -1.5141
# logit(0.160) = -1.6582; logit(0.360) = -0.5754
#
# Beta1 SD=0.20: autumn growth well-constrained by 6 data points
# Beta2 SD=0.25: trough phase — partial measures had variable compliance
# Beta3 SD=0.20: re-acceleration — short 3-week phase, tighten prior
# Beta4 SD=0.25: lockdown effectiveness uncertain (household transmission)
# rho_y SD=0.20: PCR detection in autumn 2020 — PIENTER-informed
# rho_o SD=0.35: old detection heterogeneous across settings

covid_dprior <- Csnippet("
  double lB1    = log(Beta1);
  double lB2    = log(Beta2);
  double lB3    = log(Beta3);
  double lB4    = log(Beta4);
  double lrho_y = log(rho_y / (1.0 - rho_y));
  double lrho_o = log(rho_o / (1.0 - rho_o));

  lik = dnorm(lB1,   -1.0788, 0.20, 1)
      + dnorm(lB2,   -1.6094, 0.25, 1)
      + dnorm(lB3,   -0.8440, 0.20, 1)
      + dnorm(lB4,   -1.5141, 0.25, 1)
      + dnorm(lrho_y,-1.6582, 0.20, 1)
      + dnorm(lrho_o,-0.5754, 0.35, 1);

  if (!give_log) lik = exp(lik);
")
covid_w2 <- pomp(covid_w2, dprior = covid_dprior,
                 paramnames = names(theta_start))

# ── 10. SIMULATION CHECK ──────────────────────────────────────────────────────
# Expected shape: two humps with trough between wk8-11.
# Hump 1 peak ~wk7 (~49k young), trough ~wk10 (~27k young),
# hump 2 peak ~wk14-15 (~60k young), decline to ~28k by wk19.
set.seed(2024)
sim_df <- covid_w2 |>
  simulate(nsim = 50, format = "data.frame", include.data = TRUE) |>
  mutate(reports_total = reports_young + reports_old)

# Phase labels pinned to Young panel only (avoids facet_wrap annotation error)
phase_labels <- data.frame(
  week   = c(3.0, 8.5, 12.0, 16.5),
  label  = c("Growth", "Trough", "Re-accel", "Lockdown"),
  series = factor("Young (<60)",
                  levels = c("Young (<60)", "Old (\u226560)", "Total"))
)

sim_df |>
  pivot_longer(c(reports_young, reports_old, reports_total),
               names_to = "series", values_to = "cases") |>
  mutate(series = recode(series,
                         "reports_young" = "Young (<60)",
                         "reports_old"   = "Old (\u226560)",
                         "reports_total" = "Total"),
         series = factor(series, levels = c("Young (<60)", "Old (\u226560)", "Total"))) |>
  ggplot(aes(x = week, y = cases, group = .id,
             colour    = (.id == "data"),
             linewidth = (.id == "data"),
             alpha     = (.id == "data"))) +
  geom_line() +
  geom_vline(xintercept = c(6.5, 10.5, 14.5), linetype = "dotted",
             colour = "grey50", linewidth = 0.4) +
  geom_text(data = phase_labels,
            aes(x = week, y = Inf, label = label),
            inherit.aes = FALSE,
            vjust = 1.4, size = 2.5, colour = "grey40") +
  facet_wrap(~series, scales = "free_y", nrow = 3) +
  scale_colour_manual(values = c("TRUE"="black","FALSE"="#2166ac"), guide="none") +
  scale_linewidth_manual(values = c("TRUE"=1.0,"FALSE"=0.3), guide="none") +
  scale_alpha_manual(values = c("TRUE"=1.0,"FALSE"=0.4), guide="none") +
  scale_y_continuous(labels = scales::comma) +
  labs(title    = "Wave 2 (age-structured): Simulation check",
       subtitle = "4-phase NPI | Full 19-week window | Both humps included",
       x = "Week", y = "Weekly reported cases") +
  theme_bw(base_size = 10)

# ── 11. PARTICLE FILTER DIAGNOSTIC ───────────────────────────────────────────
Np_use <- 3000
set.seed(999)
ll_check <- replicate(10, logLik(pfilter(covid_w2, Np = Np_use)))
ll_check  <- ll_check[is.finite(ll_check)]
cat(sprintf("\npfilter (Np=%d): mean=%.2f  SD=%.3f\n",
            Np_use, mean(ll_check), sd(ll_check)))
if (sd(ll_check) > 0.5) {
  cat("  CAUTION: SD > 0.5 — increase Np_use\n")
  Np_use <- Np_use * 2
  cat(sprintf("  Retrying with Np=%d\n", Np_use))
  ll2 <- replicate(10, logLik(pfilter(covid_w2, Np = Np_use)))
  cat(sprintf("  New SD=%.3f\n", sd(ll2[is.finite(ll2)])))
} else {
  cat("  OK — proceed to PMCMC\n")
}

# ── 12. PHASE 1 PMCMC ────────────────────────────────────────────────────────
pmcmc_p1 <- pmcmc(
  covid_w2,
  Nmcmc    = 5000,
  Np       = Np_use,
  proposal = mvn_diag_rw(c(Beta1=0.08, Beta2=0.10,
                           Beta3=0.08, Beta4=0.10,
                           rho_y=0.08, rho_o=0.12)^2)
)
acc_p1 <- pmcmc_p1@accepts / 5000
cat(sprintf("Phase 1 acceptance rate: %.3f\n", acc_p1))

# ── 13. PHASE 2 PMCMC ────────────────────────────────────────────────────────
chain_p1  <- as.matrix(as.data.frame(traces(pmcmc_p1)))
post_cov  <- cov(chain_p1[2001:5000, params_est])
opt_scale <- (2.38^2) / length(params_est)   # d=6 → 0.9443

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
    Rt_growth  = Beta1 * dom_eigen,   # autumn growth
    Rt_trough  = Beta2 * dom_eigen,   # partial suppression trough
    Rt_reaccel = Beta3 * dom_eigen,   # Christmas re-acceleration
    Rt_lockdown= Beta4 * dom_eigen    # hard lockdown
  )

cat("\nPosterior means:\n")
print(round(colMeans(post_df[, params_est]), 5))
cat("\n95% CIs:\n")
print(round(apply(post_df[, params_est], 2, quantile, c(0.025, 0.975)), 5))
cat("\nDerived R values:\n")
for (v in c("Rt_growth","Rt_trough","Rt_reaccel","Rt_lockdown"))
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
  labs(title="Wave 2 (age-structured): Posterior marginal distributions",
       x="Parameter value", y="Density") +
  theme_bw(base_size=11)

# ── 17. PAIRWISE SCATTER ──────────────────────────────────────────────────────
tryCatch(
  print(
    ggpairs(post_df |> select(all_of(params_est)),
            lower = list(continuous = wrap("points", alpha=0.04, size=0.4)),
            upper = list(continuous = wrap("cor", size=3)),
            diag  = list(continuous = wrap("densityDiag"))) +
      theme_bw(base_size=9) +
      labs(title="Wave 2 (age-structured): Pairwise posterior scatter")
  ),
  error = function(e) pairs(as.matrix(prod_chain[, params_est]),
                            pch=".", col=rgb(0,0,0,.04))
)

# ── 18. POSTERIOR PREDICTIVE CHECK — ALL THREE PANELS ─────────────────────────
set.seed(42)
pp_sims <- lapply(sample(nrow(post_df), 500), function(i) {
  th            <- coef(covid_w2)
  th[params_est] <- as.numeric(post_df[i, params_est])
  simulate(covid_w2, params=th, nsim=1, format="data.frame") |>
    mutate(draw=i,
           reports_total = reports_young + reports_old)
})

pp_long <- bind_rows(pp_sims) |>
  pivot_longer(c(reports_young, reports_old, reports_total),
               names_to="series", values_to="cases") |>
  mutate(series = recode(series,
                         "reports_young" = "Young (<60)",
                         "reports_old"   = "Old (\u226560)",
                         "reports_total" = "Total"))

pp_band <- pp_long |>
  group_by(week, series) |>
  summarise(lo  = quantile(cases, 0.025, na.rm=TRUE),
            hi  = quantile(cases, 0.975, na.rm=TRUE),
            med = median(cases,   na.rm=TRUE),
            .groups="drop") |>
  mutate(series = factor(series, levels=c("Young (<60)","Old (\u226560)","Total")))

obs_long <- meas_w2 |>
  mutate(reports_total = reports_young + reports_old) |>
  pivot_longer(c(reports_young, reports_old, reports_total),
               names_to="series", values_to="cases") |>
  mutate(series = recode(series,
                         "reports_young" = "Young (<60)",
                         "reports_old"   = "Old (\u226560)",
                         "reports_total" = "Total"),
         series = factor(series, levels=c("Young (<60)","Old (\u226560)","Total")))

series_cols <- c("Young (<60)"="#2166ac",
                 "Old (\u226560)"  ="#d73027",
                 "Total"      ="#1a7a1a")

# Phase bands as background shading
phase_rects <- data.frame(
  xmin  = c(0.5,  6.5, 10.5, 14.5),
  xmax  = c(6.5, 10.5, 14.5, 19.5),
  phase = c("Growth","Trough","Re-accel","Lockdown"),
  fill  = c("#e8f4fd","#fef9e7","#fde8e8","#e8f4fd")
)

ggplot() +
  geom_rect(data = phase_rects,
            aes(xmin=xmin, xmax=xmax, ymin=-Inf, ymax=Inf, fill=fill),
            alpha=0.25, inherit.aes=FALSE) +
  scale_fill_identity() +
  geom_ribbon(data=pp_band,
              aes(x=week, ymin=lo, ymax=hi, fill=series),
              alpha=0.30) +
  geom_line(data=pp_band,
            aes(x=week, y=med, colour=series),
            linewidth=1) +
  geom_point(data=obs_long,
             aes(x=week, y=cases, colour=series),
             size=2, shape=16) +
  geom_vline(xintercept=c(6.5,10.5,14.5), linetype="dotted",
             colour="grey50", linewidth=0.3) +
  facet_wrap(~series, scales="free_y", nrow=3) +
  scale_colour_manual(values=series_cols, guide="none") +
  scale_fill_manual(values=series_cols, guide="none") +
  scale_y_continuous(labels=scales::comma) +
  labs(title    = "Wave 2 (age-structured): Posterior predictive check",
       subtitle = "Shaded=95% interval; line=posterior median; dots=data | 4-phase NPI",
       x="Week", y="Weekly reported cases") +
  theme_bw(base_size=11)

# ── 19. AGGREGATE PPC — TOTAL CASES ONLY ─────────────────────────────────────
pp_total <- bind_rows(pp_sims) |>
  mutate(reports_total = reports_young + reports_old) |>
  group_by(week) |>
  summarise(
    lo  = quantile(reports_total, 0.025, na.rm = TRUE),
    hi  = quantile(reports_total, 0.975, na.rm = TRUE),
    med = median(reports_total,   na.rm = TRUE),
    .groups = "drop"
  )

obs_total <- meas_w2 |>
  mutate(reports_total = reports_young + reports_old)

ggplot() +
  geom_rect(data = phase_rects,
            aes(xmin=xmin, xmax=xmax, ymin=-Inf, ymax=Inf, fill=fill),
            alpha=0.25, inherit.aes=FALSE) +
  scale_fill_identity() +
  geom_ribbon(data = pp_total,
              aes(x = week, ymin = lo, ymax = hi),
              fill = "#2166ac", alpha = 0.25) +
  geom_line(data = pp_total,
            aes(x = week, y = med),
            colour = "#2166ac", linewidth = 1.2) +
  geom_point(data = obs_total,
             aes(x = week, y = reports_total),
             colour = "black", size = 2.5, shape = 16) +
  geom_vline(xintercept=c(6.5,10.5,14.5), linetype="dotted",
             colour="grey50", linewidth=0.3) +
  scale_y_continuous(labels = scales::comma) +
  labs(
    title    = "Wave 2 (age-structured): Posterior predictive check — total cases",
    subtitle = "Total = Young (<60) + Old (\u226560)  |  Shaded=95% interval; line=median; dots=data",
    x        = "Week",
    y        = "Weekly reported cases (total)"
  ) +
  theme_bw(base_size = 12)

# Make sure results/wave2 folder exists
dir.create("results/wave2", showWarnings = FALSE)
dir.create("results/diagnostics/baseline_wave2", showWarnings = FALSE)

# The plots were already printed to your screen during the wave script run.
# To regenerate and save them, just re-run the plot sections of the wave script.
# Or save the current session objects:
saveRDS(prod_chain, "results/baseline_wave2.rds")  

setwd("C:/Users/Usuario/Desktop/CBS - Internship/Analysis/DATA")
source("run_wave2_sensitivity.R")
source("sensitivity_plots_wave2.R")   # fills results/wave2/
# Trace plots — fast, runs in seconds
source("sensitivity_traceplots_wave2.R")

# Posteriors, correlations, PPC — slow, ~30-60 min
source("sensitivity_diagnostics_wave2.R")
