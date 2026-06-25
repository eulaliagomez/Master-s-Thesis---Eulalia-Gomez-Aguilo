# =============================================================================
# covid19_age_wave4.R
#
# Wave 4: 21 Jun – 3 Oct 2021  |  Delta variant (B.1.617.2)
# Two age groups: young (<60) and old (>=60)
#
# Key parameter decisions:
#
#   Beta1 = 1.426  → R0 = 5.81  (MLE)
#     In the contact-matrix model, R0 = Beta * dominant_eigenvalue(C/mu_IR)
#     = Beta * 4.074.  Delta R0 ≈ 5.8 (Liu & Rocklöv 2021, J Travel Med).
#     Note: the single-population Wave 4 script used Beta1=14.0/wk, which
#     equals R0=5.81 in that simpler model (R0 = Beta/mu_IR = 14/2.41).
#     The correct age-structured equivalent is Beta1 = R0/dom_eigenvalue =
#     5.81/4.074 = 1.426.  Using Beta1=14 in the contact-matrix model
#     produces an effective transmission amplified by C_yy=9.49, which
#     depletes S_y to zero within two weeks — destroying the plateau.
#
#   t_lock = 2  (NPI transition at week 2→3)
#     Emergency NPI package announced 9 July 2021. Rapid behaviour change
#     and nightclub closures hit just as the epidemic peaked at wk3.
#
#   t_open = 7  (NPI transition at week 7→8)
#     Summer holiday period ends late August; schools reopen 6 Sep 2021.
#     Phase 2 (B3) captures the sustained plateau from the school-term.
#
#   Beta2 = 0.240  → Rt = 0.98  (MLE; slight Rt<1 during summer holidays)
#     Post-emergency-NPI + voluntary behaviour change → rapid decline.
#
#   Beta3 = 0.450  → Rt = 1.83  (MLE; school-term plateau)
#     Beta3 > Beta2: plateau phase has HIGHER transmission than the decline
#     phase, driven by school-reopening contact rates.  This is the same
#     structural feature as the original single-population Wave 4 script
#     (B3=4.0 > B2=3.0 there).  With S_y ≈ 6.2M at wk8 (not depleted,
#     because Beta1 is correctly scaled), Beta3=0.450 sustains ~14k/wk young.
#
#   eta_y = 0.64  (FIXED)
#     36% of young effectively immune at wave start:
#     W1+W2+W3 cumulative attack rate ≈ 25–30% (PIENTER-Corona rd 3);
#     partial vaccination coverage (18–35 yr eligible from mid-June 2021,
#     ~15–25% covered by 21 June start).  Combined ≈ 35–40%; eta=0.64 is
#     the best-fitting value within this range (LL improves by 47 units
#     vs eta=0.55 from the original script).
#     Source: RIVM vaccinatiecijfers Jun 2021; PIENTER-Corona rd 3.
#
#   eta_o = 0.25  (FIXED)
#     ~75% of old effectively immune: 80%+ 2-dose coverage for 60+
#     × VE≈90% → ~72–77%; plus prior-wave attack rate ~5–8%.
#     Conservative eta=0.25 accounts for partial waning.
#     Source: RIVM vaccinatiecijfers; Backer et al. 2021 Euro Surveill.
#
#   rho_y = 0.140  (MLE)
#     14% detection rate young.  Lower than Wave 3 (24%) because Delta
#     wave was primarily mild community infections; young people were less
#     likely to present for testing.  Consistent with PIENTER-Corona
#     seroprevalence estimates for summer 2021.
#
#   rho_o = 0.320  (MLE)
#     32% detection rate old.  Lower than Wave 3 (40%) because vaccinated
#     old had milder symptoms and lower presentation to testing services.
#
#   k_y = 30, k_o = 10  (FIXED, age-specific)
#     Same justification as Wave 3: young group had stable broad community
#     PCR testing (k_y=30 tighter); old group had more heterogeneous
#     ascertainment (k_o=10 more dispersed).
#
# NPI phases (3-phase, matching Wave 1 structure):
#   Phase 0 (t < 2):  fully open — restrictions lifted 26 Jun, Delta
#                     spreading freely → Beta1 (highest)
#   Phase 1 (2≤t<7):  emergency NPI (9 Jul) + summer holidays → Beta2
#                     (lowest; rapid post-peak decline)
#   Phase 2 (t≥7):    schools reopen (6 Sep) → Beta3
#                     (intermediate; sustained plateau, Beta3 > Beta2)
#
# Contact matrix: Prem et al. 2021 NLD, all settings, 16x16 → 2x2 balanced
#   C_yy=9.4901, C_yo=0.7523, C_oy=2.8524, C_oo=3.2609
#   Source: https://doi.org/10.1371/journal.pcbi.1009098
#
# MLE grid search (Python): LL = -268.73
#   wk3 peak young: 63,971 predicted vs 66,970 observed (-4.5%)
#   wk4 young:      54,262 predicted vs 54,632 observed (-0.7%)
#   plateau wk8-12: 13,500-14,100 predicted vs 15,000-15,400 observed (-9%)
# =============================================================================

library(tidyverse)
library(pomp)
library(coda)
library(GGally)

set.seed(123456)

# ── 1. DATA ──────────────────────────────────────────────────────────────────
meas_y <- read_csv("covid19_wave4_young_NL.csv") |>
  rename(reports_young = reports) |> select(week, reports_young)
meas_o <- read_csv("covid19_wave4_old_NL.csv") |>
  rename(reports_old = reports) |> select(week, reports_old)
meas_w4 <- left_join(meas_y, meas_o, by = "week")

cat("Wave 4 data:\n")
cat(sprintf("  T = %d weeks\n", nrow(meas_w4)))
cat(sprintf("  Young: peak=%s (wk%d), total=%s\n",
            format(max(meas_w4$reports_young), big.mark=","),
            which.max(meas_w4$reports_young),
            format(sum(meas_w4$reports_young), big.mark=",")))
cat(sprintf("  Old:   peak=%s (wk%d), total=%s\n",
            format(max(meas_w4$reports_old), big.mark=","),
            which.max(meas_w4$reports_old),
            format(sum(meas_w4$reports_old), big.mark=",")))

# ── 2. CONTACT MATRIX & POPULATION ───────────────────────────────────────────
C_yy <- 9.4901; C_yo <- 0.7523; C_oy <- 2.8524; C_oo <- 3.2609
N_y  <- 12736000; N_o <- 3359000

# ── 3. NPI COVARIATE ─────────────────────────────────────────────────────────
# Phase 0 (t < 2):  fully open — Delta spreading freely — Beta1
# Phase 1 (2≤t<7):  emergency NPI (9 Jul) + summer holidays — Beta2 (lowest)
# Phase 2 (t≥7):   school reopening (6 Sep) — Beta3 (Beta3 > Beta2 = plateau)
#
# Structural note: Beta3 > Beta2 is the same feature as the non-age-structured
# Wave 4 script (B3=4 > B2=3). The plateau arises because:
#   (a) S_y is not depleted (Beta1 correctly scaled → epidemic burns <20% of S)
#   (b) Schools reopening raises effective contact rates above summer holiday level

npi_df    <- data.frame(week      = 0:18,
                        npi_phase = c(rep(0, 2), rep(1, 5), rep(2, 12)))
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
# Age-specific NegBin dispersion:
#   k_y = 30: young group — stable community PCR, tighter observation model
#   k_o = 10: old group — heterogeneous ascertainment, more overdispersion
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
# I0_y = 4149 / (0.140 * 2.41) = 12,297;  E0_y = 12,297 * 2.41/1.37 = 21,632
# I0_o =  228 / (0.320 * 2.41) =    296;  E0_o =    296 * 2.41/1.37 =    520

seir_rinit <- Csnippet("
  double I0_y = 12297.0;
  double E0_y = 21632.0;
  double I0_o =   296.0;
  double E0_o =   520.0;
  S_y = nearbyint(eta_y * N_y) - nearbyint(I0_y) - nearbyint(E0_y);
  E_y = nearbyint(E0_y); I_y = nearbyint(I0_y);
  R_y = nearbyint((1.0 - eta_y) * N_y); H_y = 0.0;
  S_o = nearbyint(eta_o * N_o) - nearbyint(I0_o) - nearbyint(E0_o);
  E_o = nearbyint(E0_o); I_o = nearbyint(I0_o);
  R_o = nearbyint((1.0 - eta_o) * N_o); H_o = 0.0;
")

# ── 7. BUILD pomp OBJECT ─────────────────────────────────────────────────────
params_est <- c("Beta1", "Beta2", "Beta3", "rho_y", "rho_o")

covid_w4 <- meas_w4 |>
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
                   "eta_y","eta_o","rho_y","rho_o","k_y","k_o",
                   "N_y","N_o","C_yy","C_yo","C_oy","C_oo"),
    statenames = c("S_y","E_y","I_y","R_y","H_y",
                   "S_o","E_o","I_o","R_o","H_o"),
    accumvars  = c("H_y","H_o"),
    obsnames   = c("reports_young","reports_old")
  )

# ── 8. STARTING PARAMETERS ────────────────────────────────────────────────────
theta_start <- c(
  Beta1 = 1.426,  # R0 = 5.81  (MLE; Delta R0 = Beta * dom_eigenvalue 4.074)
  Beta2 = 0.240,  # Rt = 0.98  (MLE; summer holiday + NPI; slight Rt<1)
  Beta3 = 0.450,  # Rt = 1.83  (MLE; school reopening plateau; Beta3>Beta2)
  mu_EI = 1.37,   # FIXED: Lauer 2020
  mu_IR = 2.41,   # FIXED: He 2020
  eta_y = 0.64,   # FIXED: 36% pre-immune young (W1+W2+W3 attack + partial vax)
  eta_o = 0.25,   # FIXED: 75% pre-immune old  (vax + prior waves)
  rho_y = 0.140,  # 14% detection young (mild Delta, less testing uptake)
  rho_o = 0.320,  # 32% detection old  (vaccinated, milder symptoms)
  k_y   = 30,     # FIXED: age-specific NegBin dispersion young
  k_o   = 10,     # FIXED: age-specific NegBin dispersion old
  N_y   = N_y, N_o = N_o,
  C_yy  = C_yy, C_yo = C_yo, C_oy = C_oy, C_oo = C_oo
)
coef(covid_w4) <- theta_start

# ── 9. PRIOR ─────────────────────────────────────────────────────────────────
# log(1.426) =  0.355;  log(0.240) = -1.427;  log(0.450) = -0.799
# logit(0.140) = -1.815;  logit(0.320) = -0.754
#
# Beta1 SD=0.35: Delta R0 range 5–7 in literature; prior covers this well
# Beta2 SD=0.40: post-NPI decline phase; wider since compliance varied
# Beta3 SD=0.40: school-term plateau; wider since September contact rates uncertain
# rho_y SD=0.20: stable community PCR; tight (same justification as Wave 3)
# rho_o SD=0.35: old detection less certain (vaccine effect on behaviour)

covid_dprior <- Csnippet("
  double lB1    = log(Beta1);
  double lB2    = log(Beta2);
  double lB3    = log(Beta3);
  double lrho_y = log(rho_y / (1.0 - rho_y));
  double lrho_o = log(rho_o / (1.0 - rho_o));

  lik = dnorm(lB1,    0.355, 0.35, 1)
      + dnorm(lB2,   -1.427, 0.40, 1)
      + dnorm(lB3,   -0.799, 0.40, 1)
      + dnorm(lrho_y,-1.815, 0.20, 1)
      + dnorm(lrho_o,-0.754, 0.35, 1);

  if (!give_log) lik = exp(lik);
")
covid_w4 <- pomp(covid_w4, dprior = covid_dprior,
                 paramnames = names(theta_start))

# ── 10. SIMULATION CHECK ──────────────────────────────────────────────────────
# Simulations should show:
#   - Explosive rise to peak wk3 (~60-70k young)
#   - Rapid decline wk4-7
#   - Sustained plateau wk8-15 (~13-15k total), NOT near-zero
# The plateau is produced naturally (not by a 4th phase) because:
#   (a) Beta1=1.426 preserves S_y — only ~15% of pool consumed at peak
#   (b) Beta3=0.450 > Beta2=0.240 — school term raises transmission above
#       summer holiday level

set.seed(2024)
sim_df <- covid_w4 |>
  simulate(nsim = 50, format = "data.frame", include.data = TRUE) |>
  mutate(reports_total = reports_young + reports_old)

# Phase label data: one row per label in Young panel only
phase_labels <- data.frame(
  week   = c(1.0, 4.5, 11.0),
  label  = c("Free", "Summer NPI", "Schools"),
  series = factor("Young (<60)", levels = c("Young (<60)", "Old (\u226560)", "Total"))
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
  geom_vline(xintercept = c(2.5, 7.5), linetype = "dotted",
             colour = "grey50", linewidth = 0.4) +
  geom_text(data = phase_labels,
            aes(x = week, y = Inf, label = label),
            inherit.aes = FALSE,
            vjust = 1.4, size = 2.8, colour = "grey40") +
  facet_wrap(~series, scales = "free_y", nrow = 3) +
  scale_colour_manual(values = c("TRUE"="black","FALSE"="#2166ac"), guide="none") +
  scale_linewidth_manual(values = c("TRUE"=1.0,"FALSE"=0.3), guide="none") +
  scale_alpha_manual(values = c("TRUE"=1.0,"FALSE"=0.4), guide="none") +
  scale_y_continuous(labels = scales::comma) +
  labs(title    = "Wave 4 (age-structured): Simulation check",
       subtitle = paste("Check: peak wk3, rapid decline wk4-7,",
                        "sustained plateau wk8-15 (not near-zero)"),
       x = "Week", y = "Weekly reported cases") +
  theme_bw(base_size = 10)

# ── 11. PARTICLE FILTER DIAGNOSTIC ───────────────────────────────────────────
Np_use <- 3000
set.seed(999)
ll_check <- replicate(10, logLik(pfilter(covid_w4, Np = Np_use)))
ll_check  <- ll_check[is.finite(ll_check)]
cat(sprintf("\npfilter (Np=%d): mean=%.2f  SD=%.3f\n",
            Np_use, mean(ll_check), sd(ll_check)))
if (sd(ll_check) > 0.5) {
  cat("  CAUTION: SD > 0.5 — increase Np_use\n")
  Np_use <- Np_use * 2
  cat(sprintf("  Retrying with Np=%d\n", Np_use))
  ll2 <- replicate(10, logLik(pfilter(covid_w4, Np = Np_use)))
  cat(sprintf("  New SD=%.3f\n", sd(ll2[is.finite(ll2)])))
} else {
  cat("  OK — proceed to PMCMC\n")
}

# ── 12. PHASE 1 PMCMC ────────────────────────────────────────────────────────
pmcmc_p1 <- pmcmc(
  covid_w4,
  Nmcmc    = 5000,
  Np       = Np_use,
  proposal = mvn_diag_rw(c(Beta1=0.15, Beta2=0.15, Beta3=0.15,
                           rho_y=0.10, rho_o=0.15)^2)
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
    R0_delta   = Beta1 * dom_eigen,   # Delta R0 — free spread
    Rt_summer  = Beta2 * dom_eigen,   # Rt summer holiday + NPI
    Rt_schools = Beta3 * dom_eigen    # Rt school-term plateau
  )

cat("\nPosterior means:\n")
print(round(colMeans(post_df[, params_est]), 5))
cat("\n95% CIs:\n")
print(round(apply(post_df[, params_est], 2, quantile, c(0.025, 0.975)), 5))
cat("\nDerived R values:\n")
for (v in c("R0_delta","Rt_summer","Rt_schools"))
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
  labs(title="Wave 4 (age-structured): Posterior marginal distributions",
       x="Parameter value", y="Density") +
  theme_bw(base_size=11)

# ── 17. PAIRWISE SCATTER ──────────────────────────────────────────────────────
ggpairs(post_df |> select(all_of(params_est)),
        lower = list(continuous = wrap("points", alpha=0.04, size=0.4)),
        upper = list(continuous = wrap("cor", size=3)),
        diag  = list(continuous = wrap("densityDiag"))) +
  theme_bw(base_size=9) +
  labs(title="Wave 4 (age-structured): Pairwise posterior scatter")

# ── 18. POSTERIOR PREDICTIVE CHECK — ALL THREE PANELS ─────────────────────────
set.seed(42)
pp_sims <- lapply(sample(nrow(post_df), 500), function(i) {
  th            <- coef(covid_w4)
  th[params_est] <- as.numeric(post_df[i, params_est])
  simulate(covid_w4, params=th, nsim=1, format="data.frame") |>
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

obs_long <- meas_w4 |>
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
  labs(title    = "Wave 4 (age-structured): Posterior predictive check",
       subtitle = "Shaded=95% interval; line=posterior median; dots=data",
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

obs_total <- meas_w4 |>
  mutate(reports_total = reports_young + reports_old)

ggplot() +
  geom_ribbon(data = pp_total,
              aes(x = week, ymin = lo, ymax = hi),
              fill = "#2166ac", alpha = 0.25) +
  geom_line(data = pp_total,
            aes(x = week, y = med),
            colour = "#2166ac", linewidth = 1.2) +
  geom_point(data = obs_total,
             aes(x = week, y = reports_total),
             colour = "black", size = 2.5, shape = 16) +
  scale_y_continuous(labels = scales::comma) +
  labs(
    title    = "Wave 4 (age-structured): Posterior predictive check — total cases",
    subtitle = "Total = Young (<60) + Old (\u226560)  |  Shaded=95% interval; line=median; dots=data",
    x        = "Week",
    y        = "Weekly reported cases (total)"
  ) +
  theme_bw(base_size = 12)


saveRDS(prod_chain, "results/baseline_wave4.rds")

# Run sensitivity analysis
source("run_wave4_sensitivity.R")

# Figures
source("sensitivity_plots_wave4.R")

# Trace plots — fast
source("sensitivity_traceplots_wave4.R")

# Diagnostics — slow
source("sensitivity_diagnostics_wave4.R")