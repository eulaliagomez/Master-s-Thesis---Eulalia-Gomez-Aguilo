# =============================================================================
# covid19_nl_wave3.R
#
# Wave 3: 1 February 2021 – 20 June 2021
# Dominant variant: Alpha (B.1.1.7)
# Context: wave begins UNDER the December 2020 lockdown; schools reopen
#          1 March (week 4-5) and Alpha spreads freely; third lockdown
#          6 April (week 10) stops growth; gradual reopening May-June 2021.
#
# Differences from Wave 1:
#   - NPI PHASES ARE INVERTED: starts under lockdown (low Beta), then opens
#     (high Beta), then re-locks (medium Beta). Phase labels map differently.
#   - Alpha variant: ~50-70% more transmissible than ancestral
#   - eta = 0.68: ~32% immune from Waves 1+2 combined plus early vaccination
#   - rho = 0.20: further expanded testing
#   - Wave starts at 24,690 cases/week (high baseline from W2 end)
# =============================================================================

library(tidyverse)
library(pomp)
library(coda)
library(GGally)

set.seed(123456)

# ─────────────────────────────────────────────────────────────────────────────
# 1.  DATA
# ─────────────────────────────────────────────────────────────────────────────

daily_raw <- read_csv("covid19_daily_cases_NL.csv") |>
  mutate(date = as.Date(date))

weekly <- daily_raw |>
  mutate(year_week = floor_date(date, unit = "week", week_start = 1)) |>
  group_by(year_week) |>
  summarise(reports = sum(cases), .groups = "drop") |>
  arrange(year_week)

meas_w3 <- weekly |>
  filter(year_week >= as.Date("2021-02-01"),
         year_week <= as.Date("2021-06-20")) |>
  mutate(week = seq_len(n())) |>
  select(week, reports)

cat("Wave 3:", nrow(meas_w3), "weekly obs | peak:", max(meas_w3$reports),
    "(wk", which.max(meas_w3$reports), ") | total:", sum(meas_w3$reports), "\n")

# ─────────────────────────────────────────────────────────────────────────────
# 2.  NPI COVARIATE
# ─────────────────────────────────────────────────────────────────────────────
# Wave 3 NPI timeline (relative to wave start, 1 Feb 2021):
#   Phase 0 (t < 4):  December lockdown still in force → Beta1 (LOW)
#                     Slow growth under suppression
#   Phase 1 (4 ≤ t < 10): schools reopen 1 March → Alpha spreads → Beta2 (HIGH)
#                          Epidemic accelerates through March/early April
#   Phase 2 (t ≥ 10): third lockdown 6 April → Beta3 (MEDIUM, declining)
#                     Gradual reopening May-June sustains slow decline
#
# IMPORTANT: unlike Wave 1, Beta2 > Beta1 here. Phase 0 is the suppressed phase.

npi_df    <- data.frame(week      = 0:23,
                        npi_phase = c(rep(0, 4), rep(1, 6), rep(2, 14)))
npi_covar <- covariate_table(npi_df, times = "week")

# ─────────────────────────────────────────────────────────────────────────────
# 3.  PROCESS MODEL
# ─────────────────────────────────────────────────────────────────────────────

seir_step_npi <- Csnippet("
  double eff_beta;
  if      (npi_phase < 0.5) eff_beta = Beta1;
  else if (npi_phase < 1.5) eff_beta = Beta2;
  else                       eff_beta = Beta3;

  double dN_SE = rbinom(S, 1 - exp(-eff_beta * I / N * dt));
  double dN_EI = rbinom(E, 1 - exp(-mu_EI * dt));
  double dN_IR = rbinom(I, 1 - exp(-mu_IR * dt));

  S -= dN_SE;
  E += dN_SE - dN_EI;
  I += dN_EI - dN_IR;
  R += dN_IR;
  H += dN_IR;
")

# ─────────────────────────────────────────────────────────────────────────────
# 4.  OBSERVATION MODEL
# ─────────────────────────────────────────────────────────────────────────────

covid_dmeas <- Csnippet("lik = dnbinom_mu(reports, k, rho * H, give_log);")
covid_rmeas <- Csnippet("reports = rnbinom_mu(k, rho * H);")

# ─────────────────────────────────────────────────────────────────────────────
# 5.  INITIAL CONDITIONS
# ─────────────────────────────────────────────────────────────────────────────
# Wave 3 starts at 24,690 reported cases in week 1.
# Data-implied I0 = 24690 / (0.20 × 2.41) ≈ 51224
# Data-implied E0 = 51224 × (2.41/1.37) ≈ 90109

seir_rinit <- Csnippet("
  double I_init = 51224.0;
  double E_init = 90109.0;
  I = nearbyint(I_init);
  E = nearbyint(E_init);
  S = nearbyint(eta * N) - I - E;
  R = nearbyint((1.0 - eta) * N);
  H = 0.0;
")

# ─────────────────────────────────────────────────────────────────────────────
# 6.  BUILD pomp OBJECT
# ─────────────────────────────────────────────────────────────────────────────
# Estimated: Beta1, Beta2, Beta3, rho, k
# Fixed:     N, mu_EI, mu_IR, eta
#   — eta = 0.68: ~32% immune from W1+W2 seroprevalence + early vaccination

params_est <- c("Beta1", "Beta2", "Beta3", "rho", "k")

covid_w3 <- meas_w3 |>
  pomp(
    times      = "week",
    t0         = 0,
    rprocess   = euler(seir_step_npi, delta.t = 1/7),
    rinit      = seir_rinit,
    rmeasure   = covid_rmeas,
    dmeasure   = covid_dmeas,
    covar      = npi_covar,
    partrans   = parameter_trans(
      log   = c("Beta1", "Beta2", "Beta3", "k"),
      logit = c("rho")
    ),
    paramnames = c("N","Beta1","Beta2","Beta3","mu_EI","mu_IR","eta","k","rho"),
    statenames = c("S","E","I","R","H"),
    accumvars  = "H"
  )

# ─────────────────────────────────────────────────────────────────────────────
# 7.  STARTING PARAMETERS
# ─────────────────────────────────────────────────────────────────────────────
# Beta1 = 3.5/wk → Rt=1.45 under Dec lockdown; Alpha baseline under suppression
# Beta2 = 4.0/wk → R0_Alpha=1.66; Alpha spreading freely after school reopening
#                  Alpha is ~50-70% more transmissible than ancestral
# Beta3 = 3.0/wk → Rt=1.24 post-re-lockdown; gradual decline through May-Jun
# rho   = 0.20   → expanded testing; ~20% of infections reported
# eta   = 0.95   → FIXED at W1+W2 immunity level (~32% immune → eta=0.68)

theta_start <- c(
  N     = 17400000,
  Beta1 = 3.5,
  Beta2 = 4.0,
  Beta3 = 3.0,
  mu_EI = 1.37,
  mu_IR = 2.41,
  eta   = 0.68,
  rho   = 0.20,
  k     = 10
)

coef(covid_w3) <- theta_start

# ─────────────────────────────────────────────────────────────────────────────
# 8.  SIMULATION CHECK
# ─────────────────────────────────────────────────────────────────────────────
# Wave 3 simulations start high, dip slightly (weeks 1-4 under lockdown),
# then rise to peak around week 10-11, then decline.
# The non-monotonic shape is correct and expected.

set.seed(2024)
covid_w3 |>
  simulate(nsim = 30, format = "data.frame", include.data = TRUE) |>
  ggplot(aes(x = week, y = reports, group = .id,
             colour    = (.id == "data"),
             linewidth = (.id == "data"),
             alpha     = (.id == "data"))) +
  geom_line() +
  scale_colour_manual(values = c("TRUE"="black","FALSE"="#2166ac"), guide="none") +
  scale_linewidth_manual(values = c("TRUE"=1.0,"FALSE"=0.3), guide="none") +
  scale_alpha_manual(values = c("TRUE"=1.0,"FALSE"=0.4), guide="none") +
  scale_y_continuous(labels = scales::comma) +
  labs(title    = "Wave 3: Simulation check",
       subtitle = "Slow start (lockdown) → rise (schools open + Alpha) → peak wk 11 → decline",
       x = "Week", y = "Weekly reported cases") +
  theme_bw(base_size = 11)

# ─────────────────────────────────────────────────────────────────────────────
# 9.  PRIOR SPECIFICATION
# ─────────────────────────────────────────────────────────────────────────────
# Beta1: log-N(1.25, 0.40²)   → centred on log(3.5); lockdown-suppressed range
# Beta2: log-N(1.39, 0.50²)   → centred on log(4.0); Alpha free transmission
# Beta3: log-N(1.10, 0.50²)   → centred on log(3.0); post-re-lockdown decline
# rho:   logit-N(-1.39, 0.60²) → centred on logit(0.20), 95% mass on (0.08,0.42)
# k:     log-N(2.30, 0.80²)

covid_dprior <- Csnippet("
  double lB1  = log(Beta1);
  double lB2  = log(Beta2);
  double lB3  = log(Beta3);
  double lrho = log(rho / (1.0 - rho));
  double lk   = log(k);

  lik = dnorm(lB1,  1.25, 0.40, 1)
      + dnorm(lB2,  1.39, 0.50, 1)
      + dnorm(lB3,  1.10, 0.50, 1)
      + dnorm(lrho,-1.39, 0.60, 1)
      + dnorm(lk,   2.30, 0.80, 1);

  if (!give_log) lik = exp(lik);
")

covid_w3 <- pomp(covid_w3, dprior = covid_dprior,
                 paramnames = c("N","Beta1","Beta2","Beta3",
                                "mu_EI","mu_IR","eta","k","rho"))

# ─────────────────────────────────────────────────────────────────────────────
# 10.  PARTICLE FILTER DIAGNOSTIC
# ─────────────────────────────────────────────────────────────────────────────

Np_use <- 3000

set.seed(999)
ll_check <- replicate(10, logLik(pfilter(covid_w3, Np = Np_use)))
ll_check <- ll_check[is.finite(ll_check)]
cat(sprintf("\npfilter (Np=%d): mean=%.2f  SD=%.3f\n",
            Np_use, mean(ll_check), sd(ll_check)))
if (sd(ll_check) > 0.5)
  cat("  CAUTION: SD > 0.5 — increase Np_use before running PMCMC\n") else
  cat("  OK — proceed to PMCMC\n")

# ─────────────────────────────────────────────────────────────────────────────
# 11.  PHASE 1 PMCMC — DIAGONAL PROPOSAL
# ─────────────────────────────────────────────────────────────────────────────

pmcmc_p1 <- pmcmc(
  covid_w3,
  Nmcmc    = 5000,
  Np       = Np_use,
  proposal = mvn_diag_rw(c(Beta1=0.15, Beta2=0.15, Beta3=0.15,
                            rho=0.15, k=0.20)^2)
)
acc_p1 <- pmcmc_p1@accepts / 5000
cat(sprintf("Phase 1 acceptance rate: %.3f\n", acc_p1))

# llançat 15:07 - finalitzat 15:25

# ─────────────────────────────────────────────────────────────────────────────
# 12.  PHASE 2 PMCMC — ADAPTIVE MVN PROPOSAL
# ─────────────────────────────────────────────────────────────────────────────

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
cat(sprintf("Phase 2 acceptance rate: %.3f  (target 0.20–0.40)\n", acc_p2))

# llançat 17:52 - finalitzat 
# ─────────────────────────────────────────────────────────────────────────────
# 13.  PRODUCTION CHAIN AND DIAGNOSTICS
# ─────────────────────────────────────────────────────────────────────────────

chain_p2   <- as.matrix(as.data.frame(traces(pmcmc_p2)))
prod_chain <- as.mcmc(chain_p2[2001:nrow(chain_p2), params_est])

cat("ESS:\n"); print(round(effectiveSize(prod_chain)))
plot(prod_chain, ask = FALSE)

# ─────────────────────────────────────────────────────────────────────────────
# 14.  POSTERIOR SUMMARIES
# ─────────────────────────────────────────────────────────────────────────────

post_df <- as.data.frame(prod_chain) |>
  mutate(Rt_lock1 = Beta1 / theta_start["mu_IR"],
         R0_alpha = Beta2 / theta_start["mu_IR"],
         Rt_lock2 = Beta3 / theta_start["mu_IR"])

cat("\nPosterior means:\n"); print(round(colMeans(post_df[, params_est]), 4))
cat("\n95% CIs:\n"); print(round(apply(post_df[, params_est], 2, quantile, c(0.025, 0.975)), 4))
for (v in c("Rt_lock1","R0_alpha","Rt_lock2"))
  cat(sprintf("  %s: %.2f [%.2f, %.2f]\n", v,
              mean(post_df[[v]]),
              quantile(post_df[[v]], 0.025),
              quantile(post_df[[v]], 0.975)))

# ─────────────────────────────────────────────────────────────────────────────
# 15.  POSTERIOR DENSITY PLOTS
# ─────────────────────────────────────────────────────────────────────────────

post_df |>
  select(all_of(params_est)) |>
  pivot_longer(everything(), names_to = "parameter", values_to = "value") |>
  ggplot(aes(x = value)) +
  geom_histogram(aes(y = after_stat(density)), bins=50, fill="#2166ac", alpha=0.55) +
  geom_density(colour = "black", linewidth = 0.8) +
  facet_wrap(~parameter, scales = "free", nrow = 2) +
  labs(title = "Wave 3: Posterior marginal distributions",
       x = "Parameter value", y = "Density") +
  theme_bw(base_size = 11)

# ─────────────────────────────────────────────────────────────────────────────
# 16.  PAIRWISE SCATTER
# ─────────────────────────────────────────────────────────────────────────────

ggpairs(post_df |> select(all_of(params_est)),
        lower = list(continuous = wrap("points", alpha=0.04, size=0.4)),
        upper = list(continuous = wrap("cor", size=3)),
        diag  = list(continuous = wrap("densityDiag"))) +
  theme_bw(base_size = 9) +
  labs(title = "Wave 3: Pairwise posterior scatter")

# ─────────────────────────────────────────────────────────────────────────────
# 17.  POSTERIOR PREDICTIVE CHECK
# ─────────────────────────────────────────────────────────────────────────────

set.seed(42)
pp_sims <- lapply(sample(nrow(post_df), 300), function(i) {
  th            <- coef(covid_w3)
  th[params_est] <- as.numeric(post_df[i, params_est])
  simulate(covid_w3, params=th, nsim=1, format="data.frame") |>
    mutate(draw = i)
})

pp_band <- bind_rows(pp_sims) |>
  group_by(week) |>
  summarise(lo=quantile(reports,0.025,na.rm=TRUE),
            hi=quantile(reports,0.975,na.rm=TRUE),
            med=median(reports,na.rm=TRUE), .groups="drop")

ggplot() +
  geom_ribbon(data=pp_band, aes(x=week,ymin=lo,ymax=hi),
              fill="#2166ac", alpha=0.25) +
  geom_line(data=pp_band, aes(x=week,y=med),
            colour="#2166ac", linewidth=1) +
  geom_point(data=as.data.frame(covid_w3)|>select(week,reports),
             aes(x=week,y=reports), colour="black", size=2) +
  scale_y_continuous(labels=scales::comma) +
  labs(title    = "Wave 3: Posterior predictive check",
       subtitle = "Shaded=95% interval; line=median; dots=data",
       x="Week", y="Weekly reported cases") +
  theme_bw(base_size=11)
