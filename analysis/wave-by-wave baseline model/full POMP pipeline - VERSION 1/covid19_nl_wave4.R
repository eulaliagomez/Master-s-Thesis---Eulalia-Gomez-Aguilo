# =============================================================================
# covid19_nl_wave4.R
#
# Wave 4: 21 June 2021 – 3 October 2021
# Dominant variant: Delta (B.1.617.2)
# Context: all restrictions lifted 26 June → explosive rise to peak in 2 weeks;
#          Delta R0 ≈ 5.8 (≈2x Alpha); ~50-55% vaccinated by July 2021;
#          wave does not resolve — settles into sustained plateau ~11-17k/wk.
#
# Differences from Waves 1-3:
#   - FASTEST dynamics: peak at week 3 (T=15 is the shortest wave)
#   - Very high Beta1 (Delta R0 ≈ 5.8) but short free-transmission window
#   - Tail is a PLATEAU (Rt > 1 persists) not a decline to near-zero
#     Beta3 > Beta2 in starting values; model holds cases up in the tail
#   - eta = 0.55: ~45% effectively immune (vaccination + prior infection)
#   - rho = 0.30: highest testing coverage of any wave
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

meas_w4 <- weekly |>
  filter(year_week >= as.Date("2021-06-21"),
         year_week <= as.Date("2021-10-03")) |>
  mutate(week = seq_len(n())) |>
  select(week, reports)

cat("Wave 4:", nrow(meas_w4), "weekly obs | peak:", max(meas_w4$reports),
    "(wk", which.max(meas_w4$reports), ") | total:", sum(meas_w4$reports), "\n")

# ─────────────────────────────────────────────────────────────────────────────
# 2.  NPI COVARIATE
# ─────────────────────────────────────────────────────────────────────────────
# Wave 4 NPI timeline (relative to wave start, 21 June 2021):
#   Phase 0 (t < 2):  fully open — all restrictions lifted 26 June → Beta1 (HIGH)
#                     Delta's R0 ≈ 5.8 → explosive growth
#   Phase 1 (2 ≤ t < 9): informal behaviour change + partial QR measures → Beta2
#                          Rapid post-peak decline through July-August
#   Phase 2 (t ≥ 9):  QR code system 25 September + sustained plateau → Beta3
#                      Cases stabilise; Rt just above 1 → slow ongoing spread

npi_df    <- data.frame(week      = 0:18,
                        npi_phase = c(rep(0, 2), rep(1, 7), rep(2, 10)))
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
# Wave 4 starts at 4,377 reported cases in week 1.
# Data-implied I0 = 4377 / (0.30 × 2.41) ≈ 6054
# Data-implied E0 = 6054 × (2.41/1.37) ≈ 10650

seir_rinit <- Csnippet("
  double I_init = 6054.0;
  double E_init = 10650.0;
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
#   — eta = 0.55: ~45% immune (vaccination + prior infection by Jul 2021)

params_est <- c("Beta1", "Beta2", "Beta3", "rho", "k")

covid_w4 <- meas_w4 |>
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
# Beta1 = 14.0/wk → R0_Delta = 5.81; consistent with Delta R0 estimates
#                   of 5-6 (Liu & Rocklöv 2021, J Travel Med)
# Beta2 = 3.0/wk  → Rt=1.24 post-peak; rapid decline July-August
# Beta3 = 4.0/wk  → Rt=1.66 plateau phase; sustained transmission Sep-Oct
#                   NOTE: Beta3 > Beta2 — plateau is LESS suppressed than peak decline
# rho   = 0.30    → highest reporting fraction; ~30% of infections captured
# eta   = 0.55    → FIXED: ~45% immune (vaccination + natural immunity)

theta_start <- c(
  N     = 17400000,
  Beta1 = 14.0,
  Beta2 = 3.0,
  Beta3 = 4.0,
  mu_EI = 1.37,
  mu_IR = 2.41,
  eta   = 0.55,
  rho   = 0.30,
  k     = 10
)

coef(covid_w4) <- theta_start

# ─────────────────────────────────────────────────────────────────────────────
# 8.  SIMULATION CHECK
# ─────────────────────────────────────────────────────────────────────────────
# Wave 4 simulations show explosive rise to peak by week 2-3, then decline
# to a PLATEAU (not near-zero). The plateau is ~10-17k/week.
# If simulations decline to zero after week 5: Beta3 is too low.
# If the peak is much higher than 68k: Beta1 needs reducing.

set.seed(2024)
covid_w4 |>
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
  labs(title    = "Wave 4: Simulation check",
       subtitle = "Explosive Delta rise → peak wk 3 → decline to PLATEAU (~10-17k), not zero",
       x = "Week", y = "Weekly reported cases") +
  theme_bw(base_size = 11)

# ─────────────────────────────────────────────────────────────────────────────
# 9.  PRIOR SPECIFICATION
# ─────────────────────────────────────────────────────────────────────────────
# Beta1: log-N(2.64, 0.40²)  → centred on log(14); Delta R0 range
# Beta2: log-N(1.10, 0.50²)  → centred on log(3.0); post-peak suppression
# Beta3: log-N(1.39, 0.50²)  → centred on log(4.0); plateau transmission
# rho:   logit-N(-0.85, 0.60²) → centred on logit(0.30), 95% mass on (0.14,0.52)
# k:     log-N(2.30, 0.80²)

covid_dprior <- Csnippet("
  double lB1  = log(Beta1);
  double lB2  = log(Beta2);
  double lB3  = log(Beta3);
  double lrho = log(rho / (1.0 - rho));
  double lk   = log(k);

  lik = dnorm(lB1,  2.64, 0.40, 1)
      + dnorm(lB2,  1.10, 0.50, 1)
      + dnorm(lB3,  1.39, 0.50, 1)
      + dnorm(lrho,-0.85, 0.60, 1)
      + dnorm(lk,   2.30, 0.80, 1);

  if (!give_log) lik = exp(lik);
")

covid_w4 <- pomp(covid_w4, dprior = covid_dprior,
                 paramnames = c("N","Beta1","Beta2","Beta3",
                                "mu_EI","mu_IR","eta","k","rho"))

# ─────────────────────────────────────────────────────────────────────────────
# 10.  PARTICLE FILTER DIAGNOSTIC
# ─────────────────────────────────────────────────────────────────────────────

Np_use <- 3000

set.seed(999)
ll_check <- replicate(10, logLik(pfilter(covid_w4, Np = Np_use)))
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
  covid_w4,
  Nmcmc    = 5000,
  Np       = Np_use,
  proposal = mvn_diag_rw(c(Beta1=0.15, Beta2=0.15, Beta3=0.15,
                            rho=0.15, k=0.20)^2)
)
acc_p1 <- pmcmc_p1@accepts / 5000
cat(sprintf("Phase 1 acceptance rate: %.3f\n", acc_p1))

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
  mutate(R0_delta  = Beta1 / theta_start["mu_IR"],
         Rt_post   = Beta2 / theta_start["mu_IR"],
         Rt_plateau = Beta3 / theta_start["mu_IR"])

cat("\nPosterior means:\n"); print(round(colMeans(post_df[, params_est]), 4))
cat("\n95% CIs:\n"); print(round(apply(post_df[, params_est], 2, quantile, c(0.025, 0.975)), 4))
for (v in c("R0_delta","Rt_post","Rt_plateau"))
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
  labs(title = "Wave 4: Posterior marginal distributions",
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
  labs(title = "Wave 4: Pairwise posterior scatter")

# ─────────────────────────────────────────────────────────────────────────────
# 17.  POSTERIOR PREDICTIVE CHECK
# ─────────────────────────────────────────────────────────────────────────────

set.seed(42)
pp_sims <- lapply(sample(nrow(post_df), 300), function(i) {
  th            <- coef(covid_w4)
  th[params_est] <- as.numeric(post_df[i, params_est])
  simulate(covid_w4, params=th, nsim=1, format="data.frame") |>
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
  geom_point(data=as.data.frame(covid_w4)|>select(week,reports),
             aes(x=week,y=reports), colour="black", size=2) +
  scale_y_continuous(labels=scales::comma) +
  labs(title    = "Wave 4: Posterior predictive check",
       subtitle = "Shaded=95% interval; line=median; dots=data",
       x="Week", y="Weekly reported cases") +
  theme_bw(base_size=11)
