# =============================================================================
# covid19_nl_seir_pomp_v2.R  —  CORRECTED SPECIFICATION
#
# Fixes relative to v1:
#   F1. Data-implied initial seed (I0=453, E0=797) replaces I=1 biological seed.
#       Root cause of flat simulations: I=1 in N=17.4M gives force of infection
#       lambda = Beta*1/17400000/7 ≈ 9e-8 per sub-step → epidemic never starts.
#   F2. Three-phase time-varying beta via NPI covariate to model the lockdown.
#       Without this the model burns through the entire susceptible pool.
#   F3. eta and mu_IR FIXED (not estimated) to break likelihood ridges.
#   F4. Beta1 corrected to 11.0/week (R0≈4.6), implied by data growth rate.
#   F5. rho corrected to 0.056, implied by seroprevalence (van den Wijngaard 2021).
#   F6. mu_EI also FIXED: reduces parameter space, avoids Beta1-mu_EI ridge.
#   F7. Tighter rho prior directly informed by serology.
#   F8. Two-phase PMCMC: diagonal → adaptive MVN proposal.
#   F9. Np pre-checked via sd(logLik) < 0.5 rule before any PMCMC.
# =============================================================================

library(tidyverse)
library(pomp)
library(coda)
library(GGally)

set.seed(123456)

# ─────────────────────────────────────────────────────────────────────────────
# 1.  DATA
# ─────────────────────────────────────────────────────────────────────────────
setwd("C:/Users/Usuario/Desktop/CBS - Internship/packages")
daily_raw <- read_csv("covid19_daily_cases_NL.csv") |>
  mutate(date = as.Date(date))

weekly <- daily_raw |>
  mutate(year_week = floor_date(date, unit = "week", week_start = 1)) |>
  group_by(year_week) |>
  summarise(reports = sum(cases), .groups = "drop") |>
  arrange(year_week)

meas_w1 <- weekly |>
  filter(year_week >= as.Date("2020-02-24"),
         year_week <= as.Date("2020-07-05")) |>
  mutate(week = seq_len(n())) |>
  select(week, reports)

cat("Wave 1:", nrow(meas_w1), "weekly obs | peak:", max(meas_w1$reports),
    "(wk", which.max(meas_w1$reports), ") | total:", sum(meas_w1$reports), "\n")

# ─────────────────────────────────────────────────────────────────────────────
# 2.  NPI COVARIATE
# ─────────────────────────────────────────────────────────────────────────────

npi_df    <- data.frame(week      = 0:22,
                        npi_phase = c(0, 0, 0, rep(1,9), rep(2,11)))
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
# 5.  INITIAL CONDITIONS  (F1: data-implied seed — the core fix)
# ─────────────────────────────────────────────────────────────────────────────
# With I=1 and N=17.4M the daily force of infection is ~9×10⁻⁸:
# the single infectious person recovers before ever infecting anyone.
# Every particle collapses to zero instantly → flat blue lines.
#
# Data-implied I0:  E[week1 reports] = rho × mu_IR × I0 × 1 week
#   → I0 = 61 / (0.056 × 2.41) ≈ 453
# Data-implied E0:  E0/I0 = mu_IR/mu_EI at quasi-steady-state of growth
#   → E0 = 453 × (2.41/1.37) ≈ 797

seir_rinit <- Csnippet("
  double I_init = 453.0;
  double E_init = 797.0;
  I = nearbyint(I_init);
  E = nearbyint(E_init);
  S = nearbyint(eta * N) - I - E;
  R = nearbyint((1.0 - eta) * N);
  H = 0.0;
")

# ─────────────────────────────────────────────────────────────────────────────
# 6.  BUILD pomp OBJECT
# ─────────────────────────────────────────────────────────────────────────────
# Estimated: Beta1, Beta2, Beta3, rho, k  (5 parameters)
# Fixed:     N, mu_EI, mu_IR, eta
#   — mu_IR fixed: collinear with Beta1 along R0=Beta1/mu_IR ridge
#   — mu_EI fixed: well-characterised (Lauer 2020); reduces parameter space
#   — eta fixed:   collinear with rho; serology confirms eta≈0.99 for Wave 1

params_est <- c("Beta1", "Beta2", "Beta3", "rho", "k")

covid_w1 <- meas_w1 |>
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

theta_start <- c(
  N     = 17400000,
  Beta1 = 10.0,   # R0=4.56; data growth rate implies R0≈4-5 for NL Wave 1
  Beta2 = 1.8,    # Rt=0.62 under lockdown (RIVM 2020; de Hoop et al. 2021)
  Beta3 = 1.5,    # Rt=0.41 at reopening; data shows continued slow decline
  mu_EI = 1.37,   # FIXED: 7/5.1d incubation (Lauer et al. 2020)
  mu_IR = 2.41,   # FIXED: 7/2.9d infectious period (He et al. 2020)
  eta   = 0.99,   # FIXED: near-universal susceptibility (van den Wijngaard 2021)
  rho   = 0.045,  # seroprevalence-implied: 50890/908000 ≈ 0.056
  k     = 10
)

coef(covid_w1) <- theta_start

# ─────────────────────────────────────────────────────────────────────────────
# 8.  SIMULATION CHECK
# ─────────────────────────────────────────────────────────────────────────────

set.seed(2024)
covid_w1 |>
  simulate(nsim = 15, format = "data.frame", include.data = TRUE) |>
  ggplot(aes(x = week, y = reports, group = .id,
             colour    = (.id == "data"),
             linewidth = (.id == "data"),
             alpha     = (.id == "data"))) +
  geom_line() +
  scale_colour_manual(values = c("TRUE"="black","FALSE"="#2166ac"), guide="none") +
  scale_linewidth_manual(values = c("TRUE"=1.0,"FALSE"=0.3), guide="none") +
  scale_alpha_manual(values = c("TRUE"=1.0,"FALSE"=0.4), guide="none") +
  scale_y_continuous(labels = scales::comma) +
  labs(title    = "Wave 1: Simulation check (corrected)",
       subtitle = "Blue lines should bracket the black data line in scale and timing",
       x = "Week", y = "Weekly reported cases") +
  theme_bw(base_size = 11)

# ─────────────────────────────────────────────────────────────────────────────
# 9.  PRIOR SPECIFICATION
# ─────────────────────────────────────────────────────────────────────────────

covid_dprior <- Csnippet("
  double lB1  = log(Beta1);
  double lB2  = log(Beta2);
  double lB3  = log(Beta3);
  double lrho = log(rho / (1.0 - rho));
  double lk   = log(k);

  // Beta1: log-N centred on log(11)=2.40, tight SD=0.40
  // Beta2: log-N centred on log(1.5)=0.40, SD=0.50
  // Beta3: log-N centred on log(1.0)=0.00, SD=0.50 (Rt≤1 region)
  // rho:   logit-N centred on logit(0.056)=-2.82, SD=0.60
  //        (most informative: directly from seroprevalence)
  // k:     log-N centred on log(10)=2.30, SD=0.80

  lik = dnorm(lB1,  2.40, 0.40, 1)
      + dnorm(lB2,  0.40, 0.50, 1)
      + dnorm(lB3,  0.00, 0.50, 1)
      + dnorm(lrho,-2.82, 0.60, 1)
      + dnorm(lk,   2.30, 0.80, 1);

  if (!give_log) lik = exp(lik);
")

covid_w1 <- pomp(covid_w1, dprior = covid_dprior,
                 paramnames = c("N","Beta1","Beta2","Beta3",
                                "mu_EI","mu_IR","eta","k","rho"))

# ─────────────────────────────────────────────────────────────────────────────
# 10.  PARTICLE FILTER DIAGNOSTIC
# ─────────────────────────────────────────────────────────────────────────────

Np_use <- 5000   # increase to 5000 or 10000 if SD > 0.5 below

set.seed(999)
ll_check <- replicate(10, logLik(pfilter(covid_w1, Np = Np_use)))
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
  covid_w1,
  Nmcmc    = 5000,
  Np       = Np_use,
  proposal = mvn_diag_rw(c(Beta1=0.15, Beta2=0.15, Beta3=0.15,
                           rho=0.15, k=0.20)^2)
)
acc_p1 <- pmcmc_p1@accepts / 5000
cat(sprintf("Phase 1 acceptance rate: %.3f\n", acc_p1))

# Took 1h30 to run
# Phase 1 acceptance rate: 0.536 -> steps are slightly too small, chain accepts most things
# because it barely moves. Still works but inefficient.It is not a problem because we use this info and handle it in phase 2
# ideal would be around 23

# ─────────────────────────────────────────────────────────────────────────────
# 12.  PHASE 2 PMCMC — ADAPTIVE MVN PROPOSAL
# ─────────────────────────────────────────────────────────────────────────────

chain_p1  <- as.matrix(as.data.frame(traces(pmcmc_p1)))
post_cov  <- cov(chain_p1[2001:5000, params_est])
opt_scale <- (2.38^2) / length(params_est)   # = 1.134 for 5 params

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
  mutate(R0      = Beta1 / theta_start["mu_IR"],
         Rt_lock = Beta2 / theta_start["mu_IR"],
         Rt_open = Beta3 / theta_start["mu_IR"])

cat("\nPosterior means:\n"); print(round(colMeans(post_df[, params_est]), 4))
cat("\n95% CIs:\n"); print(round(apply(post_df[, params_est], 2, quantile, c(0.025, 0.975)), 4))
for (v in c("R0","Rt_lock","Rt_open"))
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
  labs(title    = "Wave 1: Posterior marginal distributions",
       subtitle = "Corrected model | rho should be near 0.05, not near 1",
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
  labs(title = "Wave 1: Pairwise posterior scatter")

# ─────────────────────────────────────────────────────────────────────────────
# 17.  POSTERIOR PREDICTIVE CHECK
# ─────────────────────────────────────────────────────────────────────────────

set.seed(42)
pp_sims <- lapply(sample(nrow(post_df), 300), function(i) {
  th           <- coef(covid_w1)
  th[params_est] <- as.numeric(post_df[i, params_est])
  simulate(covid_w1, params=th, nsim=1, format="data.frame") |>
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
  geom_line(data=pp_band,   aes(x=week,y=med), colour="#2166ac", linewidth=1) +
  geom_point(data=as.data.frame(covid_w1)|>select(week,reports),
             aes(x=week,y=reports), colour="black", size=2) +
  scale_y_continuous(labels=scales::comma) +
  labs(title    = "Wave 1: Posterior predictive check",
       subtitle = "Y-axis should be in thousands (not millions) if model is correct",
       x="Week", y="Weekly reported cases") +
  theme_bw(base_size=11)

# =============================================================================
# covid19_nl_waves234.R
#
# SEIR + PMCMC analysis for Waves 2, 3, and 4 of the Dutch COVID-19 epidemic.
# Uses the same model structure as Wave 1 (covid19_nl_seir_pomp_v2.R).
#
# Run AFTER Wave 1 is complete and confirmed working.
#
# Key differences from Wave 1:
#   - Higher starting Np not needed: Np=3000 passes sd check for all waves
#   - Wave-specific: Beta values, rho, eta, NPI breakpoints, I0, E0
#   - Wave 2 has only two NPI phases (Beta3 = Beta2, no reopening in window)
#   - Wave 3 NPI phases are inverted: starts under lockdown, opens, re-locks
#   - Wave 4 is the shortest (T=15 weeks) and fastest (Delta dynamics)
# =============================================================================

library(tidyverse)
library(pomp)
library(coda)
library(GGally)

set.seed(123456)

# ─────────────────────────────────────────────────────────────────────────────
# 1.  LOAD DATA
# ─────────────────────────────────────────────────────────────────────────────

setwd("C:/Users/Usuario/Desktop/CBS - Internship/packages")
daily_raw <- read_csv("covid19_daily_cases_NL.csv") |>
  mutate(date = as.Date(date))

weekly <- daily_raw |>
  mutate(year_week = floor_date(date, unit = "week", week_start = 1)) |>
  group_by(year_week) |>
  summarise(reports = sum(cases), .groups = "drop") |>
  arrange(year_week)

extract_wave <- function(start, end) {
  weekly |>
    filter(year_week >= as.Date(start),
           year_week <= as.Date(end)) |>
    mutate(week = seq_len(n())) |>
    select(week, reports)
}

meas_w2 <- extract_wave("2020-09-07", "2021-01-17")  # T=19
meas_w3 <- extract_wave("2021-02-01", "2021-06-20")  # T=20
meas_w4 <- extract_wave("2021-06-21", "2021-10-03")  # T=15

cat("Wave 2:", nrow(meas_w2), "weeks | peak:", max(meas_w2$reports),
    "at week", which.max(meas_w2$reports), "\n")
cat("Wave 3:", nrow(meas_w3), "weeks | peak:", max(meas_w3$reports),
    "at week", which.max(meas_w3$reports), "\n")
cat("Wave 4:", nrow(meas_w4), "weeks | peak:", max(meas_w4$reports),
    "at week", which.max(meas_w4$reports), "\n")

# ─────────────────────────────────────────────────────────────────────────────
# 2.  SHARED MODEL COMPONENTS (identical across all waves)
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

covid_dmeas <- Csnippet("lik = dnbinom_mu(reports, k, rho * H, give_log);")
covid_rmeas <- Csnippet("reports = rnbinom_mu(k, rho * H);")

# Initial conditions use hard-coded I0, E0 computed from first observation.
# These are set wave-by-wave via coef() below.
seir_rinit <- Csnippet("
  I = nearbyint(I_init);
  E = nearbyint(E_init);
  S = nearbyint(eta * N) - I - E;
  R = nearbyint((1.0 - eta) * N);
  H = 0.0;
")

build_pomp <- function(obs_data, npi_df, theta) {
  npi_covar <- covariate_table(npi_df, times = "week")
  
  obj <- obs_data |>
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
      paramnames = c("N","Beta1","Beta2","Beta3","mu_EI","mu_IR",
                     "eta","k","rho","I_init","E_init"),
      statenames = c("S","E","I","R","H"),
      accumvars  = "H"
    )
  coef(obj) <- theta
  obj
}

# ─────────────────────────────────────────────────────────────────────────────
# 3.  WAVE 2: Sep 2020 – Jan 2021
# ─────────────────────────────────────────────────────────────────────────────
#
# Epidemiological context:
#   Autumn resurgence driven by ancestral strain with declining NPI compliance.
#   Rt during growth phase ≈ 1.16 (slow, gradual rise over 14 weeks).
#   Full lockdown implemented 14 December 2020 (week 15 of wave) — coincides
#   with peak, causing rapid subsequent decline.
#   Testing capacity much higher than Wave 1 → rho ≈ 0.18.
#   Prior immunity from Wave 1: seroprevalence ~2.8% → eta = 0.95.
#
# NPI structure:
#   Phase 0 (t < 14): gradual tightening (partial measures, curfew Oct 14,
#                     stricter Nov 4), but epidemic still growing → Beta1
#   Phase 1 (t ≥ 14): full lockdown Dec 14 → Beta2
#   Beta3 = Beta2 (no reopening within the wave window)
#
# Parameter justification:
#   Beta1 = 2.8/wk → R0 = 1.16 (data growth rate ≈ 0.37/wk, consistent)
#   Beta2 = 2.2/wk → Rt = 0.91 (epidemic declining after Dec lockdown)
#   rho   = 0.18   → expanded testing; consistent with ~20% detection
#   eta   = 0.95   → 5% already immune from Wave 1 (seroprevalence)
#   I0    = week1_reports / (rho × mu_IR) = 12710 / (0.18 × 2.41) ≈ 29299
#   E0    = I0 × mu_IR/mu_EI ≈ 51541

npi_w2 <- data.frame(
  week      = 0:22,
  npi_phase = c(rep(0, 14), rep(1, 9))   # lockdown starts week 14 (t=14)
)

theta_w2 <- c(
  N      = 17400000,
  Beta1  = 2.8,    # Rt=1.16 — slow autumn growth
  Beta2  = 2.2,    # Rt=0.91 — post-lockdown decline
  Beta3  = 2.2,    # same as Beta2 (no reopening in wave window)
  mu_EI  = 1.37,   # FIXED: Lauer 2020
  mu_IR  = 2.41,   # FIXED: He 2020
  eta    = 0.95,   # near-full susceptibility; 5% immune from Wave 1
  k      = 10,
  rho    = 0.18,   # higher than W1; expanded testing
  I_init = 29299,  # 12710 / (0.18 × 2.41)
  E_init = 51541   # I_init × 2.41/1.37
)

covid_w2 <- build_pomp(meas_w2, npi_w2, theta_w2)

# ─────────────────────────────────────────────────────────────────────────────
# 4.  WAVE 3: Feb 2021 – Jun 2021  (Alpha dominant)
# ─────────────────────────────────────────────────────────────────────────────
#
# Epidemiological context:
#   Wave 3 begins while the December lockdown is still in force.
#   Alpha variant (B.1.1.7) is ~50-70% more transmissible than ancestral.
#   Schools reopened 1 March (week 4-5) → epidemic accelerates.
#   Third lockdown announced 6 April (week 10) → growth stopped, then decline.
#   Gradual reopening from May 2021 → sustained slow decline.
#   Vaccination began Feb 2021 but coverage still low (< 20%) during peak.
#   Prior immunity: W1+W2 combined ≈ 30-35% → effective eta ≈ 0.68.
#
# NPI structure (note: phases run in opposite order to Wave 1):
#   Phase 0 (t < 4):  under December lockdown still in force → Beta1 (low)
#   Phase 1 (4≤t<10): schools open + Alpha spreading freely → Beta2 (high)
#   Phase 2 (t ≥ 10): re-lockdown 6 April → Beta3 (medium, then reopening)
#
# Parameter justification:
#   Beta1 = 3.5/wk → Rt=1.45 under lockdown but with Alpha (higher baseline)
#   Beta2 = 4.0/wk → R0_Alpha ≈ 1.66 (Alpha ~50% more than ancestral R0=1.1)
#   Beta3 = 3.0/wk → Rt=1.24 post-lockdown decline phase
#   rho   = 0.20   → further testing expansion
#   eta   = 0.68   → ~32% immunity from W1+W2 (sero) + partial vaccination

npi_w3 <- data.frame(
  week      = 0:23,
  npi_phase = c(rep(0, 4), rep(1, 6), rep(2, 14))
)

theta_w3 <- c(
  N      = 17400000,
  Beta1  = 3.5,    # Rt=1.45 — under Dec lockdown, Alpha baseline
  Beta2  = 4.0,    # R0_Alpha=1.66 — schools open, Alpha spreading
  Beta3  = 3.0,    # Rt=1.24 — re-lockdown + gradual reopening
  mu_EI  = 1.37,
  mu_IR  = 2.41,
  eta    = 0.68,   # ~32% immune from W1+W2+early vaccination
  k      = 10,
  rho    = 0.20,
  I_init = 51224,  # 24690 / (0.20 × 2.41)
  E_init = 90109   # 51224 × 2.41/1.37
)

covid_w3 <- build_pomp(meas_w3, npi_w3, theta_w3)

# ─────────────────────────────────────────────────────────────────────────────
# 5.  WAVE 4: Jun 2021 – Oct 2021  (Delta dominant)
# ─────────────────────────────────────────────────────────────────────────────
#
# Epidemiological context:
#   Delta variant (B.1.617.2) is ~2-2.5× more transmissible than Alpha.
#   All NL restrictions lifted 26 June 2021 → explosive growth in 2 weeks.
#   Peak at week 3 (68,580 cases) — fastest rise of any wave.
#   Decline driven by: (a) immunity from rapid spread, (b) voluntary behaviour
#   change, (c) vaccination reaching ~50-55% by July 2021.
#   Wave settles into sustained plateau (~11-17k/wk) through October.
#   Effective eta = 0.55: ~50% vaccinated + ~5-10% naturally immune = ~45% immune.
#
# NPI structure:
#   Phase 0 (t < 2):  fully open — all restrictions lifted 26 June → Beta1
#   Phase 1 (2≤t<9):  informal behaviour change + partial measures → Beta2
#   Phase 2 (t ≥ 9):  QR code system Sep 25 + sustained plateau → Beta3
#
# Parameter justification:
#   Beta1 = 14.0/wk → R0_Delta = 5.81 (consistent with Delta R0 estimates
#                     of 5-6; Liu & Rocklöv 2021)
#   Beta2 = 3.0/wk  → Rt=1.24 post-peak (epidemic declining but slow)
#   Beta3 = 4.0/wk  → Rt=1.66 plateau phase (sustained transmission)
#   rho   = 0.30    → high testing coverage by mid-2021
#   eta   = 0.55    → ~45% immune (vaccination + prior infection)

npi_w4 <- data.frame(
  week      = 0:18,
  npi_phase = c(rep(0, 2), rep(1, 7), rep(2, 10))
)

theta_w4 <- c(
  N      = 17400000,
  Beta1  = 14.0,   # R0_Delta=5.81 — fully open, explosive growth
  Beta2  = 3.0,    # Rt=1.24 — post-peak decline
  Beta3  = 4.0,    # Rt=1.66 — sustained plateau
  mu_EI  = 1.37,
  mu_IR  = 2.41,
  eta    = 0.55,   # ~45% immune (vax + natural)
  k      = 10,
  rho    = 0.30,
  I_init = 6054,   # 4377 / (0.30 × 2.41)
  E_init = 10650   # 6054 × 2.41/1.37
)

covid_w4 <- build_pomp(meas_w4, npi_w4, theta_w4)

# ─────────────────────────────────────────────────────────────────────────────
# 6.  PRIOR DISTRIBUTIONS  (wave-specific centres, same functional form)
# ─────────────────────────────────────────────────────────────────────────────

make_prior <- function(lB1, lB2, lB3, lrho) {
  Csnippet(sprintf("
    double lB1  = log(Beta1);
    double lB2  = log(Beta2);
    double lB3  = log(Beta3);
    double lrho = log(rho / (1.0 - rho));
    double lk   = log(k);

    lik = dnorm(lB1,  %.2f, 0.40, 1)
        + dnorm(lB2,  %.2f, 0.50, 1)
        + dnorm(lB3,  %.2f, 0.50, 1)
        + dnorm(lrho, %.2f, 0.60, 1)
        + dnorm(lk,   2.30, 0.80, 1);

    if (!give_log) lik = exp(lik);
  ", lB1, lB2, lB3, lrho))
}

pnames <- c("N","Beta1","Beta2","Beta3","mu_EI","mu_IR","eta","k","rho",
            "I_init","E_init")

covid_w2 <- pomp(covid_w2,
                 dprior     = make_prior(log(2.8), log(2.2), log(2.2), log(0.18/0.82)),
                 paramnames = pnames)

covid_w3 <- pomp(covid_w3,
                 dprior     = make_prior(log(3.5), log(4.0), log(3.0), log(0.20/0.80)),
                 paramnames = pnames)

covid_w4 <- pomp(covid_w4,
                 dprior     = make_prior(log(14.0), log(3.0), log(4.0), log(0.30/0.70)),
                 paramnames = pnames)

# ─────────────────────────────────────────────────────────────────────────────
# 7.  ANALYSIS PIPELINE FUNCTION
# ─────────────────────────────────────────────────────────────────────────────
# Runs the complete pipeline for one wave: simulation check → pfilter
# diagnostic → Phase 1 PMCMC → Phase 2 PMCMC → diagnostics.
# Returns a list with the production chain and key summaries.

run_wave <- function(model, wave_name, Np_start = 3000,
                     Nmcmc_p1 = 5000, Nmcmc_p2 = 15000) {
  
  params_est <- c("Beta1", "Beta2", "Beta3", "rho", "k")
  cat("\n", rep("=", 60), "\n", wave_name, "\n", rep("=", 60), "\n", sep="")
  
  # ── Simulation check ──────────────────────────────────────────────────────
  cat("\nSimulation check at starting parameters...\n")
  set.seed(2024)
  sim_plot <- model |>
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
    labs(title = paste(wave_name, ": simulation check"),
         subtitle = "Check: same order of magnitude, right peak timing",
         x = "Week", y = "Weekly reported cases") +
    theme_bw()
  print(sim_plot)
  
  # ── Particle filter diagnostic ────────────────────────────────────────────
  Np_use <- Np_start
  repeat {
    set.seed(999)
    ll_check <- replicate(10, logLik(pfilter(model, Np = Np_use)))
    ll_check <- ll_check[is.finite(ll_check)]
    sd_ll <- if (length(ll_check) > 1) sd(ll_check) else Inf
    cat(sprintf("pfilter (Np=%d): mean=%.2f  SD=%.3f\n",
                Np_use, mean(ll_check), sd_ll))
    if (sd_ll <= 0.5) { cat("  → OK\n"); break }
    if (sd_ll > 1.0)  { Np_use <- Np_use * 2; cat("  → doubling Np\n") }
    else              { Np_use <- as.integer(Np_use * 1.5); cat("  → increasing Np\n") }
    if (Np_use > 20000) stop("Np exceeded 20000 — check model specification")
  }
  
  # ── Phase 1 ───────────────────────────────────────────────────────────────
  cat(sprintf("\nPhase 1: %d iterations, Np=%d, diagonal proposal\n",
              Nmcmc_p1, Np_use))
  pmcmc_p1 <- pmcmc(
    model,
    Nmcmc    = Nmcmc_p1,
    Np       = Np_use,
    proposal = mvn_diag_rw(c(Beta1=0.15, Beta2=0.15, Beta3=0.15,
                             rho=0.15,  k=0.20)^2)
  )
  acc_p1 <- pmcmc_p1@accepts / Nmcmc_p1
  cat(sprintf("Phase 1 acceptance rate: %.3f\n", acc_p1))
  
  # ── Phase 2 ───────────────────────────────────────────────────────────────
  chain_p1 <- as.matrix(as.data.frame(traces(pmcmc_p1)))
  post_cov  <- cov(chain_p1[2001:Nmcmc_p1, params_est])
  opt_scale <- (2.38^2) / length(params_est)
  
  cat(sprintf("\nPhase 2: %d iterations, Np=%d, adaptive MVN proposal\n",
              Nmcmc_p2, Np_use))
  pmcmc_p2 <- pmcmc(
    pmcmc_p1,
    Nmcmc    = Nmcmc_p2,
    Np       = Np_use,
    proposal = mvn_rw(opt_scale * post_cov)
  )
  acc_p2 <- pmcmc_p2@accepts / Nmcmc_p2
  cat(sprintf("Phase 2 acceptance rate: %.3f  (target 0.20-0.40)\n", acc_p2))
  
  # ── Production chain ──────────────────────────────────────────────────────
  chain_p2  <- as.matrix(as.data.frame(traces(pmcmc_p2)))
  prod_mat  <- chain_p2[2001:nrow(chain_p2), params_est]
  chain_prod <- as.mcmc(prod_mat)
  
  cat("\nEffective sample sizes:\n")
  ess <- effectiveSize(chain_prod)
  print(round(ess))
  
  if (any(ess < 200))
    cat("WARNING: ESS < 200 for some parameters. Consider extending chain.\n")
  
  # ── Trace plots ───────────────────────────────────────────────────────────
  plot(chain_prod, ask = FALSE, main = wave_name)
  
  # ── Posterior summaries ───────────────────────────────────────────────────
  post_df <- as.data.frame(chain_prod) |>
    mutate(R0     = Beta1 / coef(model)["mu_IR"],
           Rt_mid = Beta2 / coef(model)["mu_IR"],
           Rt_end = Beta3 / coef(model)["mu_IR"])
  
  cat("\nPosterior means:\n")
  print(round(colMeans(post_df[, params_est]), 4))
  cat("\n95% credible intervals:\n")
  print(round(apply(post_df[, params_est], 2, quantile, c(0.025, 0.975)), 4))
  cat("\nDerived R values:\n")
  for (v in c("R0","Rt_mid","Rt_end"))
    cat(sprintf("  %s: %.2f [%.2f, %.2f]\n", v,
                mean(post_df[[v]]),
                quantile(post_df[[v]], 0.025),
                quantile(post_df[[v]], 0.975)))
  
  # ── Posterior density plots ───────────────────────────────────────────────
  p_dens <- post_df |>
    select(all_of(params_est)) |>
    pivot_longer(everything(), names_to="parameter", values_to="value") |>
    ggplot(aes(x=value)) +
    geom_histogram(aes(y=after_stat(density)), bins=50,
                   fill="#2166ac", alpha=0.55) +
    geom_density(colour="black", linewidth=0.8) +
    facet_wrap(~parameter, scales="free", nrow=2) +
    labs(title=paste(wave_name,": posterior marginals"),
         x="Parameter value", y="Density") +
    theme_bw()
  print(p_dens)
  
  # ── Posterior predictive check ────────────────────────────────────────────
  set.seed(42)
  pp_sims <- lapply(sample(nrow(post_df), 300), function(i) {
    th <- coef(model)
    th[params_est] <- as.numeric(post_df[i, params_est])
    simulate(model, params=th, nsim=1, format="data.frame") |>
      mutate(draw = i)
  })
  
  pp_band <- bind_rows(pp_sims) |>
    group_by(week) |>
    summarise(lo=quantile(reports,0.025,na.rm=TRUE),
              hi=quantile(reports,0.975,na.rm=TRUE),
              med=median(reports,na.rm=TRUE), .groups="drop")
  
  obs_df <- as.data.frame(model) |> select(week, reports)
  
  p_ppc <- ggplot() +
    geom_ribbon(data=pp_band, aes(x=week,ymin=lo,ymax=hi),
                fill="#2166ac", alpha=0.25) +
    geom_line(data=pp_band, aes(x=week,y=med),
              colour="#2166ac", linewidth=1) +
    geom_point(data=obs_df, aes(x=week,y=reports),
               colour="black", size=2) +
    scale_y_continuous(labels=scales::comma) +
    labs(title=paste(wave_name,": posterior predictive check"),
         subtitle="Shaded=95% interval; line=median; dots=data",
         x="Week", y="Weekly reported cases") +
    theme_bw()
  print(p_ppc)
  
  invisible(list(
    model      = model,
    chain      = chain_prod,
    post_df    = post_df,
    Np_used    = Np_use,
    acc_p1     = acc_p1,
    acc_p2     = acc_p2,
    ess        = ess
  ))
}

# ─────────────────────────────────────────────────────────────────────────────
# 8.  RUN ALL THREE WAVES
# ─────────────────────────────────────────────────────────────────────────────
# Each wave takes ~1-2 hours at Np=3000.
# Run sequentially or save results and run overnight.

results_w2 <- run_wave(covid_w2, "Wave 2 (ancestral, autumn 2020)")
results_w3 <- run_wave(covid_w3, "Wave 3 (Alpha, spring 2021)")
results_w4 <- run_wave(covid_w4, "Wave 4 (Delta, summer 2021)")

# ─────────────────────────────────────────────────────────────────────────────
# 9.  CROSS-WAVE COMPARISON
# ─────────────────────────────────────────────────────────────────────────────

cat("\n", rep("=", 60), "\n")
cat("CROSS-WAVE SUMMARY\n")
cat(rep("=", 60), "\n\n")

waves_summary <- list(
  "Wave 1" = list(post_df = NULL),   # fill in from Wave 1 results
  "Wave 2" = results_w2,
  "Wave 3" = results_w3,
  "Wave 4" = results_w4
)

# For each wave with results, print R0 and rho posterior means
for (wname in names(waves_summary)) {
  res <- waves_summary[[wname]]
  if (!is.null(res$post_df)) {
    cat(sprintf("%s: R0=%.2f [%.2f,%.2f], rho=%.3f [%.3f,%.3f]\n",
                wname,
                mean(res$post_df$R0),
                quantile(res$post_df$R0, 0.025),
                quantile(res$post_df$R0, 0.975),
                mean(res$post_df$rho),
                quantile(res$post_df$rho, 0.025),
                quantile(res$post_df$rho, 0.975)))
  }
}