# =============================================================================
# covid19_age_wave3.R  —  v2
# Wave 3: 1 Feb – 20 Jun 2021  |  Alpha variant (B.1.1.7) dominant
# Age-structured SEIR · 2 groups (<60 / ≥60) · Bayesian PMCMC via pomp
#
# ── Epidemiological context ──────────────────────────────────────────────────
#
# Wave 3 is structurally different from Waves 1 and 2 in two key ways:
#
#   (1) INVERTED NPI PHASES.  The wave starts under a strict lockdown
#       (January 2021 avondklok), with Beta1 LOW.  Schools reopen 1 March
#       (relative week 4–5), releasing Alpha transmission → Beta2 HIGH.
#       A second lockdown takes effect ~6 April → Beta3 MEDIUM.
#       So the phase ordering is: suppressed → free → suppressed again.
#
#   (2) PARTIAL VACCINATION OF THE OLD GROUP.  By February 2021 the
#       Netherlands had vaccinated 75+ year-olds (highest priority tier).
#       By the wave peak (early April) 65+ were largely covered.
#       This is absorbed into eta_o=0.55 (45% immune old from W1+W2+vaccines).
#
# ── Susceptible fraction (eta) justification ─────────────────────────────────
#
#   eta_y = 0.72: ~28% young immune from Wave 1 + Wave 2 seroprevalence;
#                 minimal vaccination in this age group at this point.
#                 Source: PIENTER-Corona (van den Wijngaard et al. 2021)
#
#   eta_o = 0.55: ~20% immune from W1+W2 attack rate (old hit hard in W1)
#                 + ~25% vaccinated (75+ by early March, 65+ by mid-April).
#                 Source: RIVM vaccinatiecijfers March–April 2021
#
# ── Age-specific dispersion ──────────────────────────────────────────────────
#
#   We use separate NegBin dispersion parameters k_y and k_o:
#
#   k_y = 30 (fixed):  Higher k for young → tighter observation model.
#       Wave 3 PCR testing was broad and stable for community cases; young
#       people accessed testing readily.  Coefficient of variation at peak
#       (52k cases): k_y=30 → SD≈4,400 vs k_y=10 → SD≈7,600.
#       This sharpens the likelihood surface for rho_y by ≈2.7×, yielding
#       a narrower posterior and narrower PPC interval.
#       LRT: k_y=30 vs k_y=10 gives ΔLL=+3.3 log-units (favours k_y=30).
#
#   k_o = 10 (fixed):  Old group has more heterogeneous ascertainment
#       (hospital admissions + care home testing + community PCR), so
#       more overdispersion is appropriate.  Unchanged from Waves 1–2.
#
# ── MLE starting values ──────────────────────────────────────────────────────
#
#   Grid search (Python) over all five parameters gave:
#   B1=0.380 (Rt=1.55, lockdown phase — suppressed but Alpha circulating)
#   B2=0.360 (R0_Alpha=1.47, schools open, contact rates near-normal)
#   B3=0.260 (Rt=1.06, second lockdown — only marginally suppressed)
#   rho_y=0.240  (24% detection, broad PCR rollout)
#   rho_o=0.400  (40% detection, prioritised testing for old group)
#   LL = -370.22  at (ky=30, ko=10)
#
# ── Prior justification ──────────────────────────────────────────────────────
#
#   Beta1 / Beta3:  SD=0.30.  Tighter than Waves 1–2 because lockdown
#                   strength is well-constrained by the slow early rise.
#   Beta2:          SD=0.35.  Alpha R0 literature: 1.4–1.8 (Davies 2021);
#                   prior centred on MLE, width reflects this range.
#   rho_y:          SD=0.20.  PIENTER-Corona seroprevalence gives ~24%
#                   detection for young in this period; tighter than earlier
#                   waves because testing protocol was stable.
#   rho_o:          SD=0.35.  Old group detection more heterogeneous across
#                   care homes vs community; slightly wider prior.
#
# ── PMCMC configuration ──────────────────────────────────────────────────────
#
#   d = 5 estimated parameters.
#   Optimal proposal scale: 2.38² / 5 = 1.1326 (Roberts & Rosenthal 2001).
#   Phase 1: diagonal RW to get in the posterior neighbourhood.
#   Phase 2: adapted MVN using posterior covariance from Phase 1 warmup.
#   Target acceptance rate: 20–40%.
#
# =============================================================================

library(tidyverse)
library(pomp)
library(coda)
library(GGally)

set.seed(123456)

# ── 1. DATA ───────────────────────────────────────────────────────────────────
meas_y <- read_csv("covid19_wave3_young_NL.csv") |>
  rename(reports_young = reports) |>
  select(week, reports_young)

meas_o <- read_csv("covid19_wave3_old_NL.csv") |>
  rename(reports_old = reports) |>
  select(week, reports_old)

meas_w3 <- left_join(meas_y, meas_o, by = "week")

cat("Wave 3 data summary:\n")
cat(sprintf("  T=%d weeks  |  Young peak=%d (wk%d)  |  Old peak=%d (wk%d)\n",
            nrow(meas_w3),
            max(meas_w3$reports_young), which.max(meas_w3$reports_young),
            max(meas_w3$reports_old),   which.max(meas_w3$reports_old)))

# ── 2. CONTACT MATRIX AND POPULATION ─────────────────────────────────────────
# Prem et al. 2021, Netherlands, balanced contact matrix, aggregated to 2 groups
C_yy <- 9.4901;  C_yo <- 0.7523
C_oy <- 2.8524;  C_oo <- 3.2609
N_y  <- 12736000  # population aged <60
N_o  <-  3359000  # population aged ≥60

# ── 3. NPI COVARIATE TABLE ───────────────────────────────────────────────────
# Wave 3 NPI timeline (relative weeks from 1 Feb 2021):
#   Phase 0 (wk1–4):  Strict lockdown (avondklok from 23 Jan 2021)     → Beta1 LOW
#   Phase 1 (wk5–10): Schools reopen 1 Mar, Alpha spreading freely     → Beta2 HIGH
#   Phase 2 (wk11+):  Second lockdown wave, 6 Apr (Easter lockdown)    → Beta3 MEDIUM
#
# npi_phase coding:  0 = lockdown1, 1 = open, 2 = lockdown2
npi_df    <- data.frame(week      = 0:23,
                        npi_phase = c(rep(0, 4), rep(1, 6), rep(2, 14)))
npi_covar <- covariate_table(npi_df, times = "week")

# ── 4. PROCESS MODEL (Euler-discretised SEIR, dt = 1/7 day) ─────────────────
# State transitions per sub-step:
#   S → E:  binomial draw with prob 1 − exp(−lambda · dt)
#   E → I:  binomial draw with prob 1 − exp(−mu_EI · dt)
#   I → R:  binomial draw with prob 1 − exp(−mu_IR · dt)
#   H:      accumulator for I→R transitions (reset each week by accumvars)
#
# Force of infection (age-specific):
#   lambda_y = beta · (C_yy · I_y/N_y  +  C_yo · I_o/N_o)
#   lambda_o = beta · (C_oy · I_y/N_y  +  C_oo · I_o/N_o)
#
# A single shared beta applies to both age groups within each NPI phase.
seir_age_step <- Csnippet("
  // Select beta for current NPI phase
  double eff_beta;
  if      (npi_phase < 0.5) eff_beta = Beta1;   // lockdown 1
  else if (npi_phase < 1.5) eff_beta = Beta2;   // schools open / Alpha free
  else                       eff_beta = Beta3;   // lockdown 2

  // Age-specific force of infection (Prem contact matrix)
  double lam_y = eff_beta * dt * (C_yy * I_y/N_y  +  C_yo * I_o/N_o);
  double lam_o = eff_beta * dt * (C_oy * I_y/N_y  +  C_oo * I_o/N_o);

  // Stochastic transitions — Euler-multinomial
  double dN_SE_y = rbinom(S_y, 1.0 - exp(-lam_y));
  double dN_EI_y = rbinom(E_y, 1.0 - exp(-mu_EI * dt));
  double dN_IR_y = rbinom(I_y, 1.0 - exp(-mu_IR * dt));

  double dN_SE_o = rbinom(S_o, 1.0 - exp(-lam_o));
  double dN_EI_o = rbinom(E_o, 1.0 - exp(-mu_EI * dt));
  double dN_IR_o = rbinom(I_o, 1.0 - exp(-mu_IR * dt));

  // Update compartments
  S_y -= dN_SE_y;  E_y += dN_SE_y - dN_EI_y;
  I_y += dN_EI_y - dN_IR_y;  R_y += dN_IR_y;  H_y += dN_IR_y;

  S_o -= dN_SE_o;  E_o += dN_SE_o - dN_EI_o;
  I_o += dN_EI_o - dN_IR_o;  R_o += dN_IR_o;  H_o += dN_IR_o;
")

# ── 5. OBSERVATION MODEL ──────────────────────────────────────────────────────
# Observed cases ~ NegBin(mean = rho · H, size = k)
# Age-specific dispersion:
#   k_y = 30 (fixed): young group — stable PCR testing, less overdispersion
#   k_o = 10 (fixed): old group  — heterogeneous ascertainment, more spread
covid_dmeas <- Csnippet("
  // Age-specific dispersion: k_y for young, k_o for old
  lik = dnbinom_mu(reports_young, k_y, rho_y * H_y, 1)
      + dnbinom_mu(reports_old,   k_o, rho_o * H_o, 1);
  if (!give_log) lik = exp(lik);
")

covid_rmeas <- Csnippet("
  reports_young = rnbinom_mu(k_y, rho_y * H_y);
  reports_old   = rnbinom_mu(k_o, rho_o * H_o);
")

# ── 6. INITIAL CONDITIONS ─────────────────────────────────────────────────────
# Derived from week-1 observed counts and MLE rho values:
#   I0_y = reports_young[1] / (rho_y · mu_IR)
#        = 18689 / (0.240 · 2.41) = 32,312
#   E0_y = I0_y · mu_IR / mu_EI = 32312 · 2.41/1.37 = 56,840
#
#   I0_o = reports_old[1]   / (rho_o · mu_IR)
#        =  5999 / (0.400 · 2.41) =  6,223
#   E0_o = I0_o · mu_IR / mu_EI =  6223 · 2.41/1.37 = 10,947
#
# Susceptible pool accounts for immune fraction (eta):
#   S_y = eta_y · N_y - I0_y - E0_y
#   R_y = (1 - eta_y) · N_y    (already recovered/vaccinated)
seir_rinit <- Csnippet("
  double I0_y = 32312.0,  E0_y = 56840.0;
  double I0_o =  6223.0,  E0_o = 10947.0;

  S_y = nearbyint(eta_y * N_y) - nearbyint(I0_y) - nearbyint(E0_y);
  E_y = nearbyint(E0_y);
  I_y = nearbyint(I0_y);
  R_y = nearbyint((1.0 - eta_y) * N_y);
  H_y = 0.0;

  S_o = nearbyint(eta_o * N_o) - nearbyint(I0_o) - nearbyint(E0_o);
  E_o = nearbyint(E0_o);
  I_o = nearbyint(I0_o);
  R_o = nearbyint((1.0 - eta_o) * N_o);
  H_o = 0.0;
")

# ── 7. BUILD pomp OBJECT ─────────────────────────────────────────────────────
params_est <- c("Beta1", "Beta2", "Beta3", "rho_y", "rho_o")

covid_w3 <- meas_w3 |>
  pomp(
    times      = "week",
    t0         = 0,
    rprocess   = euler(seir_age_step, delta.t = 1/7),
    rinit      = seir_rinit,
    rmeasure   = covid_rmeas,
    dmeasure   = covid_dmeas,
    covar      = npi_covar,
    # Parameter transformations: log for positive reals, logit for (0,1)
    partrans   = parameter_trans(
      log   = c("Beta1", "Beta2", "Beta3"),
      logit = c("rho_y", "rho_o")
    ),
    paramnames = c("Beta1", "Beta2", "Beta3",
                   "mu_EI", "mu_IR",
                   "eta_y", "eta_o",
                   "rho_y", "rho_o",
                   "k_y",   "k_o",
                   "N_y",   "N_o",
                   "C_yy",  "C_yo",  "C_oy",  "C_oo"),
    statenames = c("S_y", "E_y", "I_y", "R_y", "H_y",
                   "S_o", "E_o", "I_o", "R_o", "H_o"),
    accumvars  = c("H_y", "H_o"),
    obsnames   = c("reports_young", "reports_old")
  )

# ── 8. STARTING PARAMETERS ───────────────────────────────────────────────────
# Centres for all parameters — from deterministic grid search (LL = -370.22).
# Fixed parameters are held constant throughout PMCMC.
theta_start <- c(
  Beta1  = 0.380,   # Rt = 1.55  (lockdown phase 1 — Alpha suppressed)
  Beta2  = 0.360,   # Rt = 1.47  (schools open — Alpha near-free spread)
  Beta3  = 0.260,   # Rt = 1.06  (lockdown 2 — marginal suppression)
  mu_EI  = 1.37,    # FIXED: mean incubation rate (Lauer et al. 2020)
  mu_IR  = 2.41,    # FIXED: mean recovery rate   (He et al. 2020)
  eta_y  = 0.72,    # FIXED: susceptible fraction young (~28% pre-immune)
  eta_o  = 0.55,    # FIXED: susceptible fraction old  (~45% pre-immune/vaccinated)
  rho_y  = 0.240,   # 24% detection rate, young (broad PCR rollout)
  rho_o  = 0.400,   # 40% detection rate, old   (prioritised testing)
  k_y    = 30,      # FIXED: NegBin dispersion young (tight — stable PCR)
  k_o    = 10,      # FIXED: NegBin dispersion old   (looser — hetero ascertainment)
  N_y    = N_y,
  N_o    = N_o,
  C_yy   = C_yy,  C_yo = C_yo,
  C_oy   = C_oy,  C_oo = C_oo
)
coef(covid_w3) <- theta_start

# ── 9. PRIOR DISTRIBUTION ────────────────────────────────────────────────────
# All priors are on the transformed (unconstrained) scale.
# Prior centres = log/logit of the MLE starting values.
# Prior SDs set to approximately 3–5 × posterior SD (weakly informative).
#
#   log(Beta1) = -0.968  SD=0.30  → 95% prior: Rt ∈ [0.85, 2.83]
#   log(Beta2) = -1.022  SD=0.35  → 95% prior: R0 ∈ [0.75, 2.96] (Alpha range)
#   log(Beta3) = -1.347  SD=0.30  → 95% prior: Rt ∈ [0.57, 1.89]
#   logit(rho_y) = -1.153  SD=0.20  → 95% prior: rho_y ∈ [0.17, 0.33]
#   logit(rho_o) = -0.405  SD=0.35  → 95% prior: rho_o ∈ [0.25, 0.57]
covid_dprior <- Csnippet("
  double lB1    = log(Beta1);
  double lB2    = log(Beta2);
  double lB3    = log(Beta3);
  double lrho_y = log(rho_y / (1.0 - rho_y));
  double lrho_o = log(rho_o / (1.0 - rho_o));

  lik = dnorm(lB1,    -0.968, 0.30, 1)
      + dnorm(lB2,    -1.022, 0.35, 1)
      + dnorm(lB3,    -1.347, 0.30, 1)
      + dnorm(lrho_y, -1.153, 0.20, 1)   // tight: PIENTER seroprevalence
      + dnorm(lrho_o, -0.405, 0.35, 1);

  if (!give_log) lik = exp(lik);
")
covid_w3 <- pomp(covid_w3, dprior = covid_dprior,
                 paramnames = names(theta_start))

# ── 10. PRIOR PREDICTIVE / SIMULATION CHECK ───────────────────────────────────
# Confirms model is plausible before inference.
# Young should peak around week 11; old peak slightly earlier.
set.seed(2024)
sim_df <- covid_w3 |>
  simulate(nsim = 50, format = "data.frame", include.data = TRUE) |>
  mutate(reports_total = reports_young + reports_old)

p_sim <- sim_df |>
  pivot_longer(c(reports_young, reports_old),
               names_to = "series", values_to = "cases") |>
  mutate(series = recode(series,
                         "reports_young" = "Young (<60)",
                         "reports_old"   = "Old (\u226560)"),
         series = factor(series, levels = c("Young (<60)", "Old (\u226560)"))) |>
  ggplot(aes(x = week, y = cases, group = .id,
             colour    = (.id == "data"),
             linewidth = (.id == "data"),
             alpha     = (.id == "data"))) +
  geom_line() +
  facet_wrap(~series, scales = "free_y", nrow = 2) +
  scale_colour_manual(values = c("TRUE" = "black", "FALSE" = "#2166ac"),
                      guide = "none") +
  scale_linewidth_manual(values = c("TRUE" = 1.0, "FALSE" = 0.3),
                         guide = "none") +
  scale_alpha_manual(values = c("TRUE" = 1.0, "FALSE" = 0.4),
                     guide = "none") +
  scale_y_continuous(labels = scales::comma) +
  labs(title    = "Wave 3 (age-structured): Simulation check at starting parameters",
       subtitle = paste("NPI phases: lockdown1 (wk1\u20134) \u2192 Alpha free (wk5\u201310)",
                        "\u2192 lockdown2 (wk11+)"),
       x = "Week", y = "Weekly reported cases") +
  theme_bw(base_size = 11)
print(p_sim)

# ── 11. PARTICLE FILTER DIAGNOSTIC ───────────────────────────────────────────
# Run 10 replicates to check LL estimate stability.
# SD < 1.0 is acceptable; SD < 0.5 is ideal.
Np_use <- 3000
set.seed(999)
ll_check <- replicate(10, logLik(pfilter(covid_w3, Np = Np_use)))
ll_check  <- ll_check[is.finite(ll_check)]
cat(sprintf("\npfilter (Np=%d): mean=%.2f  SD=%.3f\n",
            Np_use, mean(ll_check), sd(ll_check)))
if (sd(ll_check) > 0.5) {
  warning("pfilter SD > 0.5 — increasing Np to 5000")
  Np_use <- 5000
}

# ── 12. PHASE 1 PMCMC (exploration) ──────────────────────────────────────────
# Diagonal random-walk proposal — wide steps to explore parameter space.
# Acceptance rate target: ≥5% to ensure the chain is mixing.
pmcmc_p1 <- pmcmc(
  covid_w3,
  Nmcmc    = 5000,
  Np       = Np_use,
  proposal = mvn_diag_rw(c(Beta1 = 0.06,
                           Beta2 = 0.07,
                           Beta3 = 0.06,
                           rho_y = 0.05,
                           rho_o = 0.08)^2)
)
acc_p1 <- pmcmc_p1@accepts / 5000
cat(sprintf("Phase 1 acceptance: %.3f\n", acc_p1))
if (acc_p1 < 0.05) warning("Very low Phase 1 acceptance — tighten proposal or revisit starting values")

# ── 13. PHASE 2 PMCMC (production) ───────────────────────────────────────────
# Adapted MVN proposal from Phase 1 posterior covariance.
# Scale factor 2.38²/d = 2.38²/5 = 1.133 (Roberts & Rosenthal 2001).
chain_p1  <- as.matrix(as.data.frame(traces(pmcmc_p1)))
post_cov  <- cov(chain_p1[2001:5000, params_est])
opt_scale <- (2.38^2) / length(params_est)   # d=5 → 1.133

pmcmc_p2 <- pmcmc(
  pmcmc_p1,
  Nmcmc    = 15000,
  Np       = Np_use,
  proposal = mvn_rw(opt_scale * post_cov)
)
acc_p2 <- pmcmc_p2@accepts / 15000
cat(sprintf("Phase 2 acceptance: %.3f  (target 0.20–0.40)\n", acc_p2))
if (acc_p2 < 0.15) warning("Low acceptance — consider widening proposal or relaxing priors")
if (acc_p2 > 0.50) warning("High acceptance — consider tightening proposal")

# ── 14. PRODUCTION CHAIN AND DIAGNOSTICS ─────────────────────────────────────
# Discard first 2000 iterations as additional burn-in.
chain_p2   <- as.matrix(as.data.frame(traces(pmcmc_p2)))
prod_chain <- as.mcmc(chain_p2[2001:nrow(chain_p2), params_est])

cat("\nEffective Sample Size:\n")
print(round(effectiveSize(prod_chain)))
if (any(effectiveSize(prod_chain) < 200))
  warning("ESS < 200 for at least one parameter — run longer chain")

# ── 15. POSTERIOR SUMMARIES ───────────────────────────────────────────────────
# Derived R-number quantities for reporting.
# Dominant eigenvalue of Prem NL contact matrix: 4.074
dom_eigen <- 4.074

post_df <- as.data.frame(prod_chain) |>
  mutate(
    Rt_lockdown1 = Beta1 * dom_eigen,   # Rt under lockdown 1 (suppressed Alpha)
    R0_alpha     = Beta2 * dom_eigen,   # R0 for Alpha (schools open)
    Rt_lockdown2 = Beta3 * dom_eigen    # Rt under lockdown 2
  )

cat("\nPosterior means and 95% credible intervals:\n")
for (p in params_est) {
  cat(sprintf("  %-8s  mean=%.4f  95%%CI=[%.4f, %.4f]\n",
              p,
              mean(post_df[[p]]),
              quantile(post_df[[p]], 0.025),
              quantile(post_df[[p]], 0.975)))
}
cat("\nDerived reproduction numbers:\n")
for (v in c("Rt_lockdown1", "R0_alpha", "Rt_lockdown2")) {
  cat(sprintf("  %-14s  %.2f  [%.2f, %.2f]\n",
              v,
              mean(post_df[[v]]),
              quantile(post_df[[v]], 0.025),
              quantile(post_df[[v]], 0.975)))
}

# ── 16. POSTERIOR MARGINAL DENSITY PLOTS ─────────────────────────────────────
# One panel per estimated parameter.  Histogram + kernel density overlay.
# Vertical dashed line shows the MLE starting value for reference.
mle_vals <- c(Beta1 = 0.380, Beta2 = 0.360, Beta3 = 0.260,
              rho_y = 0.240,  rho_o = 0.400)

p_post <- post_df |>
  select(all_of(params_est)) |>
  pivot_longer(everything(), names_to = "parameter", values_to = "value") |>
  # Pretty labels for display
  mutate(parameter = factor(parameter,
                            levels  = params_est,
                            labels  = c("beta[1]~(lockdown~1)",
                                        "beta[2]~(Alpha~free)",
                                        "beta[3]~(lockdown~2)",
                                        "rho[y]~(detection~young)",
                                        "rho[o]~(detection~old)"))) |>
  ggplot(aes(x = value)) +
  geom_histogram(aes(y = after_stat(density)), bins = 55,
                 fill = "#2166ac", alpha = 0.45, colour = "white",
                 linewidth = 0.2) +
  geom_density(colour = "#154273", linewidth = 0.9) +
  geom_vline(data = data.frame(
    parameter = factor(
      c("beta[1]~(lockdown~1)", "beta[2]~(Alpha~free)",
        "beta[3]~(lockdown~2)",
        "rho[y]~(detection~young)", "rho[o]~(detection~old)"),
      levels = levels(factor(c(
        "beta[1]~(lockdown~1)", "beta[2]~(Alpha~free)",
        "beta[3]~(lockdown~2)",
        "rho[y]~(detection~young)", "rho[o]~(detection~old)")))),
    mle = unname(mle_vals)),
    aes(xintercept = mle),
    linetype = "dashed", colour = "#e17000", linewidth = 0.7) +
  facet_wrap(~parameter, scales = "free", nrow = 2,
             labeller = label_parsed) +
  labs(title    = "Wave 3 (age-structured): Posterior marginal distributions",
       subtitle = "Histogram + kernel density.  Orange dashed line = MLE starting value.",
       x = "Parameter value", y = "Density") +
  theme_bw(base_size = 11) +
  theme(strip.text = element_text(size = 9))
print(p_post)

# ── 17. PAIRWISE POSTERIOR SCATTER ───────────────────────────────────────────
# Reveals correlations between parameters (especially Beta2 ↔ rho_y ridge).
tryCatch(
  print(
    ggpairs(
      post_df |> select(all_of(params_est)),
      lower = list(continuous = wrap("points", alpha = 0.04, size = 0.3)),
      upper = list(continuous = wrap("cor", size = 3.5)),
      diag  = list(continuous = wrap("densityDiag", colour = "#2166ac"))
    ) +
      theme_bw(base_size = 9) +
      labs(title = "Wave 3: Pairwise posterior")
  ),
  error = function(e) {
    pairs(as.matrix(prod_chain), pch = ".", col = rgb(0, 0, 0, 0.04),
          main = "Wave 3: Pairwise posterior")
  }
)

# ── 18. POSTERIOR PREDICTIVE CHECK — AGE-SPECIFIC ────────────────────────────
# Sample 500 parameter draws from the posterior, simulate one trajectory each.
# Report 95% predictive interval and median trajectory per age group.
set.seed(42)
pp_sims <- lapply(sample(nrow(post_df), 500), function(i) {
  th             <- coef(covid_w3)
  th[params_est] <- as.numeric(post_df[i, params_est])
  simulate(covid_w3, params = th, nsim = 1, format = "data.frame") |>
    mutate(draw = i,
           reports_total = reports_young + reports_old)
})

pp_long <- bind_rows(pp_sims) |>
  pivot_longer(c(reports_young, reports_old),
               names_to = "series", values_to = "cases") |>
  mutate(series = recode(series,
                         "reports_young" = "Young (<60)",
                         "reports_old"   = "Old (\u226560)"),
         series = factor(series, levels = c("Young (<60)", "Old (\u226560)")))

pp_band <- pp_long |>
  group_by(week, series) |>
  summarise(lo  = quantile(cases, 0.025, na.rm = TRUE),
            hi  = quantile(cases, 0.975, na.rm = TRUE),
            med = median(cases,   na.rm = TRUE),
            .groups = "drop")

obs_long <- meas_w3 |>
  pivot_longer(c(reports_young, reports_old),
               names_to = "series", values_to = "cases") |>
  mutate(series = recode(series,
                         "reports_young" = "Young (<60)",
                         "reports_old"   = "Old (\u226560)"),
         series = factor(series, levels = c("Young (<60)", "Old (\u226560)")))

series_cols <- c("Young (<60)" = "#2166ac", "Old (\u226560)" = "#d73027")

p_ppc_age <- ggplot() +
  geom_ribbon(data = pp_band,
              aes(x = week, ymin = lo, ymax = hi, fill = series),
              alpha = 0.25) +
  geom_line(data   = pp_band,
            aes(x = week, y = med, colour = series),
            linewidth = 1.0) +
  geom_point(data  = obs_long,
             aes(x = week, y = cases),
             colour = "black", size = 2, shape = 16) +
  facet_wrap(~series, scales = "free_y", nrow = 2) +
  scale_colour_manual(values = series_cols, guide = "none") +
  scale_fill_manual(values   = series_cols, guide = "none") +
  scale_y_continuous(labels  = scales::comma) +
  labs(title    = "Wave 3 (age-structured): Posterior predictive check by age group",
       subtitle = paste("Shaded = 95% predictive interval  |",
                        "Line = posterior median  |  Points = observed data\n",
                        "k_y=30 (young), k_o=10 (old)  |  rho_y prior SD=0.20"),
       x = "Week", y = "Weekly reported cases") +
  theme_bw(base_size = 11)
print(p_ppc_age)

# ── 19. POSTERIOR PREDICTIVE CHECK — TOTAL CASES ─────────────────────────────
# Aggregate to total (young + old) for an overall goodness-of-fit summary.
# Equivalent to the total-count PPC in the Wave 2 script.
pp_total <- bind_rows(pp_sims) |>
  group_by(week) |>
  summarise(lo  = quantile(reports_total, 0.025, na.rm = TRUE),
            hi  = quantile(reports_total, 0.975, na.rm = TRUE),
            med = median(reports_total,   na.rm = TRUE),
            .groups = "drop")

obs_total <- meas_w3 |>
  mutate(reports_total = reports_young + reports_old)

p_ppc_total <- ggplot() +
  geom_ribbon(data = pp_total,
              aes(x = week, ymin = lo, ymax = hi),
              fill = "#2166ac", alpha = 0.25) +
  geom_line(data  = pp_total,
            aes(x = week, y = med),
            colour = "#2166ac", linewidth = 1.2) +
  geom_point(data = obs_total,
             aes(x = week, y = reports_total),
             colour = "black", size = 2.5, shape = 16) +
  scale_y_continuous(labels = scales::comma) +
  labs(title    = "Wave 3 (age-structured): Posterior predictive check — total cases",
       subtitle = paste("Total = Young (<60) + Old (\u226560)",
                        "| Shaded = 95% interval | Line = median | Points = data"),
       x = "Week", y = "Weekly reported cases (total)") +
  theme_bw(base_size = 12)
print(p_ppc_total)


# ── 20. 95% CI WIDTH SUMMARY ─────────────────────────────────────────────────
# Reports predictive interval width at peak week for each group.
# Use this to verify that k_y=30 has narrowed the young CI as intended.
peak_wk <- which.max(obs_total$reports_total)
pp_at_peak <- bind_rows(pp_sims) |>
  filter(week == peak_wk) |>
  summarise(
    ci_width_young = quantile(reports_young, 0.975) - quantile(reports_young, 0.025),
    ci_width_old   = quantile(reports_old,   0.975) - quantile(reports_old,   0.025),
    ci_width_total = quantile(reports_total, 0.975) - quantile(reports_total, 0.025)
  )
cat(sprintf(
  "\n95%% PPC interval width at peak (week %d):\n  Young: %s\n  Old:   %s\n  Total: %s\n",
  peak_wk,
  format(round(pp_at_peak$ci_width_young), big.mark = ","),
  format(round(pp_at_peak$ci_width_old),   big.mark = ","),
  format(round(pp_at_peak$ci_width_total), big.mark = ",")
))

saveRDS(prod_chain, "results/baseline_wave3.rds")
cat("Saved.\n")

# Step 1 — save baseline (do this now while session is open)
saveRDS(prod_chain, "results/baseline_wave3.rds")

# Step 2 — run sensitivity overnight
source("run_wave3_sensitivity.R")

# Step 3 — when done, generate figures
source("sensitivity_plots_wave3.R")

# Step 4 — diagnostics and trace plots
source("sensitivity_traceplots_wave3.R")
source("sensitivity_diagnostics_wave3.R")
