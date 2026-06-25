# =============================================================================
# sensitivity_plots_wave4.R
# Saves all figures to results/wave4/
# setwd("C:/Users/Usuario/Desktop/CBS - Internship/Analysis/DATA")
# source("sensitivity_plots_wave4.R")
# =============================================================================

library(tidyverse)

if (!dir.exists("results/wave4")) dir.create("results/wave4", recursive=TRUE)

metrics_df <- readRDS("results/wave4_metrics.rds") |>
  mutate(degradation = tolower(as.character(degradation)))
summary_df <- readRDS("results/wave4_summary.rds") |>
  mutate(degradation = tolower(as.character(degradation)))

col_sys <- "#154273"; col_sto <- "#00a1d5"; col_agg <- "#e17000"
param_levels <- c("Beta1","Beta2","Beta3","rho_y","rho_o")
param_labeller <- as_labeller(c(
  Beta1="beta[1]", Beta2="beta[2]", Beta3="beta[3]",
  rho_y="rho[y]",  rho_o="rho[o]"), label_parsed)

save_fig <- function(p, name, w=12, h=8) {
  ggsave(file.path("results/wave4", paste0(name,".pdf")), p, width=w, height=h)
  ggsave(file.path("results/wave4", paste0(name,".png")), p, width=w, height=h, dpi=150)
  cat(name, "saved\n")
}

baseline_vals <- summary_df |>
  filter(degradation=="baseline") |>
  select(parameter, base_mean=post_mean, base_lo=ci_lo, base_hi=ci_hi)

# Fig 1
fig1 <- summary_df |>
  filter(degradation!="baseline") |>
  left_join(baseline_vals, by="parameter") |>
  mutate(scenario=case_when(
    degradation=="systematic" ~ paste0("Sys\nf=",level),
    degradation=="stochastic" ~ paste0("Sto\nf=",level),
    degradation=="aggregate"  ~ paste0("Agg\n",level,"wk"),
    TRUE ~ paste0(degradation,"\n",level)),
    parameter=factor(parameter, levels=param_levels)) |>
  ggplot(aes(x=scenario, y=post_mean, colour=degradation, ymin=ci_lo, ymax=ci_hi)) +
  geom_hline(aes(yintercept=base_mean), linetype="dashed", colour="grey40", linewidth=0.6) +
  geom_hline(aes(yintercept=base_lo),   linetype="dotted", colour="grey60", linewidth=0.4) +
  geom_hline(aes(yintercept=base_hi),   linetype="dotted", colour="grey60", linewidth=0.4) +
  geom_errorbar(width=0.25, linewidth=0.7) + geom_point(size=3) +
  facet_wrap(~parameter, scales="free_y", nrow=2, labeller=param_labeller) +
  scale_colour_manual(values=c("systematic"=col_sys,"stochastic"=col_sto,"aggregate"=col_agg), name="Type") +
  labs(title="Wave 4: Posterior means and 95% CIs under data degradation",
       subtitle="Dashed=baseline mean | Dotted=baseline 95% CI", x="Scenario", y="Value") +
  theme_bw(base_size=11) + theme(legend.position="bottom", panel.grid.minor=element_blank())
save_fig(fig1, "fig1_posterior_means")

# Fig 2
fig2_data <- metrics_df |>
  mutate(scenario=case_when(
    degradation=="systematic" ~ paste0("Sys f=",level),
    degradation=="stochastic" ~ paste0("Sto f=",level),
    degradation=="aggregate"  ~ paste0("Agg ",level,"wk"),
    TRUE ~ paste0(degradation," ",level)),
    parameter=factor(parameter, levels=param_levels))

fig2 <- fig2_data |>
  ggplot(aes(x=scenario, y=ci_width_ratio, fill=degradation)) +
  geom_hline(yintercept=1, linetype="dashed", colour="grey40", linewidth=0.6) +
  geom_col(alpha=0.7, width=0.6) +
  facet_wrap(~parameter, nrow=2, labeller=param_labeller) +
  scale_fill_manual(values=c("systematic"=col_sys,"stochastic"=col_sto,"aggregate"=col_agg), name=NULL) +
  labs(title="Wave 4: CI width ratio vs baseline", subtitle="Ratio > 1 = wider than baseline",
       x="Scenario", y="CI width ratio") +
  theme_bw(base_size=11) +
  theme(legend.position="none", axis.text.x=element_text(angle=45, hjust=1, size=8),
        panel.grid.minor=element_blank())
save_fig(fig2, "fig2_ci_ratio")

# Fig 3
fig3 <- fig2_data |>
  ggplot(aes(x=scenario, y=mean_shift_pct, fill=degradation)) +
  geom_hline(yintercept=0, linetype="dashed", colour="grey40", linewidth=0.6) +
  geom_col(alpha=0.7, width=0.6) +
  facet_wrap(~parameter, nrow=2, scales="free_y", labeller=param_labeller) +
  scale_fill_manual(values=c("systematic"=col_sys,"stochastic"=col_sto,"aggregate"=col_agg), name=NULL) +
  labs(title="Wave 4: Posterior mean shift vs baseline",
       subtitle="Positive=estimate higher than baseline", x="Scenario", y="Mean shift (%)") +
  theme_bw(base_size=11) +
  theme(legend.position="none", axis.text.x=element_text(angle=45, hjust=1, size=8),
        panel.grid.minor=element_blank())
save_fig(fig3, "fig3_mean_shift")

# Fig 4
sys_data <- summary_df |>
  filter(degradation=="systematic") |>
  left_join(baseline_vals, by="parameter") |>
  mutate(f=as.numeric(level), observed_ratio=post_mean/base_mean,
         param_type=ifelse(grepl("rho",parameter),"Detection rate","Transmission rate"),
         predicted=ifelse(grepl("rho",parameter),(1-f),1.0),
         parameter=factor(parameter, levels=param_levels))
theory_lines <- sys_data |> distinct(f, param_type, predicted)
fig4 <- ggplot() +
  geom_hline(yintercept=1, colour="grey80", linewidth=0.4) +
  geom_line(data=theory_lines, aes(x=f,y=predicted,linetype=param_type), colour="grey30", linewidth=0.8) +
  geom_line(data=sys_data, aes(x=f,y=observed_ratio,colour=parameter,group=parameter), linewidth=0.9) +
  geom_point(data=sys_data, aes(x=f,y=observed_ratio,colour=parameter,shape=parameter), size=3) +
  scale_x_continuous(breaks=c(0.10,0.30,0.50), labels=c("10%","30%","50%")) +
  scale_colour_brewer(palette="Set1", name="Parameter") +
  scale_shape_manual(values=c(16,17,15,3,4), name="Parameter") +
  scale_linetype_manual(values=c("Detection rate"="dashed","Transmission rate"="dotted"),
                        name="Theoretical\nprediction") +
  labs(title="Wave 4: Theory vs observed — systematic reduction",
       subtitle="Dashed grey=predicted rho: (1-f) | Dotted grey=predicted beta: no change",
       x="Fraction removed (f)", y="Ratio: degraded / baseline") +
  theme_bw(base_size=11) + theme(legend.position="right", panel.grid.minor=element_blank())
save_fig(fig4, "fig4_theory_vs_observed", w=10, h=6)

# Fig 5
fig5_data <- metrics_df |>
  mutate(scenario=case_when(
    degradation=="systematic" ~ paste0("Sys\nf=",level),
    degradation=="stochastic" ~ paste0("Sto\nf=",level),
    degradation=="aggregate"  ~ paste0("Agg\n",level,"wk"),
    TRUE ~ paste0(degradation,"\n",level)),
    scenario=factor(scenario, levels=unique(scenario)),
    coverage_label=ifelse(coverage,"\u2713","\u2717"),
    parameter=factor(parameter, levels=param_levels))
fig5 <- ggplot(fig5_data, aes(x=scenario, y=parameter, fill=coverage)) +
  geom_tile(colour="white", linewidth=0.8) +
  geom_text(aes(label=coverage_label), size=5, colour="white", fontface="bold") +
  scale_fill_manual(values=c("TRUE"="#154273","FALSE"="#c0392b"),
                    labels=c("TRUE"="Covered","FALSE"="Not covered"), name="Baseline CI") +
  scale_y_discrete(labels=c(Beta1=expression(beta[1]),Beta2=expression(beta[2]),
                             Beta3=expression(beta[3]),rho_y=expression(rho[y]),
                             rho_o=expression(rho[o]))) +
  labs(title="Wave 4: Coverage heatmap",
       subtitle="\u2713=robust | \u2717=shift exceeds baseline uncertainty",
       x="Scenario", y="Parameter") +
  theme_bw(base_size=11) +
  theme(panel.grid=element_blank(), legend.position="right", axis.text.x=element_text(size=9))
save_fig(fig5, "fig5_coverage", w=11, h=5)

# Fig 6
fig6_data <- metrics_df |>
  filter(degradation %in% c("systematic","stochastic")) |>
  mutate(f=as.numeric(level),
         degradation=factor(degradation, levels=c("systematic","stochastic"),
                            labels=c("Systematic","Stochastic")),
         parameter=factor(parameter, levels=param_levels))
fig6 <- ggplot(fig6_data, aes(x=f,y=ci_width_ratio,colour=degradation,
                               linetype=degradation,group=degradation)) +
  geom_hline(yintercept=1, linetype="dashed", colour="grey70", linewidth=0.4) +
  geom_line(linewidth=1.0) + geom_point(size=3) +
  facet_wrap(~parameter, nrow=2, labeller=param_labeller) +
  scale_colour_manual(values=c("Systematic"=col_sys,"Stochastic"=col_sto), name=NULL) +
  scale_linetype_manual(values=c("Systematic"="solid","Stochastic"="dashed"), name=NULL) +
  scale_x_continuous(breaks=c(0.10,0.30,0.50), labels=c("10%","30%","50%")) +
  labs(title="Wave 4: Systematic vs stochastic — CI width ratio",
       subtitle="Gap=extra uncertainty from binomial noise alone",
       x="Fraction removed (f)", y="CI width ratio") +
  theme_bw(base_size=11) +
  theme(legend.position="bottom", panel.grid.minor=element_blank())
save_fig(fig6, "fig6_sys_vs_sto")

cat("\nAll Wave 4 figures saved to results/wave4/\n")
