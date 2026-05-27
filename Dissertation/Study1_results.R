## =====================================================================
## Study 1 - Results exhibits for the paper
##   Table 0  specification map
##   Table 1  plim vs analytic target (estimand / consistency)
##   Figure 1 finite-sample -> asymptotic bias convergence (N incl. Inf)
##   Table 2  finite-N performance, stratified by TT (Panel A/B)
##   Table 3  population misfit (asymptotic fit-index paradox)
##
## Inputs : study1_T35.csv (Monte Carlo), study1_plim_T35.csv (plim)
## Outputs: tableN_*.csv + figure1_bias_convergence.{png,pdf}
##
## Aggregation rule: average over J and loading_type (the latter is
## provably irrelevant: sum & EC depend on loadings only via s_t = sum,
## held equal across tau_eq/congeneric). Reliability (vgood/poor) and
## TT x rho are kept because they interact. Full per-cell results -> supplement.
## =====================================================================

suppressMessages({ library(dplyr); library(tidyr); library(ggplot2); library(scales) })

mc <- read.csv("study1_T35.csv")
pl <- read.csv("study1_plim_T35.csv")

relab <- function(d) within(d, {
  reliability <- factor(loading_quality, levels = c("vgood","poor"),
                        labels = c("Reliable (R²=.45)","Poor (R²=.16)"))
})
mc <- relab(mc); pl <- relab(pl)

## ---------------------------------------------------------------------
## Table 0 - specification map
## ---------------------------------------------------------------------
tab0 <- data.frame(
  `Data CU`     = c("No (rho=0)", "No (rho=0)", "Yes (rho=.15)", "Yes (rho=.15)"),
  `Model CU`    = c("-", "-", "Omitted", "Modeled"),
  `Design cells`= c("rho=0, TT=3", "rho=0, TT=5",
                    "rho=.15, TT=3", "rho=.15, TT=5"),
  `Identified?` = c("yes", "yes",
                    "CU not identified in sum 1-LGCM (df = -1)", "yes (df = 6)"),
  `Specification` = c("correct", "correct", "MISSPECIFIED", "correct"),
  check.names = FALSE)
write.csv(tab0, "table0_spec_map.csv", row.names = FALSE)

## ---------------------------------------------------------------------
## Table 1 - plim vs Phi_ec  (estimand validation / consistency)
## ---------------------------------------------------------------------
## (1a) correctly specified: report the worst-case |arbias| across all cells
t1_correct <- pl %>%
  filter(spec_status != "MISSPEC (CU omitted)") %>%
  summarise(across(c(phi1.arbias, phi2.arbias, phi12.arbias, cor12.arbias),
                   ~ max(abs(.x)))) %>%
  mutate(block = "Correctly specified (all cells)", .before = 1)

## (1b) misspecified: detail by reliability x model
t1_misspec <- pl %>%
  filter(spec_status == "MISSPEC (CU omitted)") %>%
  group_by(model, reliability) %>%
  summarise(across(c(phi1.arbias, phi2.arbias, phi12.arbias, cor12.arbias),
                   ~ mean(.x)), .groups = "drop") %>%
  mutate(block = "MISSPEC: TT=3, rho=.15", .before = 1)

write.csv(t1_correct, "table1a_consistency.csv", row.names = FALSE)
write.csv(t1_misspec, "table1b_asymptotic_bias.csv", row.names = FALSE)

cat("\n================ TABLE 1a  (consistency: max |asymptotic rel.bias|) ================\n")
print(t1_correct, digits = 2)
cat("\n================ TABLE 1b  (asymptotic rel.bias, MISSPEC cells) ================\n")
print(as.data.frame(t1_misspec), digits = 3)

## ---------------------------------------------------------------------
## Figure 1 - bias convergence: N = 200, 500, 800, Inf(plim)
## ---------------------------------------------------------------------
## finite-N relative bias (avg over J, loading_type)
mc_fig <- mc %>%
  group_by(model, TT, rho, reliability, N_size) %>%
  summarise(phi2 = mean(phi22.rbias), cor12 = mean(cor12.rbias), .groups = "drop") %>%
  rename(N = N_size)

## plim (N = Inf)
pl_fig <- pl %>%
  group_by(model, TT, rho, reliability) %>%
  summarise(phi2 = mean(phi2.arbias), cor12 = mean(cor12.arbias), .groups = "drop") %>%
  mutate(N = Inf)

fig <- bind_rows(mc_fig, pl_fig) %>%
  pivot_longer(c(phi2, cor12), names_to = "param", values_to = "rbias") %>%
  mutate(
    Nlab  = factor(ifelse(is.infinite(N), "∞", as.character(N)),
                   levels = c("200","500","800","∞")),
    xpos  = as.numeric(Nlab),
    param = factor(param, levels = c("phi2","cor12"),
                   labels = c("phi[2]~(slope~variance)", "rho[IS]~(intercept-slope~cor.)")),
    cond  = factor(paste0("TT=", TT, ", rho=", rho),
                   levels = c("TT=3, rho=0","TT=3, rho=0.15","TT=5, rho=0","TT=5, rho=0.15")),
    model = factor(model, levels = c("EC_latent","sum"),
                   labels = c("EC 2-LGCM","Sum-score 1-LGCM")))

p <- ggplot(fig, aes(xpos, rbias, colour = model, linetype = reliability,
                     group = interaction(model, reliability))) +
  geom_hline(yintercept = 0, linewidth = .3, colour = "grey50") +
  geom_line(linewidth = .6) +
  geom_point(size = 1.7) +
  geom_vline(xintercept = 3.5, linetype = "dotted", colour = "grey60") +
  facet_grid(param ~ cond, scales = "free_y", labeller = labeller(param = label_parsed)) +
  scale_x_continuous(breaks = 1:4, labels = c("200","500","800","∞"),
                     expand = expansion(mult = .08)) +
  scale_colour_manual(values = c("EC 2-LGCM" = "#1b6ca8", "Sum-score 1-LGCM" = "#c0392b")) +
  labs(x = "Sample size N  (∞ = analytic plim)", y = "Relative bias",
       colour = NULL, linetype = NULL,
       title = "Finite-sample bias vanishes only where the plim is unbiased",
       subtitle = "Dotted line separates Monte Carlo (N≤800) from the analytic probability limit") +
  theme_bw(base_size = 11) +
  theme(legend.position = "top", panel.grid.minor = element_blank(),
        strip.background = element_rect(fill = "grey92"))

ggsave("figure1_bias_convergence.png", p, width = 10, height = 6, dpi = 200)
ggsave("figure1_bias_convergence.pdf", p, width = 10, height = 6)

## ---------------------------------------------------------------------
## Table 2 - finite-N performance, stratified by TT
## ---------------------------------------------------------------------
t2 <- mc %>%
  group_by(TT, model, rho, reliability, N_size) %>%
  summarise(
    phi2.rbias  = mean(phi22.rbias),
    phi2.relrmse= mean(phi22.rel.rmse),
    cor12.rbias = mean(cor12.rbias),
    phi2.cov    = mean(phi22.cov),
    cor12.cov   = mean(cor12.cov),
    pct.admiss  = mean(admissible.m),
    .groups = "drop") %>%
  arrange(TT, model, rho, reliability, N_size)

t2_A <- filter(t2, TT == 3); t2_B <- filter(t2, TT == 5)
write.csv(t2_A, "table2A_finiteN_TT3.csv", row.names = FALSE)
write.csv(t2_B, "table2B_finiteN_TT5.csv", row.names = FALSE)

cat("\n================ TABLE 2  Panel B (TT=5): EC near-unbiased ================\n")
print(as.data.frame(t2_B %>% filter(rho == 0.15) %>%
        select(model, reliability, N_size, phi2.rbias, cor12.rbias, phi2.cov, cor12.cov, pct.admiss)),
      digits = 3)
cat("\n================ TABLE 2  Panel A (TT=3): admissibility collapse ================\n")
print(as.data.frame(t2_A %>% filter(rho == 0.15) %>%
        select(model, reliability, N_size, phi2.rbias, cor12.rbias, pct.admiss)),
      digits = 3)

## ---------------------------------------------------------------------
## Table 3 - asymptotic fit-index paradox (from plim CSV)
## ---------------------------------------------------------------------
t3 <- pl %>%
  mutate(cell = ifelse(spec_status == "MISSPEC (CU omitted)",
                       "MISSPEC: TT=3, rho=.15", "Correct: TT=5, rho=.15")) %>%
  filter((spec_status == "MISSPEC (CU omitted)") |
         (TT == 5 & rho == 0.15)) %>%
  group_by(cell, model, reliability) %>%
  summarise(rmsea.pop = mean(rmsea.pop), srmr.pop = mean(srmr.pop),
            cfi.pop = mean(cfi.pop), df.pop = mean(df.pop),
            phi2.arbias = mean(phi2.arbias), cor12.arbias = mean(cor12.arbias),
            .groups = "drop") %>%
  arrange(desc(cell), model, reliability)
write.csv(t3, "table3_fit_paradox.csv", row.names = FALSE)

cat("\n================ TABLE 3  (good population fit, large bias) ================\n")
print(as.data.frame(t3), digits = 3)

cat("\n\nWrote: table0_spec_map.csv, table1a_consistency.csv, table1b_asymptotic_bias.csv,",
    "\n       table2A_finiteN_TT3.csv, table2B_finiteN_TT5.csv, table3_fit_paradox.csv,",
    "\n       figure1_bias_convergence.{png,pdf}\n")
