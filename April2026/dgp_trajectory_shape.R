## DGP-only trajectory shape plot
## Path 1 from B+C discussion: composite (sum) vs marker (j=1) population
## mean trajectories under each MNI condition, normalized to p_t in [0,1].

library(ggplot2)
library(dplyr)
library(tidyr)

## --- Design (from function_new.R) ---
TT <- 5; J <- 5
mu_alpha <- 1; mu_beta <- 0.125
tau_un <- 1; lambda_un <- 1
t_idx <- 0:(TT - 1)
u <- t_idx / (TT - 1)
drift_t <- u^2

mni_grid <- list(
  list(loc = "none",      sz = "none",   n = 0, d_tau = 0,    d_lam = 0),
  list(loc = "intercept", sz = "weak",   n = 1, d_tau = 0.25, d_lam = 0),
  list(loc = "intercept", sz = "strong", n = 3, d_tau = 0.50, d_lam = 0),
  list(loc = "loading",   sz = "weak",   n = 1, d_tau = 0,    d_lam = 0.15),
  list(loc = "loading",   sz = "strong", n = 3, d_tau = 0,    d_lam = 0.30)
)

## --- Per-condition trajectory ---
traj <- function(spec) {
  tau_mat <- matrix(tau_un,   nrow = J, ncol = TT)
  lam_mat <- matrix(lambda_un, nrow = J, ncol = TT)
  if (spec$n > 0) {
    noninv <- tail(seq_len(J), spec$n)
    for (t in seq_len(TT)) {
      tau_mat[noninv, t] <- tau_un + spec$d_tau * drift_t[t]
      lam_mat[noninv, t] <- lambda_un + spec$d_lam * drift_t[t]
    }
  }
  E_eta <- mu_alpha + t_idx * mu_beta
  E_St  <- colSums(tau_mat) + colSums(lam_mat) * E_eta   # composite (sum)
  E_y1  <- tau_mat[1, ] + lam_mat[1, ] * E_eta            # marker j=1

  p_norm <- function(x) (x - x[1]) / (x[TT] - x[1])
  data.frame(
    mni_location = spec$loc,
    mni_size     = spec$sz,
    t            = seq_len(TT),
    p_composite  = p_norm(E_St),
    p_marker     = p_norm(E_y1),
    p_linear     = u
  )
}

df <- bind_rows(lapply(mni_grid, traj))

df_long <- df |>
  pivot_longer(c(p_composite, p_marker, p_linear),
               names_to = "series", values_to = "p") |>
  mutate(
    series = recode(series,
                    p_composite = "composite (sum)",
                    p_marker    = "marker (j=1)",
                    p_linear    = "linear null"),
    cond = factor(paste(mni_location, mni_size, sep = " — "),
                  levels = c("none — none",
                             "intercept — weak", "intercept — strong",
                             "loading — weak",   "loading — strong"))
  )

## --- Plot ---
p <- ggplot(df_long, aes(t, p, color = series, linetype = series, shape = series)) +
  geom_line(linewidth = 0.8) +
  geom_point(size = 2.2) +
  facet_wrap(~ cond, nrow = 1) +
  scale_color_manual(values = c("composite (sum)" = "#E69F00",
                                "marker (j=1)"    = "#0072B2",
                                "linear null"     = "grey50")) +
  scale_linetype_manual(values = c("composite (sum)" = "solid",
                                   "marker (j=1)"    = "solid",
                                   "linear null"     = "dashed")) +
  scale_shape_manual(values = c("composite (sum)" = 16,
                                "marker (j=1)"    = 17,
                                "linear null"     = NA)) +
  scale_x_continuous(breaks = 1:TT) +
  labs(
    title    = "DGP trajectory shape: composite vs marker under MNI",
    subtitle = "Population p_t = (μ_t − μ_1) / (μ_T − μ_1); rho/target_loading do not affect means",
    x = "time", y = expression(p[t]),
    color = NULL, linetype = NULL, shape = NULL
  ) +
  theme_bw(base_size = 11) +
  theme(legend.position = "bottom",
        panel.grid.minor = element_blank())

ggsave("/Users/hazel/Desktop/dgp_trajectory_shape.png",
       p, width = 13, height = 4, dpi = 150)

## --- Companion: deviation from linear null (magnified) ---
df_dev <- df |>
  mutate(dev_composite = p_composite - p_linear,
         dev_marker    = p_marker    - p_linear) |>
  pivot_longer(c(dev_composite, dev_marker),
               names_to = "series", values_to = "dev") |>
  mutate(
    series = recode(series,
                    dev_composite = "composite (sum)",
                    dev_marker    = "marker (j=1)"),
    cond = factor(paste(mni_location, mni_size, sep = " — "),
                  levels = c("none — none",
                             "intercept — weak", "intercept — strong",
                             "loading — weak",   "loading — strong"))
  )

p2 <- ggplot(df_dev, aes(t, dev, color = series, shape = series)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
  geom_line(linewidth = 0.8) +
  geom_point(size = 2.2) +
  facet_wrap(~ cond, nrow = 1) +
  scale_color_manual(values = c("composite (sum)" = "#E69F00",
                                "marker (j=1)"    = "#0072B2")) +
  scale_x_continuous(breaks = 1:TT) +
  labs(
    title    = "Deviation from linear null: p_t − (t−1)/(T−1)",
    subtitle = "Magnified view; negative = composite trajectory sags below linear",
    x = "time", y = expression(p[t] - p[t]^{linear}),
    color = NULL, shape = NULL
  ) +
  theme_bw(base_size = 11) +
  theme(legend.position = "bottom",
        panel.grid.minor = element_blank())

ggsave("/Users/hazel/Desktop/dgp_trajectory_deviation.png",
       p2, width = 13, height = 4, dpi = 150)

print(df)
cat("\nSaved: ~/Desktop/dgp_trajectory_shape.png\n")
cat("Saved: ~/Desktop/dgp_trajectory_deviation.png\n")
