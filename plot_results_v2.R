library(dplyr)
library(tidyr)
library(ggplot2)

# ── 1. Load data ──────────────────────────────────────────────────────────────
files <- list.files("results", pattern = "\\.RData$", full.names = TRUE)

all_summary <- lapply(files, function(f) {
  env <- new.env()
  load(f, envir = env)
  tbl <- env$result[["summary_table"]]
  tbl$model <- env$row$model
  tbl$rho   <- env$row$rho_val
  tbl
}) |> bind_rows()

# ── 2. Reshape to long form ───────────────────────────────────────────────────
params <- c("phi11", "phi22", "phi12", "k1", "k2")
metrics <- c("rbias", "rmse", "cov", "power")

long <- all_summary |>
  mutate(
    N     = factor(N, levels = c(200, 500, 800)),
    model = factor(model, levels = c("latent", "composite"),
                   labels = c("Latent", "Composite")),
    rho   = factor(rho, levels = c(0, 0.1), labels = c("rho == 0", "rho == 0.1"))
  ) |>
  pivot_longer(
    cols = matches("^(phi11|phi22|phi12|k1|k2)\\.(rbias|rmse|cov|power)$"),
    names_to  = c("param", "metric"),
    names_sep = "\\.",
    values_to = "value"
  ) |>
  mutate(
    param  = factor(param,  levels = params),
    metric = factor(metric, levels = metrics)
  )

# ── 3. Shared theme ───────────────────────────────────────────────────────────
base_theme <- theme_bw(base_size = 11) +
  theme(
    strip.background = element_rect(fill = "grey92"),
    legend.position  = "bottom",
    panel.grid.minor = element_blank()
  )

param_labeller <- as_labeller(
  c(phi11 = "phi[11]~(var~alpha)",
    phi22 = "phi[22]~(var~beta)",
    phi12 = "phi[12]~(cov)",
    k1    = "kappa[1]~(mu[alpha])",
    k2    = "kappa[2]~(mu[beta])"),
  default = label_parsed
)

# ── 4. Plot function ──────────────────────────────────────────────────────────
make_plot <- function(metric_name, y_label, hline = NULL, y_limits = NULL) {
  d <- filter(long, metric == metric_name)

  p <- ggplot(d, aes(x = N, y = value,
                     color = model, group = model,
                     shape = model, linetype = model)) +
    geom_line(linewidth = 0.7) +
    geom_point(size = 2.5) +
    facet_grid(param ~ rho,
               labeller = labeller(param = param_labeller,
                                   rho   = label_parsed)) +
    scale_color_manual(
      name   = "Model",
      values = c(Latent = "#1b7837", Composite = "#762a83")
    ) +
    scale_shape_manual(
      name   = "Model",
      values = c(Latent = 16, Composite = 17)
    ) +
    scale_linetype_manual(
      name   = "Model",
      values = c(Latent = "solid", Composite = "dashed")
    ) +
    labs(x = "Sample size (N)", y = y_label, title = y_label) +
    base_theme

  if (!is.null(hline)) {
    p <- p + geom_hline(yintercept = hline, linetype = "dashed",
                        color = "grey40", linewidth = 0.5)
  }
  if (!is.null(y_limits)) {
    p <- p + coord_cartesian(ylim = y_limits)
  }
  p
}

# ── 5. Generate plots ─────────────────────────────────────────────────────────
p_rbias  <- make_plot("rbias",  "RelativeBias",         hline = 0)
p_rmse  <- make_plot("rmse",  "RMSE")
p_cov   <- make_plot("cov",   "95% Coverage", hline = 0.95, y_limits = c(0, 1))
p_power <- make_plot("power", "Power",        hline = 0.80, y_limits = c(0, 1))

# ── 6. Save ───────────────────────────────────────────────────────────────────
dir.create("results/plots", showWarnings = FALSE)

ggsave("results/plots/rbias.png",  p_rbias,  width = 7, height = 9, dpi = 150)
ggsave("results/plots/rmse.png",  p_rmse,  width = 7, height = 9, dpi = 150)
ggsave("results/plots/cov.png",   p_cov,   width = 7, height = 9, dpi = 150)
ggsave("results/plots/power.png", p_power, width = 7, height = 9, dpi = 150)

cat("Saved 4 plots to results/plots/\n")

# ── 7. Print to viewer (for interactive use) ──────────────────────────────────
print(p_rbias)
print(p_rmse)
print(p_cov)
print(p_power)
