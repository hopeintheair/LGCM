library(dplyr)
library(tidyr)
library(ggplot2)

# ── 1. Load data
files <- list.files("results2", pattern = "\\.RData$", full.names = TRUE)

all_summary <- lapply(files, function(f) {
  env <- new.env()
  load(f, envir = env)
  tbl <- env$result[["summary_table"]]
  tbl$model <- env$row$model
  tbl$rho   <- env$row$rho_val
  tbl
}) |> bind_rows()

# ── 2. Reshape to long form 
fit_params <- c("pvalue.chi.m", "chisq.reject", "cfi.m", "tli.m", "rmsea.m", "srmr.m")

long <- all_summary |>
  mutate(
    N     = factor(N, levels = c(200, 500, 800)),
    model = factor(model, levels = c("latent", "composite"),
                   labels = c("Latent", "Composite")),
    rho   = factor(rho, levels = c(0, 0.1), labels = c("rho == 0", "rho == 0.1"))
  ) |>
  pivot_longer(
    cols      = all_of(fit_params), 
    names_to  = "fit_index",
    values_to = "value"
  ) |>
  mutate(
    fit_index = factor(fit_index, levels = fit_params)
  )

# ── 3. Shared theme 
base_theme <- theme_bw(base_size = 11) +
  theme(
    strip.background = element_rect(fill = "grey92"),
    legend.position  = "bottom",
    panel.grid.minor = element_blank()
  )

# ── 4. Plot function 
make_plot <- function(index_name, y_label, hline = NULL, y_limits = NULL) {
  d <- filter(long, fit_index == index_name)
  
  p <- ggplot(d, aes(x = N, y = value,
                     color = model, group = model,
                     shape = model, linetype = model)) +
    geom_line(linewidth = 0.7) +
    geom_point(size = 2.5) +
    
    facet_wrap(~ rho, labeller = labeller(rho = label_parsed)) +
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
  
  # Cut-off
  if (!is.null(hline)) {
    p <- p + geom_hline(yintercept = hline, linetype = "dashed",
                        color = "grey40", linewidth = 0.5)
  }

  if (!is.null(y_limits)) {
    p <- p + coord_cartesian(ylim = y_limits)
  }
  p
}

# ── 5. Generate plots 
p_pval  <- make_plot("pvalue.chi.m", "Chi-square p-value")
p_chisq <- make_plot("chisq.reject", "Chi-square Rejection Rate", hline = 0.05, y_limits = c(0, 1))
p_cfi   <- make_plot("cfi.m",        "CFI",   hline = 0.95, y_limits = c(0.8, 1))
p_tli   <- make_plot("tli.m",        "TLI",   hline = 0.95, y_limits = c(0.8, 1))
p_rmsea <- make_plot("rmsea.m",      "RMSEA", hline = 0.05, y_limits = c(0, 0.15))
p_srmr  <- make_plot("srmr.m",       "SRMR",  hline = 0.08, y_limits = c(0, 0.15))

# ── 6. Save
dir.create("results2/plots_fit", showWarnings = FALSE)

ggsave("results2/plots_fit/pvalue.png",  p_pval,  width = 7, height = 4.5, dpi = 150)
ggsave("results2/plots_fit/chisq.png",   p_chisq, width = 7, height = 4.5, dpi = 150)
ggsave("results2/plots_fit/cfi.png",     p_cfi,   width = 7, height = 4.5, dpi = 150)
ggsave("results2/plots_fit/tli.png",     p_tli,   width = 7, height = 4.5, dpi = 150)
ggsave("results2/plots_fit/rmsea.png",   p_rmsea, width = 7, height = 4.5, dpi = 150)
ggsave("results2/plots_fit/srmr.png",    p_srmr,  width = 7, height = 4.5, dpi = 150)

# ── 7. Print to viewer (for interactive use) 
print(p_pval)
print(p_chisq)
print(p_cfi)
print(p_tli)
print(p_rmsea)
print(p_srmr)