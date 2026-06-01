#--- Setup ---------------------

suppressPackageStartupMessages({
  library(tidyverse)
  library(conflicted)
})
conflicts_prefer(dplyr::select)
conflicts_prefer(dplyr::filter)


#--- Load results --------------
# NEW (may): derive the design labels the plot function needs directly from the
# raw numeric condition columns (target_loading / rho / fit_cu), so it works
# regardless of whether the saved `row` already carried label columns.

load_results <- function(dir = "study1_results_Fit_CU") {
  files <- list.files(dir, pattern = "\\.RData$", full.names = TRUE)
  rows  <- lapply(files, function(f) {
    e <- new.env(); load(f, envir = e)
    keep <- setdiff(names(e$row), names(e$result$summary_table))
    bind_cols(e$row[, keep, drop = FALSE], e$result$summary_table)
  })
  out <- bind_rows(rows) %>%
    mutate(admissible.m = 1 - err.perc)

  # x-axis alias: summary_table stores N, plots default to N_size
  if ("N" %in% names(out) && !"N_size" %in% names(out))
    out <- out %>% mutate(N_size = N)

  # R^2 / measurement-quality label from target_loading (.6708 -> vgood, .4 -> poor)
  if ("target_loading" %in% names(out))
    out <- out %>%
      mutate(loading_quality = factor(
        if_else(target_loading > 0.5, "vgood", "poor"),
        levels = c("poor", "vgood")))

  # CU on/off in the DATA
  if ("rho" %in% names(out))
    out <- out %>%
      mutate(rho_lab = factor(
        recode(as.character(rho), "0" = "Basic", "0.15" = "Cor.Unique"),
        levels = c("Basic", "Cor.Unique")))

  # CU in the MODEL (fit_cu) and the data x model spec cross
  if ("fit_cu" %in% names(out))
    out <- out %>%
      mutate(cu_model = factor(if_else(fit_cu, "modeled", "ignored"),
                               levels = c("ignored", "modeled")))

  if (all(c("rho", "fit_cu") %in% names(out)))
    out <- out %>%
      mutate(spec = factor(case_when(
        rho == 0 & !fit_cu ~ "Correct (no CU)",
        rho != 0 & !fit_cu ~ "CU ignored",
        rho != 0 &  fit_cu ~ "CU modeled",
        TRUE               ~ "CU overfit"),
        levels = c("Correct (no CU)", "CU ignored", "CU modeled", "CU overfit")))

  # Drop Study 2 residual cells if column exists
  if ("mni_location" %in% names(out))
    out <- out %>% filter(mni_location != "residual")
  out
}

results <- load_results()


#--- Config --------------------
# Okabe-Ito, colorblind-friendly
PALETTE      <- c(sum = "#E69F00", EC_latent = "#0072B2")
MODEL_LABELS <- c(sum = "Sum 1-LGCM", EC_latent = "EC 2-LGCM")

METRIC_LABELS <- c(
  rbias      = "Relative bias",
  rmse       = "Relative RMSE",
  power      = "Power",
  cov        = "Coverage",
  admissible = "Admissible rate",
  fit        = "Fit index"
)

# Reference lines (NULL = none)
REF_LINES <- list(
  rbias      = c(-0.10, 0, 0.10),
  power      = 0.80,
  cov        = 0.90,
  admissible = 0.90,
  cfi        = c(0.90, 0.95),
  tli        = c(0.90, 0.95),
  rmsea      = c(0.06, 0.08),
  srmr       = 0.08
)


# metric -> column suffix in df (so "rmse" pulls relative RMSE columns)
.METRIC_SUFFIX <- c(rbias = "rbias", rmse = "rel.rmse",
                    power = "power", cov = "cov")

# k1/k2 are named kappa1/kappa2 only in rel.rmse columns
.PARAM_ALIAS <- list("rel.rmse" = c(k1 = "kappa1", k2 = "kappa2"))


#--- Column picker -------------

.pick_cols <- function(df, metric, params, indices) {
  if (metric == "admissible")
    return(tibble(col = "admissible.m", panel = "admissible"))

  if (metric == "fit") {
    if (is.null(indices)) indices <- c("cfi", "tli", "rmsea", "srmr")
    return(tibble(col = paste0(indices, ".m"), panel = indices) %>%
             filter(col %in% names(df)))
  }

  if (metric == "chisq_reject")
    return(tibble(col = "chisq.reject", panel = "chisq rejection"))

  suffix <- .METRIC_SUFFIX[metric] %||% metric
  alias  <- .PARAM_ALIAS[[suffix]]

  if (is.null(params)) {
    # pull all params available for this metric
    pat    <- paste0("\\.", suffix, "$")
    cols   <- grep(pat, names(df), value = TRUE)
    panels <- sub(pat, "", cols)
    if (!is.null(alias)) {
      rev_alias <- setNames(names(alias), alias)
      panels <- ifelse(panels %in% names(rev_alias), rev_alias[panels], panels)
    }
    return(tibble(col = cols, panel = panels))
  }

  # explicit params: resolve aliases when needed (e.g. k1 -> kappa1)
  resolved <- if (!is.null(alias))
    ifelse(params %in% names(alias), alias[params], params) else params
  tibble(col = paste0(resolved, ".", suffix), panel = params) %>%
    filter(col %in% names(df))
}

# Pretty facet labels
# NEW (may): added spec / cu_model entries; unlisted vars fall back to label_value
.facet_labeller <- labeller(
  J               = function(x) paste0("J = ", x),
  TT              = function(x) paste0("T = ", x),
  N_size          = function(x) paste0("N = ", x),
  rho             = function(x) paste0("rho = ", x),
  sd.lambda       = function(x) c("0" = "tau-equiv", "0.05" = "congeneric")[as.character(x)],
  loading_quality = function(x) c(poor = "Low R^2", vgood = "High R^2")[x],
  cu_model        = function(x) c(ignored = "CU ignored", modeled = "CU modeled")[x],
  rho_lab         = label_value,    # already "Basic" / "Cor.Unique"
  spec            = label_value,    # already "Correct (no CU)" / "CU ignored" / ...
  panel           = label_value
)

#--- Single plot function ------
plot_metric <- function(df, metric,
                        params  = NULL,           # growth params (rbias/rmse/power/cov)
                        indices = NULL,           # fit indices (metric = "fit")
                        models  = c("sum", "EC_latent"),
                        x       = "N_size",
                        x_lab   = "Sample size",
                        y_range = NULL,
                        facets  = c("sd.lambda"),
                        filters = list(),         # e.g. list(rho = 0, J = 5)
                        title   = NULL) {

  # 1. restrict to chosen models
  df <- df %>% filter(model %in% models) %>%
    mutate(model = factor(model, levels = models))

  # 2. ad-hoc filters
  for (nm in names(filters)) df <- df %>% filter(.data[[nm]] %in% filters[[nm]])
  if (!nrow(df)) stop("filter produced empty data")

  # 3. resolve target columns and reshape
  cols <- .pick_cols(df, metric, params, indices)
  if (!nrow(cols)) stop(sprintf("no columns matched metric '%s'", metric))

  grp       <- intersect(c(x, "model", "loading_quality", facets), names(df))
  panel_map <- setNames(cols$panel, cols$col)

  long <- df %>%
    select(all_of(grp), all_of(cols$col)) %>%
    pivot_longer(all_of(cols$col), names_to = "panel", values_to = "value") %>%
    mutate(panel = factor(panel_map[panel], levels = cols$panel)) %>%
    group_by(across(all_of(c(grp, "panel")))) %>%
    summarise(value = mean(value, na.rm = TRUE), .groups = "drop")

  # 4. reference lines
  ref <- if (metric == "fit") {
    bind_rows(lapply(levels(long$panel), function(p) {
      v <- REF_LINES[[p]]
      if (is.null(v)) NULL else tibble(panel = p, y = v)
    }))
  } else if (!is.null(REF_LINES[[metric]])) {
    tibble(y = REF_LINES[[metric]])
  } else NULL

  # 5. build plot
  p <- ggplot(long, aes(.data[[x]], value,
                        color = model, shape = model,
                        linetype = loading_quality,
                        group = interaction(model, loading_quality))) +
    geom_line(linewidth = 0.6) + geom_point(size = 2) +
    scale_color_manual(values = PALETTE[models],
                       labels = MODEL_LABELS[models], name = "Model") +
    scale_shape_manual(values = c(16, 17)[seq_along(models)],
                       labels = MODEL_LABELS[models], name = "Model") +
    scale_linetype_manual(values = c(poor = "dashed", vgood = "solid"),
                          labels = c(poor = "Low R^2", vgood = "High R^2"),
                          name   = expression(R^2)) +
    labs(x = x_lab, y = METRIC_LABELS[metric] %||% metric,
         title = title %||% METRIC_LABELS[metric] %||% metric) +
    theme_bw(base_size = 10, base_family = "Times New Roman") +
    theme(panel.grid.minor = element_blank(),
          strip.background = element_rect(fill = "grey92"),
          axis.text.x      = element_text(angle = 45, hjust = 1))

  if (!is.null(ref) && nrow(ref)) {
    p <- if ("panel" %in% names(ref))
      p + geom_hline(data = ref, aes(yintercept = y),
                     linetype = 2, color = "grey60", inherit.aes = FALSE)
    else
      p + geom_hline(yintercept = ref$y, linetype = 2, color = "grey60")
  }

  # shared vs free y-axis
  free_scale <- if (is.null(y_range)) "free_y" else "fixed"

  rhs <- if (length(facets)) paste(facets, collapse = " + ") else "."
  p <- p + facet_grid(as.formula(
    if (nrow(cols) > 1) paste("panel ~", rhs) else paste(". ~", rhs)
  ), scales = free_scale, labeller = .facet_labeller)

  # fix y range to the data max (include reference lines so they stay visible)
  if (!is.null(y_range)) {
    yl <- if (identical(y_range, "shared")) {
      vals <- long$value
      if (!is.null(ref) && nrow(ref)) vals <- c(vals, ref$y)
      range(vals, na.rm = TRUE)
    } else y_range
    p <- p + coord_cartesian(ylim = yl)   # clips view only, never drops points
  }

  if (is.numeric(long[[x]]))
    p <- p + scale_x_continuous(breaks = sort(unique(long[[x]])))
  p
}


#====================================================================
#  REPORTING COMPOSITION  (Study 1, may design)
#--------------------------------------------------------------------
#  Encoding held constant across every figure:
#    x        = N (200/500/800)
#    color/shape = model   (sum vs EC 2-LGCM)
#    linetype = R^2        (solid = Low, dashed = High)
#  What changes per figure = the facet grid (rows = growth params,
#  columns = the design contrast that answers a given RQ).
#
#  Design dimensions available for facets/filters:
#    spec   : "Correct (no CU)" / "CU ignored" / "CU modeled" / "CU overfit"
#             (= rho x fit_cu).  All 4 exist only at T=5; T=3 has the first two.
#    rho_lab: "Basic" / "Cor.Unique"   (CU in the DATA)
#    cu_model: "ignored" / "modeled"   (CU in the MODEL, fit_cu)
#    sd.lambda (tau-equiv/congeneric), J (3/5), TT (3/5)
#====================================================================


#--- FIG 0. Convergence first (gate everything else) ----------------
# Worst-case J=3 & T=3 cell is the one to watch. Read admissible rate
# before trusting any bias/coverage panel.
# plot_metric(results, "admissible",
#             facets  = c("J", "TT",  "spec"),
#             filters = list(sd.lambda = 0),
#             title   = "Admissible solution rate")

plot_metric(results, "admissible",
            facets  = c("J", "spec"),
            filters = list(sd.lambda = 0, TT=3),
            title   = "T = 3, Admissible solution rate")

plot_metric(results, "admissible",
            facets  = c("spec", "J"),
            filters = list(sd.lambda = 0, TT=5),
            title   = "T = 5, Admissible solution rate")


#--- FIG 1. RQ2  Which parameters break when CU is ignored? ---------
# Core robustness map at T=5: rows = growth params, cols = spec.
# Read left->right: Correct(no CU) -> CU ignored shows the bias CU induces;
# CU modeled shows whether fitting it repairs the bias.
plot_metric(results, "rbias",
            params  = c("phi11", "phi22", "phi12","k1", "k2"),
            facets  = c("spec", "J"),
            filters = list(TT = 5, J=3,sd.lambda = 0),
            y_range = "shared",
            title   = "Relative bias by specification (T = 5)")

plot_metric(results, "rmse",
            params  = c("phi11", "phi22", "phi12", "k1", "k2"),
            facets  = c("spec", "J"),
            filters = list(TT = 5,J=3, sd.lambda = 0),
            y_range = "shared",
            title   = "Relative RMSE by specification (T = 5)")

# Coverage companion (same layout) — the inferential consequence of FIG 1.
plot_metric(results, "cov",
            params  = c("phi11", "phi22", "phi12", "k1", "k2"),
            facets  = c("spec", "J"),
            filters = list( TT = 5, sd.lambda = 0),
            title   = "95% CI coverage by specification (T = 5)")


#--- FIG 2. RQ1  The 3-wave dilemma --------------------------------
# At T=3 nobody can model CU. Compare sum vs EC under CU present/absent.
# cols = rho_lab (Basic vs Cor.Unique); both models are forced no-CU here.
plot_metric(results, "rbias",
            params  = c("phi11", "phi22", "phi12", "k1", "k2"),
            facets  = c("rho_lab", "J"),
            filters = list( sd.lambda=0, TT=3),
            y_range = "shared",
            title   = "Relative bias by specification (T = 3)")

plot_metric(results, "rmse",
            params  = c("phi11", "phi22", "phi12", "k1", "k2"),
            facets  = c("rho_lab", "J"),
            filters = list( sd.lambda=0, TT=3),
            y_range = "shared",
            title   = "Relative RMSE by specification (T = 3)")


plot_metric(results, "cov",
            params  = c("phi11", "phi22", "phi12", "k1", "k2"),
            facets  = c("rho_lab", "J"),
            filters = list( sd.lambda=0, TT=3),
            y_range = "shared",
            title   = "Coverage by specification (T = 3)")


#--- FIG 3. RQ3  Does modeling CU pay off at T=5? ------------------
# Data has CU (rho = .15). cols = cu_model (ignored vs modeled), de-confounded
# from T. Gap = the value added by being able to specify CU.
plot_metric(results, "rbias",
            params  = c("phi11", "phi22", "phi12"),
            facets  = c("cu_model", "J"),
            y_range = "shared",
            filters = list(rho = 0.15, TT = 5, sd.lambda=0),
            title   = "T = 5, CU present: relative bias ignoring vs modeling CU")

# Efficiency side: does the extra CU machinery cost RMSE?
plot_metric(results, "rmse",
            params  = c("phi11", "phi22", "phi12"),
            facets  = c("cu_model", "J"),
            filters = list(rho = 0.15, TT = 5, sd.lambda=0),
            y_range = "shared",
            title   = "T = 5, CU present: relative RMSE, ignoring vs modeling")


# Efficiency side: does the extra CU machinery cost RMSE?
plot_metric(results, "cov",
            params  = c("phi11", "phi22", "phi12"),
            facets  = c("cu_model", "J"),
            filters = list(rho = 0.15, TT = 5, sd.lambda=0),
            y_range = "shared",
            title   = "T = 5, CU present: coverage, ignoring vs modeling")


#--- FIG 4. Scale-length buffer (J = 3 vs 5) -----------------------
# Does a longer scale buffer the sum score against CU bias?
# Fix the misspecified arm (CU ignored), cols = J.
# plot_metric(results, "rbias",
#             params  = c("phi11", "phi22", "phi12"),
#             facets  = c("J"),
#             filters = list(rho = 0.15, TT=5, fit_cu = FALSE, sd.lambda=0),
#             title   = "CU ignored (T = 5): does scale length buffer bias?")

plot_metric(results, "rbias",
            params  = c("phi11", "phi22", "phi12"),
            facets  = c("TT","J"),
            filters = list(rho = 0.15, fit_cu = FALSE, sd.lambda=0),
            title   = "CU ignored: does scale length buffer bias?")

#--- FIG 6. Fit-index behavior -------------------------------------
# How fit indices respond to the misspecification (mainly diagnostic for EC).
plot_metric(results, "fit",
            indices = c("cfi", "tli","rmsea", "srmr"),
            facets  = c("spec","J"),
            filters = list(TT = 5, sd.lambda = 0),
            title   = "Fit indices by specification (T = 5)")


#--- FIG 7. Power (slope mean & intercept-slope cov) ---------------
plot_metric(results, "power",
            params  = c("phi11", "phi22", "phi12", "k1", "k2"),
            facets  = c("spec", "J"),
            filters = list(TT = 5, sd.lambda = 0),
            title   = "Power by CU specification (T = 5)")


#--- Saving helper -------------------------------------------------
# save_fig(p, "fig1_rbias_spec")  ->  figs/fig1_rbias_spec.png
save_fig <- function(p, name, w = 9, h = 6, dir = "figs") {
  if (!dir.exists(dir)) dir.create(dir, recursive = TRUE)
  ggsave(file.path(dir, paste0(name, ".png")), p,
         width = w, height = h, dpi = 320)
}
# Example:
# p1 <- plot_metric(results, "rbias", params = c("k1","k2","phi11","phi22","phi12"),
#                   facets = "spec", filters = list(J = 5, TT = 5, sd.lambda = 0))
# save_fig(p1, "fig1_rbias_spec")
