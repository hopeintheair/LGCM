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

load_results <- function(dir = "study1_results_draft") {
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
  fit        = "Fit index",
  m          = "Mean estimate"
)
# Reference lines (NULL = none)
REF_LINES <- list(
  rbias      = c(-0.10, 0, 0.10),
  power      = 0.80,
  cov        = 0.95,
  admissible = 0.90,
  cfi        = c(0.90, 0.95),
  tli        = c(0.90, 0.95),
  rmsea      = c(0.06, 0.08),
  srmr       = 0.08
)


# metric -> column suffix in df (so "rmse" pulls relative RMSE columns)
.METRIC_SUFFIX <- c(rbias = "rbias", rmse = "rel.rmse",
                    power = "power", cov = "cov", m = "m")

# k1/k2 are named kappa1/kappa2 only in rel.rmse columns
.PARAM_ALIAS <- list("rel.rmse" = c(k1 = "kappa1", k2 = "kappa2"),
                     "m"        = c(k1 = "kappa1", k2 = "kappa2"))


#--- Column picker -------------

.pick_cols <- function(df, metric, params, indices) {
  if (metric == "admissible")
    return(tibble(col = "admissible.m", panel = "admissible"))
  
  if (metric == "fit") {
    if (is.null(indices)) indices <- c("cfi","tli","rmsea","srmr")
    return(tibble(col = paste0(indices, ".m"), panel = indices) %>%
             filter(col %in% names(df)))
  }
  
  if (metric == "chisq_reject")
    return(tibble(col = "chisq.reject", panel = "chisq rejection"))
  
  suffix <- unname(.METRIC_SUFFIX[metric])
  if (is.na(suffix)) suffix <- metric          # m 미등록이어도 "m" 으로 동작
  alias  <- .PARAM_ALIAS[[suffix]]
  
  if (is.null(params)) {
    pat    <- paste0("\\.", suffix, "$")
    cols   <- grep(pat, names(df), value = TRUE)
    panels <- sub(pat, "", cols)
    if (!is.null(alias)) {
      rev_alias <- setNames(names(alias), alias)
      panels <- ifelse(panels %in% names(rev_alias), rev_alias[panels], panels)
    }
    return(tibble(col = cols, panel = panels))
  }
  
  resolved <- if (!is.null(alias))
    ifelse(params %in% names(alias), unname(alias[params]), params) else params
  out <- tibble(col = paste0(resolved, ".", suffix), panel = params)
  dplyr::filter(out, col %in% names(df))        # 존재하는 컬럼만, tibble 반환 보장
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
                        params  = NULL,
                        indices = NULL,
                        models  = c("sum", "EC_latent"),
                        x       = "N_size",
                        x_lab   = "Sample size",
                        y_range = NULL,
                        facets  = c("sd.lambda"),
                        filters = list(),
                        title   = NULL,
                        #max_N   = 1e4,
                        jump_ratio = 3,
                        drop_facets = TRUE,
                        row_strip   = c("right", "left")) {
  row_strip <- match.arg(row_strip)
  
  # 1. restrict to chosen models
  df <- df %>% filter(model %in% models) %>%
    mutate(model = factor(model, levels = models))
  
  # 2. ad-hoc filters
  for (nm in names(filters)) df <- df %>% filter(.data[[nm]] %in%
                                                   filters[[nm]])
  #if (x == "N_size") df <- df %>% filter(N_size <= max_N)
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
  
  # 3b. NEW: equal-spacing x positions + jump detection
  use_discrete <- is.numeric(long[[x]])
  if (use_discrete) {
    xlev <- sort(unique(long[[x]]))
    long$.xpos <- match(long[[x]], xlev)       # 1..k, evenly spaced
    gaps    <- diff(xlev)
    jump_at <- if (length(gaps) > 1)
      which(gaps > jump_ratio * median(gaps)) else integer(0)
    x_aes <- ".xpos"
  } else {
    x_aes <- x
    jump_at <- integer(0)
  }
  
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
  p <- ggplot(long, aes(.data[[x_aes]], value,
                        color = model, shape = model,
                        linetype = loading_quality,
                        group = interaction(model, loading_quality))) +
    geom_line(linewidth = 0.6) + geom_point(size = 2) +
    scale_color_manual(values = PALETTE[models],
                       labels = MODEL_LABELS[models], name = "Model") +
    scale_shape_manual(values = c(16, 17)[seq_along(models)],
                       labels = MODEL_LABELS[models], name = "Model") +
    scale_linetype_manual(values = c(poor = "dashed", vgood = "solid"),
                          labels = c(poor = "Low", vgood = "High"),
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
  
  # 5b. NEW: mark the scale jump between equally-spaced positions
  if (length(jump_at)) {
    brk <- tibble(x = jump_at + 0.5)
    p <- p +
      geom_vline(data = brk, aes(xintercept = x),
                 linetype = "dotted", color = "grey50",
                 inherit.aes = FALSE) +
      geom_text(data = brk, aes(x = x, y = Inf, label = "//"),
                inherit.aes = FALSE, vjust = 1.3, size = 3.2,
                color = "grey50")
  }
  
  # shared vs free y-axis
  free_scale <- if (is.null(y_range)) "free_y" else "fixed"
  
  rhs <- if (length(facets)) paste(facets, collapse = " + ") else "."
  p <- p + facet_grid(as.formula(
    if (nrow(cols) > 1) paste("panel ~", rhs) else paste(". ~", rhs)
  ), scales = free_scale, labeller = .facet_labeller, drop = drop_facets,
     switch = if (row_strip == "left") "y" else NULL)
  # keep the moved row strip looking like the others (grey, outside the axis)
  if (row_strip == "left")
    p <- p + theme(strip.placement = "outside")
  
  # fix y range to the data max (include reference lines so they stay visible)
  if (!is.null(y_range)) {
    yl <- if (identical(y_range, "shared")) {
      vals <- long$value
      if (!is.null(ref) && nrow(ref)) vals <- c(vals, ref$y)
      range(vals, na.rm = TRUE)
    } else y_range
    p <- p + coord_cartesian(ylim = yl, clip = "off") 
  }
  
  # 6. NEW: discrete-position axis with real value labels
  if (use_discrete)
    p <- p + scale_x_continuous(breaks = seq_along(xlev), labels = xlev,
                                expand = expansion(mult = 0.05))
  else if (is.numeric(long[[x]]))
    p <- p + scale_x_continuous(breaks = sort(unique(long[[x]])))
  
  p
}


#--- Blank empty facet cells + aligned stacking -------------------
# Stacking a T=5 plot (8 spec*J columns) over a T=3 plot (only 4 of those
# columns have data) via patchwork distorts cell size: the T=3 panels get
# stretched to the full width. We instead keep T=3 at the full 8 columns
# (drop_facets = FALSE) and physically *erase* the empty cells -- panel,
# top strip header and bottom x-axis are removed while the column width is
# preserved -- then rbind the two gtables so kept cells stay identical in
# size and line up column-for-column. patchwork can't carry an 8-column
# gtable (wrap_ggplot_grob rejects it), so we assemble gtables directly.

# Build the gtable for `p` with its empty cells blanked. `legend` controls the
# legend position ("none" for the top block, "bottom" for the block that should
# carry the shared legend). Emptiness is read from the data geoms only
# (line/point = layers 1-2), so reference / jump lines that span every panel
# don't keep empty cells alive.
blank_gt <- function(p, legend = "none") {
  p   <- p + theme(legend.position = legend)
  b   <- ggplot_build(p)
  lay <- b$layout$layout
  has <- rep(FALSE, nrow(lay))
  for (li in seq_len(min(2, length(b$data)))) {
    d <- b$data[[li]]
    if (!is.null(d) && nrow(d) && "PANEL" %in% names(d))
      has[as.integer(as.character(d$PANEL))] <- TRUE
  }
  g <- ggplot_gtable(b)
  blank <- function(nm) {
    j <- which(g$layout$name == nm)
    if (length(j)) g$grobs[[j]] <<- ggplot2::zeroGrob()   # <<- edits the outer g
  }
  # panel grobs are named "panel-{COL}-{ROW}" in this ggplot build
  for (i in which(!has)) blank(sprintf("panel-%d-%d", lay$COL[i], lay$ROW[i]))
  # for columns empty in every row, also drop the top strip and bottom x-axis
  for (cc in sort(unique(lay$COL))) if (all(!has[lay$COL == cc])) {
    blank(sprintf("strip-t-%d", cc))
    blank(sprintf("axis-b-%d",  cc))
  }
  g
}

# Stack p_top over p_bot with identical, column-aligned cells and a single
# shared legend at the bottom. p_bot should be built with drop_facets = FALSE.
# The legend is ggplot's own bottom legend on p_bot (so it lands at the very
# bottom of the figure) -- this avoids hand-placing a fixed-width legend grob,
# which clipped to stray "M"/"igh" fragments on a narrow device.
# Returns a gtable -- print it, or pass to ggsave()/save_fig() directly.
stack_aligned <- function(p_top, p_bot) {
  g_top <- ggplotGrob(p_top + theme(legend.position = "none"))
  g_bot <- blank_gt(p_bot, legend = "bottom")
  gtable:::rbind_gtable(g_top, g_bot, size = "max")  # "max" equalises widths
}

# Convenience: clear the device first, then draw. grid.draw() does NOT start a
# new page, so calling it repeatedly overlays plots -- always newpage first.
draw_stack <- function(p_top, p_bot) {
  grid::grid.newpage()
  grid::grid.draw(stack_aligned(p_top, p_bot))
}

results_1m <- results %>%
  filter(N_size == 100000)

results_fin <- results %>%
  filter(N_size != 100000)

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

adm_t3 <- plot_metric(results_fin, "admissible",
                      facets  = c("spec", "J"),
                      filters = list(sd.lambda = 0, TT=3),
            title   = "(Time points = 3)")

adm_t5 <- plot_metric(results_fin, "admissible",
            facets  = c("spec", "J"),
            filters = list(sd.lambda = 0, TT=5),
            title   = "(Time points = 5)")

adm_t5 <- adm_t5 + facet_wrap(vars(spec, J), ncol = 4,
                              labeller = .facet_labeller)

(adm_t5 / adm_t3) +
  plot_layout(guides = "collect", heights = c(2, 1)) &
  theme(legend.position = "bottom")




rbias_t3 <- plot_metric(results_fin, "rbias",
            params  = c("phi11", "phi22", "phi12"),
            facets  = c("spec", "J"),
            filters = list( sd.lambda=0, TT=3),
            #y_range = "shared",
            drop_facets = FALSE,          # keep CU modeled / CU overfit as blank columns
            title   = "(Time points = 3)")

rbias_t5 <- plot_metric(results_fin, "rbias",
            params  = c("phi11", "phi22", "phi12"),
            facets  = c("spec", "J"),
            filters = list( sd.lambda=0, TT=5),
            #y_range = "shared",
            title   = "(Time points = 5)")


draw_stack(rbias_t5, rbias_t3)




rmse_t3 <- plot_metric(results_fin, "rmse",
            params  = c("phi11", "phi22", "phi12"),
            facets  = c("spec", "J"),
            filters = list( sd.lambda=0, TT=3),
            y_range = "shared",
            drop_facets = FALSE,          # keep CU modeled / CU overfit as blank columns
            title   = "Relative RMSE by specification (T = 3)")

rmse_t5 <- plot_metric(results_fin, "rmse",
            params  = c("phi11", "phi22", "phi12"),
            facets  = c("spec", "J"),
            filters = list( sd.lambda=0, TT=5),
            y_range = "shared",
            title   = "Relative RMSE by specification (T = 5)")


draw_stack(rmse_t5, rmse_t3)



cov_t3 <- plot_metric(results_fin, "cov",
                       params  = c("phi11", "phi22", "phi12"),
                       facets  = c("spec", "J"),
                       filters = list( sd.lambda=0, TT=3),
                       y_range = "shared",
                       drop_facets = FALSE,   # keep CU modeled / CU overfit as blank columns
                       title   = "95% Coverage by specification (T = 3)")

cov_t5 <- plot_metric(results_fin, "cov",
                       params  = c("phi11", "phi22", "phi12"),
                       facets  = c("spec", "J"),
                       filters = list( sd.lambda=0, TT=5),
                       y_range = "shared",
                       title   = "95% Coverage by specification (T = 5)")


draw_stack(cov_t5, cov_t3)







power_t3 <- plot_metric(results_fin, "power",
                      params  = c("phi11", "phi22", "phi12"),
                      facets  = c("spec", "J"),
                      filters = list( sd.lambda=0, TT=3),
                      #y_range = "shared",
                      drop_facets = FALSE,    # keep CU modeled / CU overfit as blank columns
                      title   = "Power by specification (T = 3)")

power_t5 <- plot_metric(results_fin, "power",
                      params  = c("phi11", "phi22", "phi12"),
                      facets  = c("spec", "J"),
                      filters = list( sd.lambda=0, TT=5),
                      #y_range = "shared",
                      title   = "Power by specification (T = 5)")


draw_stack(power_t5, power_t3)


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
            params  = c("phi11", "phi22", "phi12"),
            facets  = c("spec", "J"),
            filters = list(TT = 5, sd.lambda = 0),
            y_range = "shared",
            title   = "Power by CU specification (T = 5)")
plot_metric(results, "power",
            params  = c("phi11", "phi22", "phi12"),
            facets  = c("spec", "J"),
            filters = list(TT = 3, sd.lambda = 0),
            y_range = "shared",
            title   = "Power by CU specification (T = 3)")

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
