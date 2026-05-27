#--- Setup ---------------------

suppressPackageStartupMessages({
  library(tidyverse)
  library(conflicted)
})
conflicts_prefer(dplyr::select)
conflicts_prefer(dplyr::filter)


#--- Load results --------------

load_results <- function(dir = "study1_results2") {
  files <- list.files(dir, pattern = "\\.RData$", full.names = TRUE)
  rows  <- lapply(files, function(f) {
    e <- new.env(); load(f, envir = e)
    keep <- setdiff(names(e$row), names(e$result$summary_table))
    bind_cols(e$row[, keep, drop = FALSE], e$result$summary_table)
  })
  out <- bind_rows(rows) %>%
    mutate(admissible.m = 1 - err.perc)
  
  # CU on/off label
  if ("rho" %in% names(out)) {
    out <- out %>%
      mutate(rho_lab = factor(
        recode(as.character(rho), "0" = "Basic", "0.15" = "Cor.Unique"),
        levels = c("Basic", "Cor.Unique")))
  }
  # Drop Study 2 residual cells if column exists
  if ("mni_location" %in% names(out))
    out <- out %>% filter(mni_location != "residual")
  out
}

results <- load_results()


#--- Config --------------------
# Okabe-Ito, colorblind-friendly
PALETTE      <- c(sum = "#0072B2", EC_latent = "#E69F00")
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
.facet_labeller <- labeller(
  J               = function(x) paste0("J = ", x),
  TT              = function(x) paste0("T = ", x),
  N_size          = function(x) paste0("N = ", x),
  rho             = function(x) paste0("rho = ", x),
  sd.lambda       = function(x) c("0" = "tau-equiv", "0.05" = "congeneric")[as.character(x)],
  loading_quality = function(x) c(poor = "Low R^2", vgood = "High R^2")[x],
  rho_lab         = label_value,    # already "Basic" / "Cor.Unique"
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
                        facets  = c("loading_quality", "sd.lambda"),
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
  
  grp       <- intersect(c(x, "model", facets), names(df))
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
                        color = model, linetype = model, shape = model,
                        group = model)) +
    geom_line(linewidth = 0.6) + geom_point(size = 2) +
    scale_color_manual(values = PALETTE[models],
                       labels = MODEL_LABELS[models], name = "Model") +
    scale_linetype_manual(values = c("solid", "dashed")[seq_along(models)],
                          labels = MODEL_LABELS[models], name = "Model") +
    scale_shape_manual(values = c(16, 17)[seq_along(models)],
                       labels = MODEL_LABELS[models], name = "Model") +
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


#--- Usage ---------------------

# 1. CU effect, side by side - one column per CU state
#    rows = growth params, columns = Basic vs Cor.Unique
plot_metric(results, "rbias",
            params = c("k1","k2","phi11","phi22","phi12"),
            facets = c("rho_lab", "loading_quality", "sd.lambda"),
            filters = list(J = 3, TT = 5),
            y_range = "shared")

plot_metric(results, "rbias",
            params = c("k1","k2","phi11","phi22","phi12"),
            facets = c("rho_lab", "loading_quality", "sd.lambda"),
            filters = list(J = 3, TT = 3))


# 2. CU x R^2 -- does CU damage worsen under low R^2?
plot_metric(results, "rmse",
            params = c("k1","k2","phi11","phi22","phi12"),
            facets = c("rho_lab", "loading_quality", "sd.lambda"),
            filters = list(J = 3, TT = 3),
            y_range = "shared")
plot_metric(results, "rmse",
            params = c("k1","k2","phi11","phi22","phi12"),
            facets = c("rho_lab", "loading_quality", "sd.lambda"),
            filters = list(J = 3, TT = 5),
            y_range = "shared")

# 3. CU x congeneric SD -- two main bias sources at once
plot_metric(results, "rbias",
            params = c("k1","k2","phi11","phi22","phi12"),
            facets  = c("rho_lab", "sd.lambda"),
            filters = list(loading_quality = "vgood", J = 5, TT = 5))


# 4. CU x design size (J, TT) -- does more indicators / timepoints
#    rescue performance when CU is present?
plot_metric(results, "cov",
            params = c("k1","k2","phi11","phi22","phi12"),
            facets  = c("rho_lab", "J"),
            filters = list(loading_quality = "poor", sd.lambda = 0, TT = 5))

plot_metric(results, "cov",
            params = c("k1","k2","phi11","phi22","phi12"),
            facets  = c("rho_lab", "TT"),
            filters = list(loading_quality = "vgood", sd.lambda = 0, J = 5))

# 5. CU effect on admissibility
plot_metric(results, "admissible",
            facets = c("rho_lab", "loading_quality","sd.lambda","J","TT"),
            filters = list(sd.lambda = 0, rho=.15))

plot_metric(results, "admissible",
            facets = c("rho_lab", "loading_quality","sd.lambda","J","TT"),
            filters = list(sd.lambda = 0, rho=0))

# 6. CU effect on fit indices -- where the EC vs sum gap should be most visible
plot_metric(results, "fit",
            indices = c("cfi","tli","rmsea","srmr"),
            facets  = c("rho_lab","loading_quality"),
            filters = list(sd.lambda = 0,
                           J = 3, TT = 5))

plot_metric(results, "fit",
            indices = c("cfi","tli","rmsea","srmr"),
            facets  = c("rho_lab","loading_quality"),
            filters = list(sd.lambda = 0,
                           J = 3, TT = 3))

plot_metric(results, "fit",
            indices = c("cfi","tli","rmsea","srmr"),
            facets  = c("rho_lab","loading_quality"),
            filters = list(sd.lambda = 0,
                           J = 5, TT = 3))

plot_metric(results, "fit",
            indices = c("cfi","tli","rmsea","srmr"),
            facets  = c("rho_lab","loading_quality"),
            filters = list(sd.lambda = 0,
                           J = 5, TT = 5))



# 7. CU effect on chi-square rejection
plot_metric(results, "chisq_reject",
            facets  = c("rho_lab", "loading_quality","J","TT"),
            filters = list(sd.lambda = 0))

# 8. Power under CU vs Basic
plot_metric(results, "power",
            params = c("k1","k2","phi11","phi22","phi12"),
            facets  = c("rho_lab", "loading_quality"),
            filters = list(sd.lambda = 0, J = 3, TT = 3),
            y_range = "shared")
