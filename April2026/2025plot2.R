#-------------------------------
#--- Load the Results File -----
#-------------------------------
# 2025 data + 2026 framework (palette / helpers / plot_sim / plot_fit)
# 디자인(폰트/테마/플랏 크기)은 2025plot.R 스타일(Times New Roman, theme_bw)을 유지.

suppressPackageStartupMessages({
  library(patchwork)
  library(scales)
  library(tidyverse)
  library(dplyr)
  library(tidyr)
  library(purrr)
  library(stringr)
  library(forcats)
  library(ggplot2)
  library(rlang)
  library(conflicted)
})

conflicts_prefer(dplyr::select)
conflicts_prefer(dplyr::filter)

#-------------------------------
# 1) 2025 데이터 로드 및 long 변환
#-------------------------------
load("~/Documents/GitHub/SOG_SIM/Qual_Code/1217/5i5t_ST_basic3.RData")
load("~/Documents/GitHub/SOG_SIM/Qual_Code/1217/rhos/5i5t_ST_rho3.RData")

want_cols <- c(
  "phi11.power","phi12.power","phi22.power","kappa1.power","kappa2.power",
  "phi11.R.bias","phi12.R.bias","phi22.R.bias","kappa1.R.bias","kappa2.R.bias",
  "phi11.rel.rmse","phi12.rel.rmse","phi22.rel.rmse","kappa1.rel.rmse","kappa2.rel.rmse",
  "phi11.coverage","phi12.coverage","phi22.coverage","kappa1.coverage","kappa2.coverage",
  "na.perc","admissible"  # convergence rate
)

obj_names <- ls(pattern = "^(latent|composite)_.+_ST_\\d+_(basic|rho|sd)$")

one_to_long <- function(nm) {
  x <- get(nm, envir = .GlobalEnv)
  sm <- x[["summary"]]
  if (is.null(sm)) {
    warning("No [[summary]] in: ", nm); return(NULL)
  }
  sm <- as.data.frame(sm)
  if (nrow(sm) > 1) sm <- sm[1, , drop = FALSE]

  cn <- colnames(sm)
  cn <- str_replace(cn, "^alpha1\\b", "kappa1")
  cn <- str_replace(cn, "^alpha2\\b", "kappa2")
  colnames(sm) <- cn

  keep <- intersect(want_cols, colnames(sm))
  sm <- sm[, keep, drop = FALSE]

  m <- str_match(
    nm,
    "^(latent|composite)_([^_]+)_([^_]+)_([^_]+)_ST_(\\d+)_(basic|rho|sd)$"
  )
  if (anyNA(m)) {
    warning("Name does not match expected pattern: ", nm); return(NULL)
  }

  meta <- tibble(
    object_name = nm,
    model      = m[,2],
    design     = m[,3],
    slope_var  = m[,4],
    quality    = m[,5],
    nsize      = as.integer(m[,6]),
    type       = m[,7]
  )

  special_keys <- c("na.perc", "admissible")
  out <- sm %>%
    tibble::as_tibble() %>%
    pivot_longer(cols = everything(), names_to = "key", values_to = "value") %>%
    mutate(
      key    = if_else(key %in% special_keys, "na.perc", key),
      param  = if_else(key == "na.perc", "admissible", str_extract(key, "^[^.]+")),
      metric = if_else(key == "na.perc", "m",          str_extract(key, "(?<=\\.).+$"))
    ) %>%
    select(-key) %>%
    bind_cols(meta, .)
  out
}

sim_long_orig <- map_dfr(obj_names, one_to_long) %>%
  select(model, type, design, slope_var, quality, nsize, param, metric, value, object_name) %>%
  arrange(model, type, slope_var, quality, nsize, param, metric)

#-------------------------------
# 2) 2026 프레임워크가 기대하는 wide 포맷으로 변환
#    metric 별 alias: rel.rmse→rmse, R.bias→rbias, coverage→cov, power→power, m→m
#-------------------------------
metric_short_map <- c(
  "rel.rmse" = "rmse",
  "R.bias"   = "rbias",
  "coverage" = "cov",
  "power"    = "power",
  "m"        = "m"
)

final_df <- sim_long_orig %>%
  mutate(metric_short = unname(metric_short_map[metric])) %>%
  filter(!is.na(metric_short)) %>%
  mutate(col = paste0(param, ".", metric_short)) %>%
  select(model, type, design, slope_var, quality, nsize, col, value) %>%
  pivot_wider(names_from = col, values_from = value) %>%
  mutate(
    quality   = recode(quality,
                       "poor"  = "Poor",
                       "good"  = "Good",
                       "vgood" = "Very Good"),
    slope_var = recode(slope_var,
                       "small"  = "Small",
                       "medium" = "Medium",
                       "large"  = "Large"),
    type      = recode(type,
                       "basic" = "Basic",
                       "rho"   = "Cor.Unique",
                       "sd"    = "SD"),
    quality   = factor(quality,   levels = c("Poor", "Good", "Very Good")),
    slope_var = factor(slope_var, levels = c("Small", "Medium", "Large")),
    type      = factor(type,      levels = c("Basic", "Cor.Unique", "SD")),
    model     = factor(model,     levels = c("latent", "composite"))
  )

#-------------------------------
# 3) 2026 SHARED HELPERS  (palette / scales / filters)
#    + 2025 디자인(Times New Roman, theme_bw) 적용
#-------------------------------

.OKABE_ITO  <- c("#E69F00", "#56B4E9", "#009E73", "#F0E442",
                 "#0072B2", "#D55E00", "#CC79A7", "#000000")
.CB3        <- c("#0072B2", "#E69F00", "#009E73")
.CB3_TRITAN <- c("#0072B2", "#E69F00", "#CC79A7")

.metric_alias <- c(
  "R.RMSE" = "rmse", "RMSE" = "rmse", "rel.rmse" = "rmse",
  "R.bias" = "rbias", "Rbias" = "rbias", "rel.bias" = "rbias",
  "bias"   = "bias",
  "coverage" = "cov", "cov.rate" = "cov",
  "power"  = "power",
  "mean"   = "m"
)
.resolve_metric <- function(metric) {
  if (metric %in% names(.metric_alias)) .metric_alias[[metric]] else metric
}

# 표시용 라벨 (facet strip / y축 / 타이틀)
.METRIC_DISPLAY <- c(
  "rmse"  = "RMSE",
  "rbias" = "Relative Bias",
  "bias"  = "Bias",
  "cov"   = "Coverage",
  "power" = "Power",
  "m"     = "Mean"
)
.metric_label <- function(metric) {
  k <- .resolve_metric(metric)
  if (k %in% names(.METRIC_DISPLAY)) unname(.METRIC_DISPLAY[k]) else metric
}

.extract_filters <- function(dots) {
  if (length(dots) == 0) return(list())
  nm <- names(dots)
  if (is.null(nm)) return(list())
  keep <- startsWith(nm, "filter_")
  if (!any(keep)) return(list())
  out <- dots[keep]
  names(out) <- sub("^filter_", "", names(out))
  out
}

.available_params <- function(df, metric, restrict = NULL) {
  cols <- grep(paste0("\\.", metric, "$"), names(df), value = TRUE)
  ps   <- sub(paste0("\\.", metric, "$"), "", cols)
  if (!is.null(restrict)) ps <- intersect(ps, restrict)
  ps
}

# 2025 데이터 기준 기본 모델
.default_models <- c("latent", "composite")

.guess_sim_factors <- function(df) {
  id_cols <- setdiff(names(df), grep("\\.", names(df), value = TRUE))
  exclude <- c("nrep", "seed", "job_id", "N", "nrep_ok", "object_name", "design")
  setdiff(id_cols, exclude)
}

.get_palette <- function(palette = "okabe_ito", k) {
  if (is.character(palette) && length(palette) == 1) {
    pal <- switch(palette,
                  "okabe_ito" = if (k <= 3) .CB3 else .OKABE_ITO,
                  "tritan"    = .CB3_TRITAN,
                  "viridis"   = if (requireNamespace("viridisLite", quietly = TRUE))
                    viridisLite::viridis(k, option = "D", end = 0.9)
                  else .OKABE_ITO,
                  "turbo"     = if (requireNamespace("viridisLite", quietly = TRUE))
                    viridisLite::viridis(k, option = "H", end = 0.95)
                  else .OKABE_ITO,
                  .OKABE_ITO
    )
  } else if (is.character(palette) && length(palette) > 1) {
    pal <- palette
  } else {
    pal <- .OKABE_ITO
  }
  if (length(pal) < k) pal <- rep_len(pal, k)
  pal[seq_len(k)]
}

.LT_SEQ    <- c("solid", "twodash", "dashed", "dotted",
                "dotdash", "longdash", "1F", "F1")
.SHAPE_SEQ <- c(16, 17, 15, 18, 4, 8, 3, 7)

.MODEL_LABELS <- c(
  "latent"    = "Latent 2-LGCM",
  "composite" = "Sum 1-LGCM"
)

.legend_title <- function(var) {
  if (is.null(var) || !nzchar(var)) return(var)
  paste0(toupper(substr(var, 1, 1)), substring(var, 2))
}

.cb_scales <- function(color_var, k, palette = "okabe_ito", labels = NULL) {
  pal  <- .get_palette(palette, k)
  name <- .legend_title(color_var)
  args_col <- list(values = pal,                    name = name)
  args_lt  <- list(values = .LT_SEQ[seq_len(k)],    name = name)
  args_sh  <- list(values = .SHAPE_SEQ[seq_len(k)], name = name)
  if (!is.null(labels)) {
    args_col$labels <- labels
    args_lt$labels  <- labels
    args_sh$labels  <- labels
  }
  list(
    do.call(scale_color_manual,    args_col),
    do.call(scale_linetype_manual, args_lt),
    do.call(scale_shape_manual,    args_sh)
  )
}

# 2025 디자인: Times New Roman 기반 theme_bw
.BASE_FAMILY <- "Times New Roman"
.theme_2025 <- function(base_size = 9) {
  theme_bw(base_size = base_size, base_family = .BASE_FAMILY) +
    theme(
      text             = element_text(family = .BASE_FAMILY),
      panel.grid.minor = element_blank(),
      strip.background = element_rect(fill = "grey92"),
      strip.text       = element_text(size = base_size - 1),
      legend.position  = "right",
      panel.spacing.x  = unit(0.8, "lines"),
      axis.text        = element_text(size = base_size - 1),
      axis.text.x      = element_text(angle = 45, hjust = 1, vjust = 1,
                                      size = base_size - 1)
    )
}

## 공통: 모델 필터 + 추가 필터 적용
.apply_filters <- function(df, models, dots) {
  if (!is.null(models) && "model" %in% names(df)) {
    miss <- setdiff(models, unique(as.character(df$model)))
    if (length(miss))
      warning(sprintf("다음 model 수준은 df 에 없어 제외됩니다: %s",
                      paste(miss, collapse = ", ")))
    models <- intersect(models, unique(as.character(df$model)))
    df <- df %>% filter(as.character(model) %in% models) %>%
      mutate(model = factor(model, levels = models))
  }
  filters <- .extract_filters(dots)
  for (fnm in names(filters)) {
    if (!fnm %in% names(df)) {
      warning(sprintf("필터 변수 '%s' 가 df 에 없어 무시합니다.", fnm))
      next
    }
    val <- filters[[fnm]]
    df  <- df %>% filter(.data[[fnm]] %in% val)
    if (is.character(df[[fnm]]) || is.factor(df[[fnm]]))
      df[[fnm]] <- factor(df[[fnm]], levels = val)
  }
  if (nrow(df) == 0) stop("필터 결과가 비어 있습니다.")
  list(df = df, filters = filters)
}

## 공통: 잔여 요인 감지 + 집약
.detect_and_aggregate <- function(long, sim_factors, controlled, aggregate) {
  leftover <- setdiff(sim_factors, controlled)
  if (length(leftover) > 0) {
    keep <- vapply(leftover, function(v) {
      v %in% names(long) && length(unique(long[[v]])) > 1
    }, logical(1))
    leftover <- leftover[keep]
  }
  if (length(leftover) > 0) {
    message(sprintf(
      "잔여 요인 %s → %s",
      paste(leftover, collapse = ", "),
      if (aggregate) "평균 집약" else "별도 선"
    ))
  }
  group_vars <- intersect(controlled, names(long))
  if (aggregate && length(leftover) > 0) {
    long <- long %>%
      group_by(across(all_of(group_vars))) %>%
      summarise(value = mean(value, na.rm = TRUE), .groups = "drop")
  }
  long
}

#-------------------------------
# 4) plot_sim / plot_sim_multi / plot_fit  (2026 framework)
#    — theme 만 2025 디자인으로 교체
#-------------------------------

plot_sim <- function(df,
                     metric,
                     params      = NULL,
                     models      = .default_models,
                     x_var       = "nsize",
                     x_lab       = "Sample Size",
                     color_var   = "model",
                     facet_rows  = "param",
                     facet_cols  = c("type", "quality"),
                     palette     = "okabe_ito",
                     aggregate   = TRUE,
                     sim_factors = NULL,
                     hline       = NULL,
                     title       = NULL,
                     ...) {

  stopifnot(is.data.frame(df), is.character(metric), length(metric) == 1)
  mkey <- .resolve_metric(metric)

  res <- .apply_filters(df, models, list(...))
  df  <- res$df; filters <- res$filters

  default_pool <- c("phi11", "phi22", "phi12", "kappa1", "kappa2")
  avail <- if (is.null(params) && mkey == "m")
    .available_params(df, mkey, restrict = default_pool)
  else
    .available_params(df, mkey)
  if (length(avail) == 0)
    stop(sprintf("metric '%s' (-> '.%s') 에 해당하는 컬럼이 없습니다.", metric, mkey))
  if (is.null(params)) { params <- avail
  } else {
    miss <- setdiff(params, avail)
    if (length(miss)) warning(sprintf("metric '%s' 에 없는 파라미터: %s",
                                      metric, paste(miss, collapse = ", ")))
    params <- intersect(params, avail)
    if (length(params) == 0) stop("유효한 파라미터가 없습니다.")
  }

  value_cols <- paste0(params, ".", mkey)
  if (is.null(sim_factors)) sim_factors <- .guess_sim_factors(df)
  id_cols <- unique(c(sim_factors, intersect(names(df),
                                             c(x_var, color_var,
                                               setdiff(c(facet_rows, facet_cols), "param"),
                                               names(filters)))))

  long <- df %>%
    select(all_of(id_cols), all_of(value_cols)) %>%
    pivot_longer(cols = all_of(value_cols),
                 names_to = "param", values_to = "value") %>%
    mutate(param = sub(paste0("\\.", mkey, "$"), "", param),
           param = factor(param, levels = params))

  controlled <- unique(c(x_var, color_var,
                         setdiff(c(facet_rows, facet_cols), "param"),
                         "param", names(filters)))
  long <- .detect_and_aggregate(long, sim_factors, controlled, aggregate)

  if (is.null(hline)) {
    hline <- switch(mkey,
                    rbias = c(-0.10, 0, 0.10),
                    bias  = 0,
                    power = 0.80,
                    cov   = 0.95,
                    NULL)
  }

  k <- length(unique(long[[color_var]]))
  p <- ggplot(long,
              aes(x = .data[[x_var]], y = value,
                  color = .data[[color_var]], linetype = .data[[color_var]],
                  shape = .data[[color_var]], group = .data[[color_var]])) +
    geom_line(linewidth = 0.65) + geom_point(size = 2) +
    labs(x = x_lab %||% x_var, y = .metric_label(metric),
         title = title %||% sprintf("Simulation results: %s",
                                    .metric_label(metric))) +
    .cb_scales(color_var, k, palette,
               labels = if (color_var == "model") .MODEL_LABELS else NULL) +
    .theme_2025()

  if (!is.null(hline))
    p <- p + geom_hline(yintercept = hline, linetype = 2,
                        color = "grey60", linewidth = 0.3)
  if (is.numeric(long[[x_var]]))
    p <- p + scale_x_continuous(breaks = sort(unique(long[[x_var]])),
                                expand = expansion(mult = c(0.12, 0.12)))

  lhs <- if (length(facet_rows)) paste(facet_rows, collapse = " + ") else "."
  rhs <- if (length(facet_cols)) paste(facet_cols, collapse = " + ") else "."
  if (!(lhs == "." && rhs == "."))
    p <- p + facet_grid(as.formula(paste(lhs, "~", rhs)), scales = "free_y")
  p
}

plot_sim_multi <- function(df,
                           metrics,
                           params      = NULL,
                           models      = .default_models,
                           x_var       = "nsize",
                           x_lab       = "Sample Size",
                           color_var   = "model",
                           facet_rows  = NULL,   # metric 앞에 붙는 추가 row facet
                           facet_cols  = c("type", "quality"),
                           palette     = "okabe_ito",
                           aggregate   = TRUE,
                           sim_factors = NULL,
                           hline       = NULL,   # named list: list(rmse = 0, power = 0.8) 또는 단일 vector
                           title       = NULL,
                           ...) {

  stopifnot(length(metrics) >= 1)
  res <- .apply_filters(df, models, list(...))
  df  <- res$df; filters <- res$filters
  if (is.null(sim_factors)) sim_factors <- .guess_sim_factors(df)
  default_pool <- c("phi11", "phi22", "phi12", "kappa1", "kappa2")

  long_all <- lapply(metrics, function(mt) {
    mkey  <- .resolve_metric(mt)
    avail <- if (is.null(params) && mkey == "m")
      .available_params(df, mkey, restrict = default_pool)
    else .available_params(df, mkey)
    use <- if (is.null(params)) avail else intersect(params, avail)
    if (length(use) == 0) return(NULL)
    vc <- paste0(use, ".", mkey)
    id <- unique(c(sim_factors, intersect(names(df),
                                          c(x_var, color_var,
                                            setdiff(c(facet_rows, facet_cols), "param"),
                                            names(filters)))))
    df %>%
      select(all_of(id), all_of(vc)) %>%
      pivot_longer(all_of(vc), names_to = "param", values_to = "value") %>%
      mutate(param = sub(paste0("\\.", mkey, "$"), "", param), metric = mt)
  }) %>% bind_rows()

  if (nrow(long_all) == 0) stop("그릴 데이터가 없습니다.")
  long_all$metric <- factor(long_all$metric, levels = metrics,
                            labels = vapply(metrics, .metric_label, character(1)))

  controlled <- unique(c(x_var, color_var,
                         setdiff(c(facet_rows, facet_cols), "param"),
                         "param", "metric", names(filters)))
  long_all <- .detect_and_aggregate(long_all, sim_factors, controlled, aggregate)

  # metric 별 hline 데이터 만들기
  default_hline_for <- function(mt) {
    mk <- .resolve_metric(mt)
    switch(mk,
           rbias = c(-0.10, 0, 0.10),
           bias  = 0,
           power = 0.80,
           cov   = 0.9,
           NULL)
  }
  hline_df <- do.call(rbind, lapply(metrics, function(mt) {
    vals <- if (is.null(hline)) default_hline_for(mt)
    else if (is.list(hline))   hline[[mt]]
    else                       hline
    if (is.null(vals)) return(NULL)
    data.frame(metric = mt, yint = vals, stringsAsFactors = FALSE)
  }))
  if (!is.null(hline_df) && nrow(hline_df) > 0)
    hline_df$metric <- factor(hline_df$metric, levels = metrics,
                              labels = vapply(metrics, .metric_label, character(1)))

  k <- length(unique(long_all[[color_var]]))
  p <- ggplot(long_all,
              aes(x = .data[[x_var]], y = value,
                  color = .data[[color_var]], linetype = .data[[color_var]],
                  shape = .data[[color_var]], group = .data[[color_var]])) +
    geom_line(linewidth = 0.65) + geom_point(size = 2) +
    labs(x = x_lab %||% x_var, y = "Value",
         title = title %||% "Simulation results") +
    .cb_scales(color_var, k, palette,
               labels = if (color_var == "model") .MODEL_LABELS else NULL) +
    .theme_2025()

  if (!is.null(hline_df) && nrow(hline_df) > 0)
    p <- p + geom_hline(data = hline_df,
                        aes(yintercept = yint),
                        linetype = 2, color = "grey60", linewidth = 0.3,
                        inherit.aes = FALSE)

  lhs <- paste(c("metric", facet_rows), collapse = " + ")
  rhs <- if (length(facet_cols)) paste(facet_cols, collapse = " + ") else "."
  p <- p + facet_grid(as.formula(paste(lhs, "~", rhs)), scales = "free_y")

  if (is.numeric(df[[x_var]]))
    p <- p + scale_x_continuous(breaks = sort(unique(long_all[[x_var]])),
                                expand = expansion(mult = c(0.12, 0.12)))
  p
}

#-------------------------------
# 5) Example calls — 2025plot.R 와 동일한 패널 구성
#    facet: metric + slope_var (행) ~ type + quality (열)
#    Power+Coverage 한 플랏, RMSE+Bias 한 플랏.
#-------------------------------

# phi11
plot_sim_multi(final_df, metrics = c( "cov","power"),
               params = "phi11",
               facet_rows = "slope_var",
               facet_cols = c("type", "quality"),
               title = "Variance of Growth Intercept - Power & Coverage")


plot_sim_multi(final_df, metrics = c("rmse", "rbias"),
               params = "phi11",
               facet_rows = "slope_var",
               facet_cols = c("type", "quality"),
               title = "Variance of Growth Intercept - RMSE & Bias")

plot_sim_multi(final_df, metrics = c("rmse", "rbias"),
               params = "phi12",
               facet_rows = "slope_var",
               facet_cols = c("type", "quality"),
               title = "Covariance of Growth Intercept & Slope — RMSE & Bias")

# phi12
plot_sim_multi(final_df, metrics = c("cov","power"),
               params = "phi12",
               facet_rows = "slope_var",
               facet_cols = c("type", "quality"),
               title = "Cov(Intercept, Slope) — Power & Coverage")

plot_sim_multi(final_df, metrics = c("rmse", "rbias"),
               params = "phi12",
               facet_rows = "slope_var",
               facet_cols = c("type", "quality"),
               title = "Cov(Intercept, Slope) — RMSE & Bias")

# phi22
plot_sim_multi(final_df, metrics = c("cov","power"),
               params = "phi22",
               facet_rows = "slope_var",
               facet_cols = c("type", "quality"),
               title = "Variance of Growth Slope — Power & Coverage")

plot_sim_multi(final_df, metrics = c("rmse", "rbias"),
               params = "phi22",
               facet_rows = "slope_var",
               facet_cols = c("type", "quality"),
               title = "Variance of Growth Slope — RMSE & Bias")

# kappa1 / kappa2 (param 도 같이 facet 으로)
plot_sim_multi(final_df, metrics = c( "cov","power"),
               params = c("kappa1", "kappa2"),
               facet_rows = c("slope_var"),
               facet_cols = c("param", "type", "quality"),
               title = "Mean of Growth Intercept & Slope — Power & Coverage")

plot_sim_multi(final_df, metrics = c("rmse", "rbias"),
               params = c("kappa1", "kappa2"),
               facet_rows = c("slope_var"),
               facet_cols = c("param", "type", "quality"),
               title = "Mean of Growth Intercept & Slope — RMSE & Bias")

plot_sim_multi(final_df, metrics = c("rmse", "rbias"),
               params = c("phi11"),
               facet_rows = c("param", "slope_var"),
               facet_cols = c("type", "quality"),
               title = "Mean of Growth Intercept & Slope — RMSE & Bias")

# Convergence (admissible 비율) — 단독 플랏
plot_sim(final_df, metric = "m", params = "admissible",
         facet_rows = "slope_var",
         facet_cols = c("type", "quality"),
         hline = 0.80,
         title = "Admissible Output Proportion")
