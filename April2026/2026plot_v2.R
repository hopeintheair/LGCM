#-------------------------------
#--- Load the Results File -----
#-------------------------------
# 2026 데이터 + 2025plot2 와 일관된 디자인:
#   - Times New Roman / theme_bw / 작은 base_size
#   - x축: "Sample Size", 양옆 expand
#   - 모델 라벨: latent → "Latent 2-LGCM", sum → "Sum 1-LGCM", sum_scaled → "Scaled-Sum 1-LGCM"
#   - rho 값: 0 → "Basic", 0.15 → "Cor.Unique"
#   - metric 표시 라벨: RMSE / Relative Bias / Coverage / Power / Mean

suppressPackageStartupMessages({
  library(patchwork)
  library(scales)
  library(tidyverse)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(rlang)
  library(stringr)
  library(conflicted)
})

conflicts_prefer(dplyr::select)
conflicts_prefer(dplyr::filter)

result_dir <- "."
my_files <- list.files(path = "results_newer", pattern = "\\.RData$", full.names = TRUE)

all_data_list <- list()
for (file in my_files) {
  temp_env <- new.env()
  load(file, envir = temp_env)
  current_row <- temp_env$row
  current_res <- temp_env$result$summary_table
  cols_to_add <- setdiff(names(current_row), names(current_res))
  combined_df <- bind_cols(current_row[, cols_to_add, drop = FALSE], current_res)
  all_data_list[[file]] <- combined_df
}

final_df <- bind_rows(all_data_list)
final_df <- final_df %>%
  mutate(admissible.m = 1 - err.perc) %>%
  filter(mni_location != "residual") %>%
  mutate(
    rho = recode(as.character(rho),
                 "0"      = "Basic",
                 "0.15"   = "Cor.Unique",
                 "no rho" = "Basic",
                 .default = as.character(rho)),
    rho = factor(rho, levels = c("Basic", "Cor.Unique"))
  )

#-------------------------------
# SHARED HELPERS
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

.default_models <- c("latent", "sum", "sum_scaled")

.guess_sim_factors <- function(df) {
  id_cols <- setdiff(names(df), grep("\\.", names(df), value = TRUE))
  exclude <- c("nrep", "seed", "job_id", "N", "nrep_ok")
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
  "latent"     = "Latent 2-LGCM",
  "sum"        = "Sum 1-LGCM",
  "sum_scaled" = "Scaled-Sum 1-LGCM",
  "composite"  = "Sum 1-LGCM"
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

## 1. plot_sim  — growth model parameters
plot_sim <- function(df,
                     metric,
                     params      = NULL,
                     models      = .default_models,
                     x_var       = "N_size",
                     x_lab       = "Sample Size",
                     color_var   = "model",
                     facet_rows  = "param",
                     facet_cols  = c("mni_location", "mni_size"),
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

  default_pool <- c("phi11", "phi22", "phi12", "cor12", "k1", "k2")
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

  if (!is.null(hline)) {
    if (is.list(hline)) {
      # param 별 named list → panel 별 hline
      hline_df <- do.call(rbind, lapply(names(hline), function(p_name) {
        vals <- hline[[p_name]]
        if (is.null(vals) || !p_name %in% levels(long$param)) return(NULL)
        data.frame(param = p_name, yint = vals, stringsAsFactors = FALSE)
      }))
      if (!is.null(hline_df) && nrow(hline_df) > 0) {
        hline_df$param <- factor(hline_df$param, levels = levels(long$param))
        p <- p + geom_hline(data = hline_df,
                            aes(yintercept = yint),
                            linetype = 2, color = "grey60", linewidth = 0.3,
                            inherit.aes = FALSE)
      }
    } else {
      p <- p + geom_hline(yintercept = hline, linetype = 2,
                          color = "grey60", linewidth = 0.3)
    }
  }
  if (is.numeric(long[[x_var]]))
    p <- p + scale_x_continuous(breaks = sort(unique(long[[x_var]])),
                                expand = expansion(mult = c(0.12, 0.12)))

  lhs <- if (length(facet_rows)) paste(facet_rows, collapse = " + ") else "."
  rhs <- if (length(facet_cols)) paste(facet_cols, collapse = " + ") else "."
  if (!(lhs == "." && rhs == "."))
    p <- p + facet_grid(as.formula(paste(lhs, "~", rhs)), scales = "free_y")
  p
}


## 2. plot_sim_multi  — 여러 metric 을 row facet 으로 쌓기
plot_sim_multi <- function(df,
                           metrics,
                           params      = NULL,
                           models      = .default_models,
                           x_var       = "N_size",
                           x_lab       = "Sample Size",
                           color_var   = "model",
                           facet_rows  = NULL,
                           facet_cols  = c("mni_location", "mni_size"),
                           palette     = "okabe_ito",
                           aggregate   = TRUE,
                           sim_factors = NULL,
                           hline       = NULL,
                           title       = NULL,
                           ...) {

  stopifnot(length(metrics) >= 1)
  res <- .apply_filters(df, models, list(...))
  df  <- res$df; filters <- res$filters
  if (is.null(sim_factors)) sim_factors <- .guess_sim_factors(df)
  default_pool <- c("phi11", "phi22", "phi12", "cor12", "k1", "k2")

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

  default_hline_for <- function(mt) {
    mk <- .resolve_metric(mt)
    switch(mk,
           rbias = c(-0.10, 0, 0.10),
           bias  = 0,
           power = 0.80,
           cov   = 0.95,
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

## 3. plot_fit  — model fit indices
.FIT_INDICES <- c("cfi", "tli", "rmsea", "srmr", "chisq")

.fit_hlines <- list(
  cfi   = c(0.90, 0.95),
  tli   = c(0.90, 0.95),
  rmsea = c(0.06, 0.08),
  srmr  = c(0.08),
  chisq = NULL,
  chisq.reject = c(0.05),
  df    = NULL,
  pvalue.chi = c(0.05)
)

plot_fit <- function(df,
                     stat        = "m",
                     indices     = NULL,
                     models      = .default_models,
                     x_var       = "N_size",
                     x_lab       = "Sample Size",
                     color_var   = "model",
                     facet_cols  = c("mni_location", "mni_size"),
                     palette     = "okabe_ito",
                     aggregate   = TRUE,
                     sim_factors = NULL,
                     hline       = NULL,
                     title       = NULL,
                     ...) {

  stopifnot(stat %in% c("m", "sd", "reject"))

  res <- .apply_filters(df, models, list(...))
  df  <- res$df; filters <- res$filters
  if (is.null(sim_factors)) sim_factors <- .guess_sim_factors(df)

  if (stat == "reject") {
    if (is.null(indices)) indices <- "chisq"
    target_cols <- paste0(indices, ".reject")
    target_cols <- intersect(target_cols, names(df))
    if (length(target_cols) == 0) stop("reject 컬럼을 찾지 못했습니다.")
    idx_labels <- sub("\\.reject$", "", target_cols)
    y_label <- "rejection rate"
  } else {
    if (is.null(indices)) indices <- c("cfi", "tli", "rmsea", "srmr")
    target_cols <- paste0(indices, ".", stat)
    target_cols <- intersect(target_cols, names(df))
    if (length(target_cols) == 0)
      stop(sprintf("stat='%s' 에 해당하는 fit 컬럼을 찾지 못했습니다.", stat))
    idx_labels <- sub(paste0("\\.", stat, "$"), "", target_cols)
    y_label <- switch(stat, "m" = "mean", "sd" = "SD")
  }

  id_cols <- unique(c(sim_factors, intersect(names(df),
                                             c(x_var, color_var, facet_cols, names(filters)))))

  long <- df %>%
    select(all_of(id_cols), all_of(target_cols)) %>%
    pivot_longer(cols      = all_of(target_cols),
                 names_to  = "index",
                 values_to = "value") %>%
    mutate(index = sub(paste0("\\.", stat, "$|\\.reject$"), "", index),
           index = factor(index, levels = idx_labels))

  controlled <- unique(c(x_var, color_var, facet_cols,
                         "index", names(filters)))
  long <- .detect_and_aggregate(long, sim_factors, controlled, aggregate)

  if (is.null(hline)) {
    hline_df <- do.call(rbind, lapply(idx_labels, function(idx) {
      vals <- .fit_hlines[[idx]]
      if (stat == "reject") vals <- .fit_hlines[["chisq.reject"]]
      if (stat == "sd") vals <- NULL
      if (is.null(vals)) return(NULL)
      data.frame(index = idx, yint = vals, stringsAsFactors = FALSE)
    }))
    if (!is.null(hline_df) && nrow(hline_df) > 0) {
      hline_df$index <- factor(hline_df$index, levels = idx_labels)
    }
  } else {
    hline_df <- data.frame(
      index = rep(idx_labels, each = length(hline)),
      yint  = rep(hline, times = length(idx_labels)),
      stringsAsFactors = FALSE
    )
    hline_df$index <- factor(hline_df$index, levels = idx_labels)
  }

  k <- length(unique(long[[color_var]]))

  p <- ggplot(long,
              aes(x = .data[[x_var]], y = value,
                  color = .data[[color_var]], linetype = .data[[color_var]],
                  shape = .data[[color_var]], group = .data[[color_var]])) +
    geom_line(linewidth = 0.65) +
    geom_point(size = 2) +
    labs(x = x_lab %||% x_var, y = y_label,
         title = title %||% sprintf("Model fit: %s", y_label)) +
    .cb_scales(color_var, k, palette,
               labels = if (color_var == "model") .MODEL_LABELS else NULL) +
    .theme_2025()

  if (!is.null(hline_df) && nrow(hline_df) > 0) {
    p <- p + geom_hline(data = hline_df,
                        aes(yintercept = yint),
                        linetype = 2, color = "grey60", linewidth = 0.3,
                        inherit.aes = FALSE)
  }

  if (is.numeric(long[[x_var]]))
    p <- p + scale_x_continuous(breaks = sort(unique(long[[x_var]])),
                                expand = expansion(mult = c(0.12, 0.12)))

  rhs <- if (length(facet_cols)) paste(facet_cols, collapse = " + ") else "."
  p <- p + facet_grid(as.formula(paste("index ~", rhs)), scales = "free_y")

  p
}


#-------------------------------
# Example calls
#-------------------------------

plot_sim(final_df, metric = "rmse", model = c("latent", "sum_scaled"),
         facet_cols = c("rho", "mni_location", "mni_size"),
         palette = "tritan", params = c("phi11", "phi22", "phi12", "cor12", "k2"))

plot_sim(final_df, metric = "rbias", model = c("latent", "sum_scaled"),
         facet_cols = c("rho", "mni_location", "mni_size"),
         palette = "tritan", params = c("phi11", "phi22", "phi12", "cor12", "k2"),
         filter_rho = "Basic")

plot_sim(final_df, metric = "rbias", model = c("latent", "sum_scaled"),
         facet_cols = c("rho", "mni_location", "mni_size"),
         palette = "tritan", params = c("phi11", "phi22", "phi12", "cor12", "k2"))

plot_sim(final_df, metric = "cov", model = c("latent", "sum_scaled"),
         hline = c(.9),
         facet_cols = c("rho", "mni_location", "mni_size"),
         palette = "tritan", params = c("phi11", "phi22", "phi12", "cor12", "k1", "k2"),
         filter_rho = "Basic")

plot_sim(final_df, metric = "power", model = c("latent", "sum_scaled"), hline = c(.8),
         facet_cols = c("mni_location", "mni_size"),
         facet_rows = c("param", "loading_quality", "rho"),
         palette = "tritan", params = c("phi11", "phi22", "phi12"),
         title = "Statistical Power")

plot_sim(final_df, metric = "power", model = c("latent", "sum_scaled"), hline = c(.8),
         facet_cols = c("mni_location", "mni_size"),
         facet_rows = c("param", "loading_quality", "rho"),
         palette = "tritan", params = c("k1", "k2"),
         title = "Statistical Power")

# k2 — bias
plot_sim(final_df, metric = "bias", params = "k2", model = c("latent", "sum_scaled"),
         facet_cols = c("rho", "mni_location", "mni_size"), hline = c(-0.02, 0.02))
# 분산 역산
final_df <- final_df %>%
  mutate(k2.var = k2.rmse^2 - k2.bias^2)
plot_sim(final_df, metric = "var", params = "k2", model = c("latent", "sum_scaled"),
         facet_cols = c("rho", "mni_location", "mni_size"))

plot_sim(final_df, metric = "m", model = c("latent", "sum"),
         params = c("phi11", "phi22", "phi12"), palette = "tritan",
         facet_rows = c("param", "loading_quality"),
         facet_cols = c("rho", "mni_location", "mni_size"),
         title = "Variance/Covariance Estimates: Growth Intercept & Slope")

plot_sim(final_df, metric = "m", model = c("latent", "sum"),
         params = c("k1", "k2"), palette = "tritan",
         facet_rows = c("param", "loading_quality"),
         facet_cols = c("rho", "mni_location", "mni_size"),
         title = "Mean Estimates: Growth Intercept & Slope")


plot_sim(final_df, metric = "m", model = c("latent", "sum"),
         facet_cols = c("rho", "mni_location", "mni_size"),
         palette = "tritan", params = c("cor12"))

plot_sim(final_df, metric = "rbias", model = c("latent", "sum"),
         facet_cols = c("rho", "mni_location", "mni_size"),
         palette = "tritan", params = c("cor12"))

plot_sim(final_df, metric = "m", model = c("latent", "sum"), palette = "tritan",
         params = c("rho_IS"),
         facet_cols = c("rho", "mni_location", "mni_size"))

plot_sim(final_df, metric = "m", model = c("latent", "sum"), palette = "tritan",
         params = c("rho_IS", "d_slope"),
         facet_rows = c("param", "loading_quality"),
         facet_cols = c("rho", "mni_location", "mni_size"),
         title = "Scale-Free Measures")

plot_sim(final_df, metric = "m", model = c("latent", "sum"), palette = "tritan",
         params = c("pvalue.chi",  "rmsea"),
         facet_rows = c("param", "loading_quality"),
         facet_cols = c("rho", "mni_location", "mni_size"),
         hline = list(cfi   = 0.95,                                           
                      tli   = 0.95,                                           
                      rmsea = 0.06,                                           
                      srmr  = 0.08),
         title = "Model Fit Indices")

plot_sim(final_df, metric = "m", model = c("latent", "sum"), palette = "tritan",
         params = c( "cfi", "tli", "srmr"),
         facet_rows = c("param", "loading_quality"),
         facet_cols = c("rho", "mni_location", "mni_size"),
         hline = list(cfi   = 0.95,                                           
                      tli   = 0.95,                                           
                      rmsea = 0.06,                                           
                      srmr  = 0.08),
         title = "Model Fit Indices")


plot_sim(final_df, metric = "m", model = c("latent",  "sum"), palette = "tritan",
         params = c("admissible"),
         facet_cols = c("rho","loading_quality"),         
         facet_rows = c("mni_location", "mni_size"),
         title = "Admissible Results Rate")

plot_fit(final_df, stat = "reject", model = c("latent", "sum_scaled", "sum"), palette = "tritan",
         facet_cols = c("rho", "mni_location", "mni_size"))
