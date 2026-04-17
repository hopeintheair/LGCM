## ═══════════════════════════════════════════════════════════════════
## Composite LGCM: Linear vs. Quadratic vs. Latent Basis
## ═══════════════════════════════════════════════════════════════════
## 이 코드는 analysis_5i5t_MI_v5.R의 composite section (model == "composite")
## 내부 lavaan fitting 부분을 확장합니다.
##
## 기존 코드의 data.comp 생성 이후 (line ~669) 부터 삽입합니다.
## ═══════════════════════════════════════════════════════════════════


## ─────────────────────────────────────────────────────────────────
## Model 1: Linear (기존 그대로)
## ─────────────────────────────────────────────────────────────────
C.linear <- '
  i =~ 1*time1 + 1*time2 + 1*time3 + 1*time4 + 1*time5
  s =~ 0*time1 + 1*time2 + 2*time3 + 3*time4 + 4*time5
  i ~~ (phi1)*i
  s ~~ (phi2)*s
  i ~~ (phi12)*s
  i ~  (a1)*1
  s ~  (a2)*1
'


## ─────────────────────────────────────────────────────────────────
## Model 2: Quadratic
## ─────────────────────────────────────────────────────────────────
## 
## S_t = alpha_I + alpha_S*t + alpha_Q*t^2 + error
##
## 주의: growth factor가 3개 (i, s, q)이므로
##   - 자유모수가 늘어남: mean 3개, variance 3개, covariance 3개 = 9개
##   - 5시점에서 observed mean 5개, observed variance 5개 → identification OK
##   - 하지만 residual variance는 5 - 9 = -4이므로
##     residual variance를 동일하게 제약해야 합니다 (1 free param)
##     또는 일부 시점만 자유추정
## ─────────────────────────────────────────────────────────────────
C.quad <- '
  i =~ 1*time1 + 1*time2 + 1*time3 + 1*time4 + 1*time5
  s =~ 0*time1 + 1*time2 + 2*time3 + 3*time4 + 4*time5
  q =~ 0*time1 + 1*time2 + 4*time3 + 9*time4 + 16*time5

  i ~~ (phi1)*i
  s ~~ (phi2)*s
  q ~~ (phi3)*q
  i ~~ (phi12)*s
  i ~~ (phi13)*q
  s ~~ (phi23)*q

  i ~ (a1)*1
  s ~ (a2)*1
  q ~ (a3)*1

  ## residual variance: 동일 제약 (identification)
  time1 ~~ (theta)*time1
  time2 ~~ (theta)*time2
  time3 ~~ (theta)*time3
  time4 ~~ (theta)*time4
  time5 ~~ (theta)*time5
'


## ─────────────────────────────────────────────────────────────────
## Model 3: Latent Basis (Freed Loading)
## ─────────────────────────────────────────────────────────────────
##
## S_t = alpha_I + alpha_S * b_t + error
##   b_1 = 0 (fixed), b_5 = 1 (fixed), b_2~b_4 자유추정
##
## Growth factor는 2개지만 basis를 3개 자유추정하므로
## mean: 2 + 3 basis = 5 → observed mean 5개 (just identified for means)
## variance: 2 var + 1 cov + 5 resid var = 8, observed var/cov = 15 → OK
## ─────────────────────────────────────────────────────────────────
C.latent.basis <- '
  i =~ 1*time1 + 1*time2 + 1*time3 + 1*time4 + 1*time5
  s =~ 0*time1 + (b2)*time2 + (b3)*time3 + (b4)*time4 + 1*time5

  i ~~ (phi1)*i
  s ~~ (phi2)*s
  i ~~ (phi12)*s

  i ~ (a1)*1
  s ~ (a2)*1
'


## ═══════════════════════════════════════════════════════════════════
## Fitting & Extraction: 세 모형을 모두 fit하는 wrapper
## ═══════════════════════════════════════════════════════════════════

fit_composite_models <- function(data.comp) {
  
  results <- list()
  
  ## --- Linear ---
  fit.lin <- tryCatch(
    lavaan::growth(C.linear, data = data.comp),
    error   = function(e) NULL,
    warning = function(w) suppressWarnings(
      lavaan::growth(C.linear, data = data.comp))
  )
  
  ## --- Quadratic ---
  fit.quad <- tryCatch(
    lavaan::growth(C.quad, data = data.comp),
    error   = function(e) NULL,
    warning = function(w) suppressWarnings(
      lavaan::growth(C.quad, data = data.comp))
  )
  
  ## --- Latent Basis ---
  fit.lb <- tryCatch(
    lavaan::growth(C.latent.basis, data = data.comp),
    error   = function(e) NULL,
    warning = function(w) suppressWarnings(
      lavaan::growth(C.latent.basis, data = data.comp))
  )
  
  ## Extract fit indices for model comparison
  extract_fit <- function(fit, name) {
    if (is.null(fit)) {
      return(data.frame(
        model = name,
        chi2 = NA, df = NA, pvalue = NA,
        cfi = NA, tli = NA, rmsea = NA,
        aic = NA, bic = NA, srmr = NA,
        converged = FALSE
      ))
    }
    fm <- lavaan::fitMeasures(fit,
            c("chisq","df","pvalue","cfi","tli",
              "rmsea","aic","bic","srmr"))
    data.frame(
      model     = name,
      chi2      = fm["chisq"],
      df        = fm["df"],
      pvalue    = fm["pvalue"],
      cfi       = fm["cfi"],
      tli       = fm["tli"],
      rmsea     = fm["rmsea"],
      aic       = fm["aic"],
      bic       = fm["bic"],
      srmr      = fm["srmr"],
      converged = lavaan::lavInspect(fit, "converged")
    )
  }
  
  results$fit_indices <- rbind(
    extract_fit(fit.lin,  "linear"),
    extract_fit(fit.quad, "quadratic"),
    extract_fit(fit.lb,   "latent_basis")
  )
  
  ## Extract growth parameters
  extract_growth <- function(fit, name) {
    if (is.null(fit) || !lavaan::lavInspect(fit, "converged")) {
      return(data.frame(
        model = name,
        a1 = NA, a1.se = NA,
        a2 = NA, a2.se = NA,
        a3 = NA, a3.se = NA,
        phi1 = NA, phi2 = NA, phi3 = NA,
        phi12 = NA, phi13 = NA, phi23 = NA
      ))
    }
    pe <- as.data.frame(lavaan::parameterEstimates(fit))
    pe <- pe[pe$label != "", ]
    gv  <- function(lbl, col) {
      val <- pe[pe$label == lbl, col, drop = TRUE]
      if (length(val) == 0L) NA_real_ else val[1L]
    }
    
    data.frame(
      model = name,
      a1    = gv("a1", "est"),   a1.se  = gv("a1", "se"),
      a2    = gv("a2", "est"),   a2.se  = gv("a2", "se"),
      a3    = gv("a3", "est"),   a3.se  = gv("a3", "se"),
      phi1  = gv("phi1", "est"), phi2   = gv("phi2", "est"),
      phi3  = gv("phi3", "est"),
      phi12 = gv("phi12","est"), phi13  = gv("phi13","est"),
      phi23 = gv("phi23","est")
    )
  }
  
  results$growth_params <- rbind(
    extract_growth(fit.lin,  "linear"),
    extract_growth(fit.quad, "quadratic"),
    extract_growth(fit.lb,   "latent_basis")
  )
  
  ## Extract latent basis estimates (b2, b3, b4)
  if (!is.null(fit.lb) && lavaan::lavInspect(fit.lb, "converged")) {
    pe.lb <- as.data.frame(lavaan::parameterEstimates(fit.lb))
    pe.lb <- pe.lb[pe.lb$label != "", ]
    gv.lb <- function(lbl) {
      val <- pe.lb[pe.lb$label == lbl, "est", drop = TRUE]
      if (length(val) == 0L) NA_real_ else val[1L]
    }
    results$basis_estimates <- c(
      b1 = 0,                  # fixed
      b2 = gv.lb("b2"),
      b3 = gv.lb("b3"),
      b4 = gv.lb("b4"),
      b5 = 1                   # fixed
    )
  } else {
    results$basis_estimates <- rep(NA_real_, 5)
  }
  
  ## Store fit objects for further inspection
  results$fits <- list(
    linear = fit.lin, quadratic = fit.quad, latent_basis = fit.lb
  )
  
  results
}


## ═══════════════════════════════════════════════════════════════════
## Population truth for latent basis: 
## "올바른 basis"가 무엇인지 계산
## ═══════════════════════════════════════════════════════════════════
##
## MI 위반 시 composite의 true mean trajectory는:
##   E[S_t] = c_t + s_t * (kappa_I + t * kappa_S)
##
## 이것을 latent basis model로 표현하면:
##   E[S_t] = alpha_I + alpha_S * b_t
##
## 따라서 true basis는:
##   b_t = (E[S_t] - E[S_1]) / (E[S_T] - E[S_1])
##
## ─────────────────────────────────────────────────────────────────

compute_true_basis <- function(lam_mat, nu_mat, kappa_I, kappa_S) {
  tlen <- ncol(lam_mat)
  tvec <- 0:(tlen - 1)
  
  s_t <- colSums(lam_mat)   # sum of loadings per time
  c_t <- colSums(nu_mat)    # sum of intercepts per time
  
  ## True expected sum score
  E_St <- c_t + s_t * (kappa_I + tvec * kappa_S)
  
  ## True basis (anchored at first=0, last=1)
  b_true <- (E_St - E_St[1]) / (E_St[tlen] - E_St[1])
  
  ## For comparison: linear basis would be
  b_linear <- tvec / (tlen - 1)
  
  list(
    E_St     = E_St,
    b_true   = b_true,
    b_linear = b_linear,
    s_t      = s_t,
    c_t      = c_t
  )
}


## ═══════════════════════════════════════════════════════════════════
## 사용 예시: 기존 시뮬레이션 루프 안에 삽입
## ═══════════════════════════════════════════════════════════════════

## --- 기존 코드에서 data.comp 생성 직후 (line ~669) ---
##
## [기존] fit <- lavaan::growth(C.strc.EC, data = data.comp)
##
## [대체] 아래와 같이 세 모형을 한꺼번에 fit
##
## comp_results <- fit_composite_models(data.comp)
##
## ## fit indices 비교
## print(comp_results$fit_indices)
##
## ## growth parameter 비교
## print(comp_results$growth_params)
##
## ## latent basis가 추정한 basis vs. true basis
## true_b <- compute_true_basis(lam_mat, nu_mat, k1, k2)
## cat("True basis:     ", round(true_b$b_true, 4), "\n")
## cat("Estimated basis:", round(comp_results$basis_estimates, 4), "\n")
## cat("Linear basis:   ", round(true_b$b_linear, 4), "\n")


## ═══════════════════════════════════════════════════════════════════
## Quick diagnostic test (단일 데이터셋으로 빠른 확인용)
## ═══════════════════════════════════════════════════════════════════

quick_test <- function() {
  
  cat("=== Quick diagnostic: 3 composite models ===\n\n")
  
  ## 간단한 DGP
  set.seed(42)
  I <- 5; tlen <- 5; nsize <- 500
  lambda <- 0.7; tau <- 1.0
  kI <- 3.0; kS <- 0.5
  phi1 <- 0.5; phi2 <- 0.2; phi12 <- -0.05
  
  phi.mtx <- matrix(c(phi1^2, phi12, phi12, phi2^2), 2, 2)
  
  ## Growth factors
  phi_draw <- MASS::mvrnorm(nsize, c(0, 0), phi.mtx)
  xi1 <- kI + phi_draw[, 1]
  xi2 <- kS + phi_draw[, 2]
  
  tvec <- 0:(tlen - 1)
  eta_mat <- outer(xi1, rep(1, tlen)) + outer(xi2, tvec)
  
  ## Non-invariant loadings (config, large)
  lam_base <- rep(lambda, I)
  lam_mat  <- matrix(lam_base, I, tlen)
  u <- tvec / (tlen - 1)
  drift <- u^2
  for (tt in seq_len(tlen)) {
    lam_mat[3:5, tt] <- lam_base[3:5] * (1 + 0.5 * drift[tt])
  }
  
  nu_mat <- matrix(tau, I, tlen)
  
  ## Generate observed data
  sigma_mu <- sqrt(1 - lambda^2)
  y_mat <- matrix(NA, nsize, I * tlen)
  for (tt in seq_len(tlen)) {
    for (ii in seq_len(I)) {
      col_idx <- (tt - 1) * I + ii
      y_mat[, col_idx] <- nu_mat[ii, tt] + lam_mat[ii, tt] * eta_mat[, tt] +
        rnorm(nsize, 0, sigma_mu)
    }
  }
  
  ## Sum scores
  data.comp <- data.frame(
    time1 = rowSums(y_mat[, 1:5]),
    time2 = rowSums(y_mat[, 6:10]),
    time3 = rowSums(y_mat[, 11:15]),
    time4 = rowSums(y_mat[, 16:20]),
    time5 = rowSums(y_mat[, 21:25])
  )
  
  ## Fit all three
  res <- fit_composite_models(data.comp)
  
  cat("── Fit Indices ──\n")
  print(res$fit_indices[, c("model","chi2","df","cfi","rmsea","aic","bic")],
        row.names = FALSE, digits = 4)
  
  cat("\n── Growth Parameters ──\n")
  print(res$growth_params, row.names = FALSE, digits = 4)
  
  cat("\n── Latent Basis Estimates ──\n")
  true_b <- compute_true_basis(lam_mat, nu_mat, kI, kS)
  cat("True basis:      ", round(true_b$b_true, 4), "\n")
  cat("Estimated basis: ", round(res$basis_estimates, 4), "\n")
  cat("Linear basis:    ", round(true_b$b_linear, 4), "\n")
  
  cat("\n── True E[S_t] ──\n")
  cat("E[S_t]:          ", round(true_b$E_St, 4), "\n")
  cat("Observed means:  ", round(colMeans(data.comp), 4), "\n")
  
  invisible(res)
}

## quick_test()  # 실행하려면 주석 해제
