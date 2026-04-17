run <- function(
    nrep = 2000,
    N_vec = c(200, 500, 800),
    model = c("latent", "composite", "scaled"),
    
    ## Design ------------------------------
    TT = 5,    # time points
    J = 5,    # number of indicators / time
    
    ## Growth Parameters ------------------
    mu_alpha  = 0,
    mu_beta = 0.5,
    var_alpha = 1,
    var_beta = 0.2,
    cov_ab = 0.1,
    
    ## Measurement Parameters -------------å
    target_loading = 0.4,   # Standardized target loading at t=1
    tau_un  = 1,            # baseline intercept (all items)
    lambda_un = 1,          # baseline loading (all items)
    sd.lambda = 0,          # loading fluctuation
    rho = 0,                # correlated uniqueness
    psi = 0.2,              # latent variable noise
    
    ## Measurement Non-Invariance ----------------
    ## location: which parameter drifts over time 
    ## size:     magnitude of drift × number of non-invariant items
    mni_location = c("none", "intercept", "loading", "residual"),
    mni_size     = c("none", "strong", "weak"),
    
    seed = 123
) {
  
  model        <- match.arg(model)
  mni_location <- match.arg(mni_location)
  mni_size     <- match.arg(mni_size)
  
  # 0. Fixed settings 
  ## naming
  col_names <- character(J * TT)
  for (t in seq_len(TT)) {
    for (j in seq_len(J)) {
      col_names[(t - 1L) * J + j] <- paste0("i", j, "t", t)
    }
  }
  

  # 1. MNI setup 
  mni_params <- switch(
    mni_size,
    none = list(n_noninv = 0L, d_tau = 0, d_lam = 0, d_res = 0),
    strong = list(n_noninv = 3L, d_tau = 0.30, d_lam = 0.30, d_res = 0.50),
    weak = list(n_noninv = 1L, d_tau = 0.10, d_lam = 0.10, d_res = 0.20)
  )
  
  ## which items are non-invariant (the last n items)
  noninv_items <- if (mni_params$n_noninv > 0L) {
    tail(seq_len(J), mni_params$n_noninv)
  } else {
    integer(0)
  }
  
  ## drift_t: Non-linear change
  t_idx    <- 0:(TT - 1)
  u        <- t_idx / (TT - 1)    # 0~1 interval
  drift_t  <- u^2                 # Non-linear drift 

  ## Parameters
  tau_vec <- rep(tau_un, times=J)
  tau_mat <- matrix(tau_un, nrow = J, ncol = TT)
  
  lam_vec <- truncnorm::rtruncnorm(
    J,
    a    = lambda_un - .2,
    b    = lambda_un + .2,
    mean = lambda_un,
    sd   = sd.lambda)
  lam_mat <- matrix(lam_vec, nrow = J, ncol = TT)

  var_eta   <- var_alpha + psi
  theta_vec <- (lam_vec)^2 * var_eta * (1 / target_loading^2 - 1)
  theta_mat <- matrix(rep(theta_vec, TT), nrow = J)

  psi_mat <- c(psi, psi+.1, psi+.2, psi+.3, psi+.4)
  
  if (length(noninv_items) > 0L) {
    if (mni_location == "intercept") {
      for (t in seq_len(TT)) {
        tau_mat[noninv_items, t] <- tau_vec[noninv_items] + mni_params$d_tau * drift_t[t]
      }
      
    } else if (mni_location == "loading") {
      for (t in seq_len(TT)) {
        lam_mat[noninv_items, t] <- lam_vec[noninv_items] + mni_params$d_lam * drift_t[t]
      }
      
    } else if (mni_location == "residual") {
      for (t in seq_len(TT)) {
        theta_mat[noninv_items, t] <- theta_vec[noninv_items] * (1 + mni_params$d_res * drift_t[t])
      }
    }
  }
  
  ## Residual covariance with correlated uniqueness
  make_resid_cov <- function(var_e, TT, rho) {
    R <- diag(TT)
    for (t in seq_len(TT - 1L)) {
      R[t, t + 1L] <- rho
      R[t + 1L, t] <- rho
    }
    D <- diag(sqrt(var_e), nrow = TT)
    D %*% R %*% D
  }
  
  # 2. Analysis
  if (model == "latent") {
    
    lav_model <- '
# Measurement part
f1 =~ 1*i1t1 + l21*i2t1 + l31*i3t1 + l41*i4t1 + l51*i5t1
f2 =~ 1*i1t2 + l22*i2t2 + l32*i3t2 + l42*i4t2 + l52*i5t2
f3 =~ 1*i1t3 + l23*i2t3 + l33*i3t3 + l43*i4t3 + l53*i5t3
f4 =~ 1*i1t4 + l24*i2t4 + l34*i3t4 + l44*i4t4 + l54*i5t4
f5 =~ 1*i1t5 + l25*i2t5 + l35*i3t5 + l45*i4t5 + l55*i5t5

# Intercepts
i1t1 ~ 1*1;   i1t2 ~ 1*1;   i1t3 ~ 1*1;   i1t4 ~ 1*1;   i1t5 ~ 1*1
i2t1 ~ int2*1; i2t2 ~ int2*1; i2t3 ~ int2*1; i2t4 ~ int2*1; i2t5 ~ int2*1
i3t1 ~ int3*1; i3t2 ~ int3*1; i3t3 ~ int3*1; i3t4 ~ int3*1; i3t5 ~ int3*1
i4t1 ~ int4*1; i4t2 ~ int4*1; i4t3 ~ int4*1; i4t4 ~ int4*1; i4t5 ~ int4*1
i5t1 ~ int5*1; i5t2 ~ int5*1; i5t3 ~ int5*1; i5t4 ~ int5*1; i5t5 ~ int5*1

# Residual variances
i1t1 ~~ e1*i1t1; i1t2 ~~ e1*i1t2; i1t3 ~~ e1*i1t3; i1t4 ~~ e1*i1t4; i1t5 ~~ e1*i1t5
i2t1 ~~ e2*i2t1; i2t2 ~~ e2*i2t2; i2t3 ~~ e2*i2t3; i2t4 ~~ e2*i2t4; i2t5 ~~ e2*i2t5
i3t1 ~~ e3*i3t1; i3t2 ~~ e3*i3t2; i3t3 ~~ e3*i3t3; i3t4 ~~ e3*i3t4; i3t5 ~~ e3*i3t5
i4t1 ~~ e4*i4t1; i4t2 ~~ e4*i4t2; i4t3 ~~ e4*i4t3; i4t4 ~~ e4*i4t4; i4t5 ~~ e4*i4t5
i5t1 ~~ e5*i5t1; i5t2 ~~ e5*i5t2; i5t3 ~~ e5*i5t3; i5t4 ~~ e5*i5t4; i5t5 ~~ e5*i5t5

# Correlated uniqueness
i1t1 ~~ i1t2; i1t2 ~~ i1t3; i1t3 ~~ i1t4; i1t4 ~~ i1t5
i2t1 ~~ i2t2; i2t2 ~~ i2t3; i2t3 ~~ i2t4; i2t4 ~~ i2t5
i3t1 ~~ i3t2; i3t2 ~~ i3t3; i3t3 ~~ i3t4; i3t4 ~~ i3t5
i4t1 ~~ i4t2; i4t2 ~~ i4t3; i4t3 ~~ i4t4; i4t4 ~~ i4t5
i5t1 ~~ i5t2; i5t2 ~~ i5t3; i5t3 ~~ i5t4; i5t4 ~~ i5t5

# Second-order LGCM
I =~ 1*f1 + 1*f2 + 1*f3 + 1*f4 + 1*f5
S =~ 0*f1 + 1*f2 + 2*f3 + 3*f4 + 4*f5

I ~~ (phi1)*I
S ~~ (phi2)*S
I ~~ (phi12)*S
I ~  (k1)*1
S ~  (k2)*1

# Factor intercepts fixed 0 
f1 ~ 0*1;  f2 ~ 0*1;  f3 ~ 0*1;  f4 ~ 0*1;  f5 ~ 0*1

# Latent residual variances
f1 ~~ p1*f1;  f2 ~~ p2*f2;  f3 ~~ p3*f3;  f4 ~~ p4*f4;  f5 ~~ p5*f5
'
  } else {
    
    lav_model <- '
i =~ 1*time1 + 1*time2 + 1*time3 + 1*time4 + 1*time5
s =~ 0*time1 + 1*time2 + 2*time3 + 3*time4 + 4*time5

i ~~ (phi1)*i
s ~~ (phi2)*s
i ~~ (phi12)*s
i ~  (k1)*1
s ~  (k2)*1
'
  }  
  # 3. Output names
  p_names <- c("phi11","phi11.SE","phi11.p",
               "phi22","phi22.SE","phi22.p",
               "phi12","phi12.SE","phi12.p",
               "kappa1","kappa1.SE","kappa1.p",
               "kappa2","kappa2.SE","kappa2.p",
               "err.mark",
               "phi11.cov","phi12.cov","phi22.cov","kappa1.cov","kappa2.cov",
               "chisq","df","pvalue.chi","cfi","tli","rmsea","srmr")
  
  # 4. For Loop
  zval <- 1.96
  
  ## CI calculation
  ci_in <- function(est, se, true) {
    as.numeric((est - zval * se) <= true & true <= (est + zval * se))
  }
  
  ## Relative bias
  rbias_fun <- function(est, true) {
    if (isTRUE(all.equal(true, 0))) return(NA_real_)
    mean((est - true) / true, na.rm = TRUE)
  }
  
  results_list <- vector("list", length(N_vec))
  names(results_list) <- paste0("N", N_vec)
  

  for (n_idx in seq_along(N_vec)) {
    N <- N_vec[n_idx]
    cat(sprintf("[N = %5d]  Simulation...\n", N))
    
    recall <- as.data.frame(matrix(NA_real_, nrow = nrep, ncol = length(p_names)))
    colnames(recall) <- p_names
    
    for (i in seq_len(nrep)) {
      
      ## seed per rep
      set.seed(seed + i)
      
      # (A) Data Generation 
      ## (1) Growth Factor
      phi.mtx <- matrix(c(var_alpha, cov_ab,
                           cov_ab,   var_beta), 2, 2, byrow = TRUE)
      phi_draw <- MASS::mvrnorm(N, mu = c(mu_alpha, mu_beta), Sigma = phi.mtx)
      w1 <- phi_draw[, 1]; w2 <- phi_draw[, 2]
      
      ## (2) Latent Factor Scores
      eta_mat <- matrix(NA_real_, N, TT)
      for (t in seq_len(TT)) {
        zeta_t <- rnorm(N, mean = 0, sd = sqrt(psi_mat[t]))
        eta_mat[, t] <- w1 + w2 * t_idx[t] + zeta_t
      }
      
      ## residual variance
      dat_mat <- matrix(NA_real_, N, J * TT)
      for (j in seq_len(J)) {
        Sigma_ej <- make_resid_cov(theta_mat[j, ], TT, rho = rho)
        E_j      <- MASS::mvrnorm(N, mu = rep(0, TT), Sigma = Sigma_ej)
        for (t in seq_len(TT)) {
          dat_mat[, (t - 1L) * J + j] <-
            tau_mat[j, t] + lam_mat[j, t] * eta_mat[, t] + E_j[, t]
        }
      }

      dat <- as.data.frame(dat_mat)
      colnames(dat) <- col_names
      
      # (B) Model Analysis 
      if (model == "composite") {
        dat_fit <- data.frame(
          time1 = rowMeans(dat[,  1: 5]),
          time2 = rowMeans(dat[,  6:10]),
          time3 = rowMeans(dat[, 11:15]),
          time4 = rowMeans(dat[, 16:20]),
          time5 = rowMeans(dat[, 21:25])
        )
      } else if (model == "scaled"){
        mean_t <- data.frame(
          time1 = rowSums(dat[,  1: 5]),
          time2 = rowSums(dat[,  6:10]),
          time3 = rowSums(dat[, 11:15]),
          time4 = rowSums(dat[, 16:20]),
          time5 = rowSums(dat[, 21:25])
        )
        
        # target latent metric: summed loadings and summed intercepts
        a_t <- colSums(lam_mat)    # = sum of loadings at each wave
        c_t <- colSums(tau_mat)    # = sum of intercepts at each wave
        
        # rescaled composite scores
        mean_t_lat <- sweep(mean_t, 2, c_t, "-")
        mean_t_lat <- sweep(mean_t_lat, 2, a_t, "/")
        
        dat_fit <- as.data.frame(mean_t_lat)
        names(dat_fit) <- paste0("time", 1:5)
      } else {
        dat_fit <- dat
        }
      
      fit <- tryCatch({
        if (model == "latent") {
          lavaan::growth(
            lav_model,
            data = dat_fit,
            meanstructure = TRUE,
            control = list(reltol = 1e-12, xtol_rel = 1e-12)
          )
        } else {
          lavaan::growth(lav_model, data = dat_fit)
        }
      }, error = function(e) NULL)
      
      if (is.null(fit)) {
        recall[i, "err.mark"] <- 1L
        next
      }
      
      pe <- as.data.frame(lavaan::parameterEstimates(fit, standardized = TRUE))
      pe <- pe[pe$label != "", c("op","label","est","se","pvalue","std.all")]
      
      pe$err.mark <- dplyr::case_when(
        Reduce(`|`, lapply(pe, is.na))                    ~ 1,
        pe$op == "~~" & pe$std.all >=  1.000001          ~ 1,
        pe$op == "~~" & pe$std.all <= -1.000001          ~ 1,
        pe$op == "~~" & pe$label != "phi12" & pe$est < 0 ~ 1,
        TRUE                                             ~ NA_real_
      )
      
      gv <- function(lbl, col) {
        v <- pe[pe$label == lbl, col, drop = TRUE]
        if (length(v) == 0L) NA_real_ else v[1L]
      }
      
      ## Save output
      recall[i,  1L] <- gv("phi1",  "est")
      recall[i,  2L] <- gv("phi1",  "se")
      recall[i,  3L] <- gv("phi1",  "pvalue")
      recall[i,  4L] <- gv("phi2",  "est")
      recall[i,  5L] <- gv("phi2",  "se")
      recall[i,  6L] <- gv("phi2",  "pvalue")
      recall[i,  7L] <- gv("phi12", "est")
      recall[i,  8L] <- gv("phi12", "se")
      recall[i,  9L] <- gv("phi12", "pvalue")
      recall[i, 10L] <- gv("k1",    "est")
      recall[i, 11L] <- gv("k1",    "se")
      recall[i, 12L] <- gv("k1",    "pvalue")
      recall[i, 13L] <- gv("k2",    "est")
      recall[i, 14L] <- gv("k2",    "se")
      recall[i, 15L] <- gv("k2",    "pvalue")
      
      ## Non-convergence flag
      recall[i, 16L] <- ifelse(sum(pe$err.mark, na.rm = TRUE) >= 1L, 1L, NA_real_)
      
      ## 95% Coverage against generating values
      recall[i, 17L] <- ci_in(recall[i,  1L], recall[i,  2L], var_alpha)
      recall[i, 18L] <- ci_in(recall[i,  7L], recall[i,  8L], cov_ab)
      recall[i, 19L] <- ci_in(recall[i,  4L], recall[i,  5L], var_beta)
      recall[i, 20L] <- ci_in(recall[i, 10L], recall[i, 11L], mu_alpha)
      recall[i, 21L] <- ci_in(recall[i, 13L], recall[i, 14L], mu_beta)
      
      fi <- tryCatch(
        lavaan::fitMeasures(fit, c("chisq","df","pvalue","cfi","tli","rmsea","srmr")),
        error = function(e) rep(NA_real_, 7)
      )
      recall[i, 22L] <- fi[1L]   #chisq
      recall[i, 23L] <- fi[2L]   #df
      recall[i, 24L] <- fi[3L]   #chisq p-value
      recall[i, 25L] <- fi[4L]   #cfi
      recall[i, 26L] <- fi[5L]   #tli
      recall[i, 27L] <- fi[6L]   #rmsea
      recall[i, 28L] <- fi[7L]   #srmr
      
    } # end nrep
    
    # (C) Summary 
    err.perc_ <- mean(!is.na(recall$err.mark))
    recall_ok <- dplyr::filter(recall, is.na(err.mark))
    
    summary_row <- data.frame(
      N = N,
      nrep_ok  = nrow(recall_ok),
      err.perc = err.perc_,
      
      phi11.m = mean(recall_ok$phi11, na.rm = TRUE),
      phi22.m = mean(recall_ok$phi22, na.rm = TRUE),
      phi12.m = mean(recall_ok$phi12, na.rm = TRUE),
      k1.m    = mean(recall_ok$kappa1, na.rm = TRUE),
      k2.m    = mean(recall_ok$kappa2, na.rm = TRUE),
      
      phi11.bias = mean(recall_ok$phi11 - var_alpha, na.rm = TRUE),
      phi22.bias = mean(recall_ok$phi22 - var_beta,  na.rm = TRUE),
      phi12.bias = mean(recall_ok$phi12 - cov_ab,    na.rm = TRUE),
      k1.bias    = mean(recall_ok$kappa1 - mu_alpha, na.rm = TRUE),
      k2.bias    = mean(recall_ok$kappa2 - mu_beta,  na.rm = TRUE),
      
      phi11.rbias = rbias_fun(recall_ok$phi11,  var_alpha),
      phi22.rbias = rbias_fun(recall_ok$phi22,  var_beta),
      phi12.rbias = rbias_fun(recall_ok$phi12,  cov_ab),
      k1.rbias    = rbias_fun(recall_ok$kappa1, mu_alpha),
      k2.rbias    = rbias_fun(recall_ok$kappa2, mu_beta),
      
      phi11.rmse = sqrt(mean((recall_ok$phi11 - var_alpha)^2, na.rm = TRUE)),
      phi22.rmse = sqrt(mean((recall_ok$phi22 - var_beta)^2,  na.rm = TRUE)),
      phi12.rmse = sqrt(mean((recall_ok$phi12 - cov_ab)^2,    na.rm = TRUE)),
      k1.rmse    = sqrt(mean((recall_ok$kappa1 - mu_alpha)^2, na.rm = TRUE)),
      k2.rmse    = sqrt(mean((recall_ok$kappa2 - mu_beta)^2,  na.rm = TRUE)),
      
      phi11.cov = mean(recall_ok$phi11.cov, na.rm = TRUE),
      phi22.cov = mean(recall_ok$phi22.cov, na.rm = TRUE),
      phi12.cov = mean(recall_ok$phi12.cov, na.rm = TRUE),
      k1.cov = mean(recall_ok$kappa1.cov, na.rm = TRUE),
      k2.cov = mean(recall_ok$kappa2.cov, na.rm = TRUE),
      
      phi11.power = mean(recall_ok$phi11.p < 0.05, na.rm = TRUE),
      phi22.power = mean(recall_ok$phi22.p < 0.05, na.rm = TRUE),
      phi12.power = mean(recall_ok$phi12.p < 0.05, na.rm = TRUE),
      k1.power = mean(recall_ok$kappa1.p < 0.05, na.rm = TRUE),
      k2.power = mean(recall_ok$kappa2.p < 0.05, na.rm = TRUE),
      
      chisq.m = mean(recall_ok$chisq, na.rm = TRUE),
      df.m  = mean(recall_ok$df, na.rm = TRUE),
      pvalue.chi.m = mean(recall_ok$pvalue.chi, na.rm = TRUE),
      cfi.m = mean(recall_ok$cfi, na.rm = TRUE),
      tli.m = mean(recall_ok$tli, na.rm = TRUE),
      rmsea.m = mean(recall_ok$rmsea, na.rm = TRUE),
      srmr.m  = mean(recall_ok$srmr, na.rm = TRUE),
      
      chisq.sd = sd(recall_ok$chisq, na.rm = TRUE),
      cfi.sd  = sd(recall_ok$cfi, na.rm = TRUE),
      tli.sd  = sd(recall_ok$tli, na.rm = TRUE),
      rmsea.sd  = sd(recall_ok$rmsea, na.rm = TRUE),
      srmr.sd = sd(recall_ok$srmr, na.rm = TRUE),
      chisq.reject = mean(recall_ok$pvalue.chi < 0.05, na.rm = TRUE)
    )
    
    results_list[[n_idx]] <- list(
      raw     = recall_ok,
      summary = summary_row
    )
    
    cat(sprintf("  -> Complete (Converged %d/%d, Non-converged %.1f%%)\n",
                nrow(recall_ok), nrep, err.perc_ * 100))
  }
  
  # 5. Aggregate 
  summary_table <- do.call(rbind, lapply(results_list, `[[`, "summary"))
  rownames(summary_table) <- NULL
  
  list(
    by_N = results_list,
    summary_table = summary_table,
    true_values = list(
      phi1 = var_alpha,
      phi2 = var_beta,
      phi12 = cov_ab,
      k1 = mu_alpha,
      k2 = mu_beta
    ),
    conditions = list(
      mni_location = mni_location,
      mni_size = mni_size,
      noninv_items = noninv_items,
      tau_mat = tau_mat,
      lam_mat = lam_mat,
      theta_mat = theta_mat
    )
  )
}
