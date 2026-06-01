library(dplyr); library(lavaan); library(MASS);

# Admissible Results 
.lav_admissible <- function(fit) {
  if (is.null(fit) || !inherits(fit, "lavaan")) return(FALSE)   # lavaan object?
  if (!isTRUE(lavaan::lavInspect(fit, "converged"))) return(FALSE)  # converged?
  if (!isTRUE(lavaan::lavTech(fit, "post.check"))) return(FALSE)  # any neg.var?
  pe <- tryCatch(
    as.data.frame(lavaan::parameterEstimates(fit, standardized = TRUE)),
    error = function(e) NULL)
  if (is.null(pe) || nrow(pe) == 0L) return(FALSE)  # no output
  if (any(is.na(pe$est))) return(FALSE)  # no estimates
  var_est <- pe$op == "~~" & pe$lhs == pe$rhs   # variance
  if (any(var_est & pe$est < 0, na.rm = TRUE)) return(FALSE)  # negative variance
  cov_est <- pe$op == "~~" & pe$lhs != pe$rhs   # covariance
  if (any(cov_est & abs(pe$std.all) > 1.000001, na.rm = TRUE)) return(FALSE)  # correlation >1, <-1
  TRUE
}
  
  
run <- function(
      nrep = 2000,
      Nsize = c(200, 500, 800),
      model = c("EC_latent", "sum"),
      
      ## Design ----
      TT = c(3,5,8),    # time points
      J = c(3,5,8),     # number of indicators / time
      
      ## Growth Parameters ----
      mu_alpha  = 1,
      mu_beta   = 0.125,
      var_alpha = 1,
      var_beta  = 0.2025,
      cov_ab    = 0.18,
      cor_ab    = cov_ab / sqrt(var_alpha * var_beta),
      
      ## Measurement Parameters ---- 
      target_loading = c(.6708, .4),  # measurement quality (R^2 = .45 / .16)
      tau_un    = .1,
      lambda_un = .2,
      sd.lambda = c(0, 0.05),      # tau-equivalent, congeneric
      rho = c(0, 0.15),            # correlated uniqueness
      fit_cu = FALSE,
      psi_mat = c(1.097^2, 0.98^2, 0.95^2, 0.93^2, 0.91^2, 0.89^2, 0.87^2, 0.85^2),  # latent disturbances
  
      seed = 123
  ) {
    
    model <- match.arg(model)
  
    col_names <- character(J * TT)
    for (t in seq_len(TT)) {
      for (j in seq_len(J)) {
        col_names[(t - 1L) * J + j] <- paste0("i", j, "t", t)
      }
    }
    
    ##### 1. MNI setup #####
  
    ## Parameters
    phi.mtx <- matrix(c(var_alpha, cov_ab,
                        cov_ab,   var_beta), 2, 2, byrow = TRUE)
    tau_mat <- matrix(tau_un, nrow = J, ncol = TT)
    
    lam_vec <- if (sd.lambda == 0) {
      rep(lambda_un, J)
    } else {
      # sd.lambda = target SD of lambdas across items
      intervals <- seq(-(J - 1) / 2, (J - 1) / 2, by = 1)  
      ratio   <- sd.lambda / sd(intervals)
      lambda_un + intervals * ratio
    }
    lam_mat <- matrix(lam_vec, nrow = J, ncol = TT)
    
    var_eta   <- var_alpha + psi_mat[1]          # targetted at t=1
    theta_vec <- (lam_vec)^2 * var_eta * (1 / target_loading^2 - 1)
    theta_mat <- matrix(rep(theta_vec, TT), nrow = J)
    
    ## Residual covariance with correlated uniqueness (AR Lag-1)
    make_resid_cov <- function(var_e, TT, rho) {
      R <- diag(TT)
      for (t in seq_len(TT - 1L)) {
        R[t, t + 1L] <- rho     # AR Lag-1 only
        R[t + 1L, t] <- rho
      }
      D <- diag(sqrt(var_e), nrow = TT)
      D %*% R %*% D             # cov = sd * cor * sd
    }
    
    ##### 2. Lavaan model syntax ######
    
    # For easy coding 
    load_t <- function(t, J)   # measurement
      paste0("f", t, " =~ ", paste(sprintf("lmb%d*i%dt%d", 1:J, 1:J, t), collapse = " + "))
    
    int_t  <- function(t, J)   # intercepts
      paste(sprintf("i%dt%d ~ nu%d*1", 1:J, t, 1:J), collapse = "; ")
    
    res_j  <- function(j, TT)  # residuals
      paste(sprintf("i%dt%d ~~ e%d*i%dt%d", j, 1:TT, j, j, 1:TT), collapse = "; ")
    
    cu_j   <- function(j, TT)  # CU (AR lag-1)
      paste(sprintf("i%dt%d ~~ i%dt%d", j, 1:(TT-1), j,  2:TT), collapse = "; ")
    
    
    ## 2.1. EC-fitted 2-LGCM #####
    if (model == "EC_latent") {
      
      # Structural 
      I_eq  <- paste0("I =~ ",  paste(sprintf("1*f%d", 1:TT), collapse = " + "))
      S_eq  <- paste0("S =~ ",  paste(sprintf("%d*f%d", 0:(TT-1), 1:TT), collapse = " + "))
      f_int <- paste(sprintf("f%d ~ 0*1", 1:TT), collapse = "; ")
      f_var <- paste(sprintf("f%d ~~ s%d*f%d", 1:TT, 1:TT, 1:TT), collapse = "; ")
      
      # Measurement / Intercepts / Residuals
      meas  <- paste(sapply(1:TT, load_t, J = J),    collapse = "\n")
      ints  <- paste(sapply(1:TT, int_t,  J = J),    collapse = "\n")
      rvars <- paste(sapply(1:J,  res_j,  TT = TT),  collapse = "\n")
      
      # EC constraints
      ec_con <- paste0(
        paste(sprintf("lmb%d", 1:J), collapse = " + "), " == 1\n",
        paste(sprintf("nu%d",  1:J), collapse = " + "), " == 0"
      )
      
      lav_model <- paste(
        "# Second-order LGCM (effect coding)",
        I_eq, S_eq,
        "I ~~ (phi1)*I",
        "S ~~ (phi2)*S",
        "I ~~ (phi12)*S",
        "I ~ (k1)*1",
        "S ~ (k2)*1",
        f_int, f_var,
        "# Measurement",  meas,
        "# Intercepts",   ints,
        "# Residual variances", rvars,
        "# EC constraints", ec_con,
        sep = "\n"
      )
      
      if (fit_cu && TT != 3) { 
    cu_items  <- paste(sapply(1:J, cu_j, TT = TT), collapse = "\n")
    lav_model <- paste(lav_model, cu_items, sep = "\n")
  }
      
    } else {
      
      ## 2.2. Sum-score 1-LGCM #####
      i_eq  <- paste0("i =~ ", paste(sprintf("1*time%d", 1:TT), collapse = " + "))
      s_eq  <- paste0("s =~ ", paste(sprintf("%d*time%d", 0:(TT-1), 1:TT), collapse = " + "))
      t_var <- paste(sprintf("time%d ~~ s%d*time%d", 1:TT, 1:TT, 1:TT), collapse = "; ")
      
      lav_model <- paste(
        i_eq, s_eq,
        "i ~~ (phi1)*i",
        "s ~~ (phi2)*s",
        "i ~~ (phi12)*s",
        "i ~ (k1)*1",
        "s ~ (k2)*1",
        t_var,
        sep = "\n"
      )
      
      if (fit_cu && TT !=3) {
        cu_sums   <- paste(sprintf("time%d ~~ time%d", 1:(TT-1), 2:TT), collapse = "; ")
        lav_model <- paste(lav_model, cu_sums, sep = "\n")
      }
    }
    
    
    ##### 3. Saving Output #####
    p_base <- c("phi11","phi11.SE","phi11.p",
                "phi22","phi22.SE","phi22.p",
                "phi12","phi12.SE","phi12.p",
                "cor12", "cor12.SE", "cor12.p",
                "kappa1","kappa1.SE","kappa1.p",
                "kappa2","kappa2.SE","kappa2.p",
                "err.mark",
                "phi11.cov","phi12.cov","phi22.cov","cor12.cov","kappa1.cov","kappa2.cov",
                "chisq","df","pvalue.chi","cfi","tli","rmsea","srmr", "sigma2_bar")
    
    ##### 4. EC-fitted population under MI #####
    ec_means <- (lambda_un*J) %*% c(mu_alpha, mu_beta) + c(tau_un*J, 0) 
    ec_covs  <- (lambda_un*J)^2 * phi.mtx
    
    k1_true    <- ec_means[1L]
    k2_true    <- ec_means[2L]
    phi1_true  <- ec_covs[1L,1L]
    phi12_true <- ec_covs[2L,1L]
    phi2_true  <- ec_covs[2L,2L]
    
    
    ##### 5. Simulation loop #####
    zval <- 1.96
    t_idx  <- 0:(TT - 1)
    resid_labs <- paste0("s", seq_len(TT))
    
    ci_in <- function(est, se, true) {     # confidence interval
      as.numeric((est - zval * se) <= true & true <= (est + zval * se))
    }
    
    rbias_fun <- function(est, true) {     # relative bias
      if (isTRUE(all.equal(true, 0))) return(NA_real_)
      mean((est - true) / true, na.rm = TRUE)
    }
    
    results_list <- vector("list", length(Nsize))
    names(results_list) <- paste0("N", Nsize)
    
    
    for (n_idx in seq_along(Nsize)) {
      N <- Nsize[n_idx]
      
      recall <- as.data.frame(matrix(NA_real_, nrow = nrep, ncol = length(p_base)))
      colnames(recall) <- p_base
      
      for (i in seq_len(nrep)) {
        
        ## seed per rep
        set.seed(seed + i)
        
        # (A) Data Generation 
        ## (A.1) Growth Factor
        phi_draw <- MASS::mvrnorm(N, mu = c(mu_alpha, mu_beta), Sigma = phi.mtx)
        w1 <- phi_draw[, 1]; w2 <- phi_draw[, 2]
        
        ## (A.2) Latent Factor Scores
        eta_mat <- matrix(NA_real_, N, TT)
        for (t in seq_len(TT)) {
          zeta_t <- rnorm(N, mean = 0, sd = sqrt(psi_mat[t]))
          eta_mat[, t] <- w1 + w2 * t_idx[t] + zeta_t
        }
        
        ## (A.3) Observed data
        dat_mat <- matrix(NA_real_, N, J * TT)
        for (j in seq_len(J)) {
          Sigma_ej <- make_resid_cov(theta_mat[j, ], TT, rho = rho)
          E_j <- MASS::mvrnorm(N, mu = rep(0, TT), Sigma = Sigma_ej)
          for (t in seq_len(TT)) {
            dat_mat[, (t - 1L) * J + j] <-
              tau_mat[j, t] + lam_mat[j, t] * eta_mat[, t] + E_j[, t]
          }
        }
        
        dat <- as.data.frame(dat_mat)
        colnames(dat) <- col_names
        
        
        # (B) Data preparation
        if (model == "sum") {      # 1-LGCM
          dat_fit <- as.data.frame(
            sapply(seq_len(TT), function(t) {
              rowSums(dat[, ((t - 1L) * J + 1L):(t * J), drop = FALSE])
            })
          )
          colnames(dat_fit) <- paste0("time", seq_len(TT))
        } else {                   # 2-LGCM
          dat_fit <- dat
        }
        
        
        # (C) Model fitting
        fit <- tryCatch({
          if (model == "EC_latent") {
            lavaan::sem(lav_model,
                        data = dat_fit,
                        auto.fix.first = FALSE, 
                        auto.fix.single = FALSE,
                        effect.coding   = FALSE,
                        meanstructure = TRUE,
                        start = "simple",
                        control = list(reltol = 1e-10, xtol_rel = 1e-10, iter.max = 5000))
          } else {
            lavaan::growth(lav_model, data = dat_fit,
                           control = list(reltol = 1e-10, xtol_rel = 1e-10, iter.max = 5000))
          }
        }, error = function(e) NULL)
        
        if (!.lav_admissible(fit)) {
          recall[i, "err.mark"] <- 1L
          next  # skip the left (-> wrong/NA estimates are excluded in results calculation)
        }
        
        recall[i, "err.mark"] <- 0L
        fit_example <- fit
        
        pe <- as.data.frame(lavaan::parameterEstimates(fit, standardized = TRUE))
        pe <- pe[pe$label != "", c("op","label","est","se","pvalue","std.all")]
        
        # sigma2_bar: average residual variance 
        recall[i, "sigma2_bar"] <- mean(pe$est[pe$label %in% resid_labs], na.rm = TRUE)
        
        gv <- function(lbl, col) {
          v <- pe[pe$label == lbl, col, drop = TRUE]
          if (length(v) == 0L) NA_real_ else v[1L]
        }
        
        ss <- as.data.frame(lavaan::standardizedSolution(fit))
        gs <- function(lbl, col) {
          v <- ss[ss$label == lbl, col, drop = TRUE]
          if (length(v) == 0L) NA_real_ else v[1L]
        }
        
        # (D) Extract common parameters
        recall[i, "phi11"]    <- gv("phi1",  "est")
        recall[i, "phi11.SE"] <- gv("phi1",  "se")
        recall[i, "phi11.p"]  <- gv("phi1",  "pvalue")
        recall[i, "phi22"]    <- gv("phi2",  "est")
        recall[i, "phi22.SE"] <- gv("phi2",  "se")
        recall[i, "phi22.p"]  <- gv("phi2",  "pvalue")
        recall[i, "phi12"]    <- gv("phi12", "est")
        recall[i, "phi12.SE"] <- gv("phi12", "se")
        recall[i, "phi12.p"]  <- gv("phi12", "pvalue")
        recall[i, "cor12"]    <- gs("phi12", "est.std")
        recall[i, "cor12.SE"] <- gs("phi12", "se")
        recall[i, "cor12.p"]   <- gs("phi12", "pvalue")
        recall[i, "kappa1"]   <- gv("k1",  "est")
        recall[i, "kappa1.SE"]<- gv("k1",  "se")
        recall[i, "kappa1.p"] <- gv("k1",  "pvalue")
        recall[i, "kappa2"]   <- gv("k2",  "est")
        recall[i, "kappa2.SE"]<- gv("k2",  "se")
        recall[i, "kappa2.p"] <- gv("k2",  "pvalue")
        
        ## 95% Coverage against generating values
        recall[i, "phi11.cov"] <- ci_in(recall[i, "phi11"], recall[i, "phi11.SE"], phi1_true)
        recall[i, "phi12.cov"] <- ci_in(recall[i, "phi12"], recall[i, "phi12.SE"], phi12_true)
        recall[i, "phi22.cov"] <- ci_in(recall[i, "phi22"], recall[i, "phi22.SE"], phi2_true)
        recall[i, "cor12.cov"] <- ci_in(recall[i, "cor12"],  recall[i, "cor12.SE"], cor_ab)
        recall[i, "kappa1.cov"] <- ci_in(recall[i, "kappa1"], recall[i, "kappa1.SE"], k1_true)
        recall[i, "kappa2.cov"] <- ci_in(recall[i, "kappa2"], recall[i, "kappa2.SE"], k2_true)
        
        ## Fit indices
        fi <- tryCatch(
          lavaan::fitMeasures(fit, c("chisq","df","pvalue","cfi","tli","rmsea","srmr")),
          error = function(e) rep(NA_real_, 7)
        )
        recall[i, "chisq"]      <- fi[1L]
        recall[i, "df"]         <- fi[2L]
        recall[i, "pvalue.chi"] <- fi[3L]
        recall[i, "cfi"]        <- fi[4L]
        recall[i, "tli"]        <- fi[5L]
        recall[i, "rmsea"]      <- fi[6L]
        recall[i, "srmr"]       <- fi[7L]
      } # end nrep loop
      
      ##### 6. Summary #####
      err.perc_ <- mean(recall$err.mark == 1, na.rm = TRUE)
      recall_ok <- dplyr::filter(recall, err.mark == 0)
      
      ## Common summary 
      summary_row <- data.frame(
        model = model,
        N = N,
        nrep_ok  = nrow(recall_ok),
        err.perc = err.perc_,
        
        ## Mean estimates
        phi11.m = mean(recall_ok$phi11, na.rm = TRUE),
        phi22.m = mean(recall_ok$phi22, na.rm = TRUE),
        phi12.m = mean(recall_ok$phi12, na.rm = TRUE),
        cor12.m = mean(recall_ok$cor12, na.rm = TRUE),
        k1.m    = mean(recall_ok$kappa1, na.rm = TRUE),
        k2.m    = mean(recall_ok$kappa2, na.rm = TRUE),
        
        ## Bias
        phi11.bias = mean(recall_ok$phi11 - phi1_true,  na.rm = TRUE),
        phi22.bias = mean(recall_ok$phi22 - phi2_true,  na.rm = TRUE),
        phi12.bias = mean(recall_ok$phi12 - phi12_true, na.rm = TRUE),
        cor12.bias = mean(recall_ok$cor12 - cor_ab,     na.rm = TRUE), 
        k1.bias    = mean(recall_ok$kappa1 - k1_true,   na.rm = TRUE),
        k2.bias    = mean(recall_ok$kappa2 - k2_true,   na.rm = TRUE),
        
        ## Relative bias
        phi11.rbias = rbias_fun(recall_ok$phi11, phi1_true),
        phi22.rbias = rbias_fun(recall_ok$phi22, phi2_true),
        phi12.rbias = rbias_fun(recall_ok$phi12, phi12_true),
        cor12.rbias = rbias_fun(recall_ok$cor12, cor_ab),
        k1.rbias    = rbias_fun(recall_ok$kappa1, k1_true),
        k2.rbias    = rbias_fun(recall_ok$kappa2, k2_true),
        
        ## RMSE
        phi11.rmse = sqrt(mean((recall_ok$phi11 - phi1_true)^2,  na.rm = TRUE)),
        phi22.rmse = sqrt(mean((recall_ok$phi22 - phi2_true)^2,  na.rm = TRUE)),
        phi12.rmse = sqrt(mean((recall_ok$phi12 - phi12_true)^2, na.rm = TRUE)),
        cor12.rmse = sqrt(mean((recall_ok$cor12 - cor_ab)^2,     na.rm = TRUE)),
        k1.rmse    = sqrt(mean((recall_ok$kappa1 - k1_true)^2,   na.rm = TRUE)),
        k2.rmse    = sqrt(mean((recall_ok$kappa2 - k2_true)^2,   na.rm = TRUE)),
        
        ## Relative RMSE
        phi11.rel.rmse  = sqrt(mean((recall_ok$phi11 - phi1_true)^2, na.rm = TRUE)) / phi1_true,
        phi12.rel.rmse  = sqrt(mean((recall_ok$phi12 - phi12_true)^2, na.rm = TRUE)) / phi12_true,
        phi22.rel.rmse  = sqrt(mean((recall_ok$phi22 - phi2_true)^2, na.rm = TRUE)) / phi2_true,
        kappa1.rel.rmse = sqrt(mean((recall_ok$kappa1 - k1_true)^2, na.rm = TRUE)) / k1_true,
        kappa2.rel.rmse = sqrt(mean((recall_ok$kappa2 - k2_true)^2, na.rm = TRUE)) / k2_true,
        
        ## Mean Coverage 
        phi11.cov = mean(recall_ok$phi11.cov, na.rm = TRUE),
        phi22.cov = mean(recall_ok$phi22.cov, na.rm = TRUE),
        phi12.cov = mean(recall_ok$phi12.cov, na.rm = TRUE),
        cor12.cov = mean(recall_ok$cor12.cov, na.rm = TRUE),
        k1.cov = mean(recall_ok$kappa1.cov, na.rm = TRUE),
        k2.cov = mean(recall_ok$kappa2.cov, na.rm = TRUE),
        
        ## Power
        phi11.power = mean(recall_ok$phi11.p < 0.05, na.rm = TRUE),
        phi22.power = mean(recall_ok$phi22.p < 0.05, na.rm = TRUE),
        phi12.power = mean(recall_ok$phi12.p < 0.05, na.rm = TRUE),
        cor12.power = mean(recall_ok$cor12.p < 0.05, na.rm = TRUE),
        k1.power = mean(recall_ok$kappa1.p < 0.05, na.rm = TRUE),
        k2.power = mean(recall_ok$kappa2.p < 0.05, na.rm = TRUE),
    
        ## Fit indices
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
    } # end Nsize loop
    
    
    ##### 7. Aggregate the results #####
    summary_table <- do.call(rbind, lapply(results_list, `[[`, "summary"))
    rownames(summary_table) <- NULL
    
    list(
      by_N = results_list,
      summary_table = summary_table,
      true_values = list(
        phi1  = phi1_true,
        phi2  = phi2_true,
        phi12 = phi12_true,
        cor12 = cor_ab,
        k1    = k1_true,
        k2    = k2_true
      ), 
      lav_out = if (exists("fit_example", inherits = FALSE))    # save the latest admissible output
        lavaan::summary(fit_example, standardized = TRUE)
      else NULL,
      conditions = list(
        model = model,
        tau_mat = tau_mat,
        lam_mat = lam_mat,
        theta_mat = theta_mat
      )
    )
  }
