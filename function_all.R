run <- function(
    nrep = 2000,
    N_vec = c(200, 500, 800),
    model = c("latent", "mean", "sum", "mean_scaled", "sum_scaled", 
              "quadratic", "latent_basis"),
    
    ## Design ----
    TT = 5,    # time points
    J = 5,    # number of indicators / time
    
    ## Growth Parameters ----
    mu_alpha  = 0,
    mu_beta = 0.5,
    var_alpha = 1,
    var_beta = 0.2,
    cov_ab = 0.1,
    cor_ab = cov_ab / sqrt(var_alpha * var_beta),
    
    ## Measurement Parameters ----
    target_loading = 0.4,   # Standardized target loading at t=1
    tau_un  = 1,            
    lambda_un = 1,          
    sd.lambda = 0,          # loading fluctuation
    rho = 0,                # correlated uniqueness
    psi = 0.2,              # latent variable noise
    
    ## Measurement Non-Invariance ----
    ## location: which parameter drifts 
    ## size:     magnitude of drift × number of non-invariant items
    mni_location = c("none", "intercept", "loading", "residual"),
    mni_size     = c("none", "strong", "weak"),
    
    seed = 123
) {
  
  model <- match.arg(model)
  mni_location <- match.arg(mni_location)
  mni_size <- match.arg(mni_size)
  
  col_names <- character(J * TT)
  for (t in seq_len(TT)) {
    for (j in seq_len(J)) {
      col_names[(t - 1L) * J + j] <- paste0("i", j, "t", t)
    }
  }
  
  # 1. MNI setup 
  mni_params <- switch(
    mni_size,
    none = list(n_noninv = 0, d_tau = 0, d_lam = 0, d_res = 0),
    weak = list(n_noninv = 1, d_tau = 0.10, d_lam = 0.10, d_res = 0.20),
    strong = list(n_noninv = 3, d_tau = 0.30, d_lam = 0.30, d_res = 0.50)
  )
  
  ## the last n Non-Invariant items
  noninv_items <- if (mni_params$n_noninv > 0L) {
    tail(seq_len(J), mni_params$n_noninv)
  } else {
    integer(0)
  }
  
  ## drift properties: Non-linear change
  t_idx  <- 0:(TT - 1)
  u  <- t_idx / (TT - 1)    # intervals in 0~1 
  drift_t  <- u^2                 # Quadratic drift 
  
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
  
  var_eta   <- var_alpha + psi    # targetted at t=1
  theta_vec <- (lam_vec)^2 * var_eta * (1 / target_loading^2 - 1)
  theta_mat <- matrix(rep(theta_vec, TT), nrow = J)
  
  psi_mat <- c(psi, psi+.1, psi+.2, psi+.3, psi+.4)
  
  
  ## Measurement Non-Invariance
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
    D %*% R %*% D  # cov = sd * cor * sd
  }
  
  # ================================================================
  # 2. Lavaan model syntax
  # ================================================================
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
i2t1 ~ int21*1; i2t2 ~ int22*1; i2t3 ~ int23*1; i2t4 ~ int24*1; i2t5 ~ int25*1
i3t1 ~ int31*1; i3t2 ~ int32*1; i3t3 ~ int33*1; i3t4 ~ int34*1; i3t5 ~ int35*1
i4t1 ~ int41*1; i4t2 ~ int42*1; i4t3 ~ int43*1; i4t4 ~ int44*1; i4t5 ~ int45*1
i5t1 ~ int51*1; i5t2 ~ int52*1; i5t3 ~ int53*1; i5t4 ~ int54*1; i5t5 ~ int55*1

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
    
  } else if (model == "quadratic") {
    lav_model <- '
i =~ 1*time1 + 1*time2 + 1*time3 + 1*time4 + 1*time5
s =~ 0*time1 + 1*time2 + 2*time3 + 3*time4 + 4*time5
q =~ 0*time1 + 1*time2 + 4*time3 + 9*time4 + 16*time5

i ~~ (phi1)*i
s ~~ (phi2)*s
q ~~ (phi3)*q
i ~~ (phi12)*s
i ~~ (phi13)*q
s ~~ (phi23)*q

i ~ (k1)*1
s ~ (k2)*1
q ~ (k3)*1
'
    
  } else if (model == "latent_basis") {
    lav_model <- '
i =~ 1*time1 + 1*time2 + 1*time3 + 1*time4 + 1*time5
s =~ 0*time1 + (b2)*time2 + (b3)*time3 + (b4)*time4 + 1*time5

i ~~ (phi1)*i
s ~~ (phi2)*s
i ~~ (phi12)*s

i ~ (k1)*1
s ~ (k2)*1
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
  
  # ================================================================
  # 3. Output column names
  # ================================================================
  p_base <- c("phi11","phi11.SE","phi11.p",
              "phi22","phi22.SE","phi22.p",
              "phi12","phi12.SE","phi12.p",
              "cor12", "cor12.SE", "cor12.p",
              "kappa1","kappa1.SE","kappa1.p",
              "kappa2","kappa2.SE","kappa2.p",
              "err.mark",
              "phi11.cov","phi12.cov","phi22.cov","cor12.cov","kappa1.cov","kappa2.cov",
              "chisq","df","pvalue.chi","cfi","tli","rmsea","srmr")
  
  ## quadratic 
  p_quad <- c("kappa3","kappa3.SE","kappa3.p",
              "phi33","phi33.SE","phi33.p",
              "phi13","phi13.SE","phi13.p",
              "phi23","phi23.SE","phi23.p",
              "kappa3.cov")
  
  ## latent basis 
  p_lb <- c("b2","b2.SE","b3","b3.SE","b4","b4.SE")
  
  if (model == "quadratic") {
    p_names <- c(p_base, p_quad)
  } else if (model == "latent_basis") {
    p_names <- c(p_base, p_lb)
  } else {
    p_names <- p_base
  }
  
  
  # ================================================================
  # 4. Population truth for composite models under MI
  # ================================================================
  
  ## Sum composite true values:
  a_t_sum <- colSums(lam_mat)    # sum of loadings per time
  c_t_sum <- colSums(tau_mat)    # sum of intercepts per time
  E_St <- c_t_sum + a_t_sum * (mu_alpha + t_idx * mu_beta)
  
  ## Mean composite true values:
  a_t_mean <- colSums(lam_mat) / J  # mean of loadings per time
  c_t_mean <- colSums(tau_mat) / J  # mean of intercepts per time
  E_Mt <- c_t_mean + a_t_mean * (mu_alpha + t_idx * mu_beta)
  
  ## Linear projection truth 
  B_lin <- cbind(1, t_idx)
  sum_lin_true  <- as.numeric(solve(crossprod(B_lin), t(B_lin) %*% E_St))
  mean_lin_true <- as.numeric(solve(crossprod(B_lin), t(B_lin) %*% E_Mt))
  
  ## Quadratic true values
  B_quad <- cbind(1, t_idx, t_idx^2)
  quad_true_coefs <- as.numeric(solve(crossprod(B_quad), t(B_quad) %*% E_St))
  
  ## Latent basis true values 
  b_true     <- (E_St - E_St[1]) / (E_St[TT] - E_St[1])
  lb_k1_true <- E_St[1]
  lb_k2_true <- E_St[TT] - E_St[1]
  
  ## k1/k2 truth by model
  k1_true <- switch(model,
                    "latent"       = mu_alpha,
                    "mean"         = mean_lin_true[1],
                    "sum"          = sum_lin_true[1],
                    "mean_scaled"  = mu_alpha,
                    "sum_scaled"   = mu_alpha,
                    "quadratic"    = quad_true_coefs[1],
                    "latent_basis" = lb_k1_true)
  k2_true <- switch(model,
                    "latent"       = mu_beta,
                    "mean"         = mean_lin_true[2],
                    "sum"          = sum_lin_true[2],
                    "mean_scaled"  = mu_beta,
                    "sum_scaled"   = mu_beta,
                    "quadratic"    = quad_true_coefs[2],
                    "latent_basis" = lb_k2_true)
  
  ## phi truth by model
  phi_scale <- switch(model,
                      "sum" =,
                      "quadratic" =,
                      "latent_basis" = mean(a_t_sum)^2,
                      "mean" = mean(a_t_mean)^2,
                      1  # latent, scaled
  )
  
  phi1_true  <- phi_scale * var_alpha
  phi2_true  <- phi_scale * var_beta
  phi12_true <- phi_scale * cov_ab
  # ================================================================
  # 4. Simulation loop
  # ================================================================
  zval <- 1.96
  
  ci_in <- function(est, se, true) {
    as.numeric((est - zval * se) <= true & true <= (est + zval * se))
  }
  
  rbias_fun <- function(est, true) {
    if (isTRUE(all.equal(true, 0))) return(NA_real_)
    mean((est - true) / true, na.rm = TRUE)
  }
  
  results_list <- vector("list", length(N_vec))
  names(results_list) <- paste0("N", N_vec)
  
  
  for (n_idx in seq_along(N_vec)) {
    N <- N_vec[n_idx]
    cat(sprintf("[%s | N = %5d]  Simulation...\n", model, N))
    
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
      
      ## (3) Observed data
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
      if (model %in% c("mean")) {
        dat_fit <- data.frame(
          time1 = rowMeans(dat[,  1: 5]),
          time2 = rowMeans(dat[,  6:10]),
          time3 = rowMeans(dat[, 11:15]),
          time4 = rowMeans(dat[, 16:20]),
          time5 = rowMeans(dat[, 21:25])
        )
      } else if (model == "mean_scaled") {
        mean_t <- data.frame(
          time1 = rowMeans(dat[,  1: 5]),
          time2 = rowMeans(dat[,  6:10]),
          time3 = rowMeans(dat[, 11:15]),
          time4 = rowMeans(dat[, 16:20]),
          time5 = rowMeans(dat[, 21:25])
        )
        a_t <- colSums(lam_mat) / J 
        c_t <- colSums(tau_mat) / J 
        mean_t_lat <- sweep(mean_t, 2, c_t, "-")
        mean_t_lat <- sweep(mean_t_lat, 2, a_t, "/")
        dat_fit <- as.data.frame(mean_t_lat)
        names(dat_fit) <- paste0("time", 1:5)
        
      } else if (model == "sum") {
        dat_fit <- data.frame(
          time1 = rowSums(dat[,  1: 5]),
          time2 = rowSums(dat[,  6:10]),
          time3 = rowSums(dat[, 11:15]),
          time4 = rowSums(dat[, 16:20]),
          time5 = rowSums(dat[, 21:25])
        )
        
      } else if (model == "sum_scaled") {
        sum_t <- data.frame(
          time1 = rowSums(dat[,  1: 5]),
          time2 = rowSums(dat[,  6:10]),
          time3 = rowSums(dat[, 11:15]),
          time4 = rowSums(dat[, 16:20]),
          time5 = rowSums(dat[, 21:25])
        )
        a_t <- colSums(lam_mat)
        c_t <- colSums(tau_mat)
        sum_t_lat <- sweep(sum_t, 2, c_t, "-")
        sum_t_lat <- sweep(sum_t_lat, 2, a_t, "/")
        dat_fit <- as.data.frame(sum_t_lat)
        names(dat_fit) <- paste0("time", 1:5)
        
      } else if (model %in% c("quadratic", "latent_basis")) {
        dat_fit <- data.frame(
          time1 = rowSums(dat[,  1: 5]),
          time2 = rowSums(dat[,  6:10]),
          time3 = rowSums(dat[, 11:15]),
          time4 = rowSums(dat[, 16:20]),
          time5 = rowSums(dat[, 21:25])
        )
      } else {
        dat_fit <- dat  
      }
      
      
      # (C) Model fitting
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
        pe$op == "~~" & pe$label != "phi12" & 
          pe$label != "phi13" & pe$label != "phi23" &
          pe$est < 0                                     ~ 1,
        TRUE                                             ~ NA_real_
      )
      
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
      recall[i,  "phi11"]    <- gv("phi1",  "est")
      recall[i,  "phi11.SE"] <- gv("phi1",  "se")
      recall[i,  "phi11.p"]  <- gv("phi1",  "pvalue")
      recall[i,  "phi22"]    <- gv("phi2",  "est")
      recall[i,  "phi22.SE"] <- gv("phi2",  "se")
      recall[i,  "phi22.p"]  <- gv("phi2",  "pvalue")
      recall[i,  "phi12"]    <- gv("phi12", "est")
      recall[i,  "phi12.SE"] <- gv("phi12", "se")
      recall[i,  "phi12.p"]  <- gv("phi12", "pvalue")
      recall[i,  "cor12"]    <- gs("phi12", "est.std")
      recall[i,  "cor12.SE"] <- gs("phi12", "se")
      recall[i, "cor12.p"]   <- gs("phi12", "pvalue")
      recall[i,  "kappa1"]   <- gv("k1",    "est")
      recall[i,  "kappa1.SE"]<- gv("k1",    "se")
      recall[i,  "kappa1.p"] <- gv("k1",    "pvalue")
      recall[i,  "kappa2"]   <- gv("k2",    "est")
      recall[i,  "kappa2.SE"]<- gv("k2",    "se")
      recall[i,  "kappa2.p"] <- gv("k2",    "pvalue")
      
      ## Non-convergence flag
      recall[i, "err.mark"] <- ifelse(
        sum(pe$err.mark, na.rm = TRUE) >= 1L, 1L, NA_real_)
      
      ## 95% Coverage against generating values
      recall[i, "phi11.cov"] <- ci_in(recall[i, "phi11"], recall[i, "phi11.SE"], phi1_true)
      recall[i, "phi12.cov"] <- ci_in(recall[i, "phi12"], recall[i, "phi12.SE"], phi12_true)
      recall[i, "phi22.cov"] <- ci_in(recall[i, "phi22"], recall[i, "phi22.SE"], phi2_true)
      recall[i, "cor12.cov"]  <- ci_in(recall[i, "cor12"],  recall[i, "cor12.SE"], cor_ab)
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
      
      # (E) Model-specific parameters
      if (model == "quadratic") {
        recall[i, "kappa3"]    <- gv("k3",    "est")
        recall[i, "kappa3.SE"] <- gv("k3",    "se")
        recall[i, "kappa3.p"]  <- gv("k3",    "pvalue")
        recall[i, "phi33"]     <- gv("phi3",  "est")
        recall[i, "phi33.SE"]  <- gv("phi3",  "se")
        recall[i, "phi33.p"]   <- gv("phi3",  "pvalue")
        recall[i, "phi13"]     <- gv("phi13", "est")
        recall[i, "phi13.SE"]  <- gv("phi13", "se")
        recall[i, "phi13.p"]   <- gv("phi13", "pvalue")
        recall[i, "phi23"]     <- gv("phi23", "est")
        recall[i, "phi23.SE"]  <- gv("phi23", "se")
        recall[i, "phi23.p"]   <- gv("phi23", "pvalue")
        ## coverage: k3 against quadratic population truth
        recall[i, "kappa3.cov"] <- ci_in(
          recall[i, "kappa3"], recall[i, "kappa3.SE"], quad_true_coefs[3])
      }
      
      if (model == "latent_basis") {
        recall[i, "b2"]    <- gv("b2", "est")
        recall[i, "b2.SE"] <- gv("b2", "se")
        recall[i, "b3"]    <- gv("b3", "est")
        recall[i, "b3.SE"] <- gv("b3", "se")
        recall[i, "b4"]    <- gv("b4", "est")
        recall[i, "b4.SE"] <- gv("b4", "se")
      }
      
    } # end nrep
    
    # ================================================================
    # (F) Summary
    # ================================================================
    err.perc_ <- mean(!is.na(recall$err.mark))
    recall_ok <- dplyr::filter(recall, is.na(err.mark))
    
    ## --- Common summary (all models) ---
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
      cor12.bias = mean(recall_ok$cor12 - cor_ab,     na.rm = TRUE),  # cor은 scale-free라 그대로
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
      
      ## Coverage
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
      
      ## Scale-free quantities
      rho_IS.m = mean( recall_ok$phi12 / (sqrt(recall_ok$phi11) * sqrt(recall_ok$phi22)), na.rm = TRUE),
      d_slope.m = mean( recall_ok$kappa2 / sqrt(recall_ok$phi22), na.rm = TRUE),
      
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
    
    ## --- Quadratic-specific summary ---
    if (model == "quadratic") {
      summary_row <- cbind(summary_row, data.frame(
        k3.m       = mean(recall_ok$kappa3, na.rm = TRUE),
        k3.bias    = mean(recall_ok$kappa3 - quad_true_coefs[3], na.rm = TRUE),
        k3.rbias   = rbias_fun(recall_ok$kappa3, quad_true_coefs[3]),
        k3.rmse    = sqrt(mean((recall_ok$kappa3 - quad_true_coefs[3])^2, na.rm = TRUE)),
        k3.cov     = mean(recall_ok$kappa3.cov, na.rm = TRUE),
        k3.power   = mean(recall_ok$kappa3.p < 0.05, na.rm = TRUE),
        phi33.m    = mean(recall_ok$phi33, na.rm = TRUE),
        phi13.m    = mean(recall_ok$phi13, na.rm = TRUE),
        phi23.m    = mean(recall_ok$phi23, na.rm = TRUE),
        ## Population truth for mean structure
        quad.true.intercept = quad_true_coefs[1],
        quad.true.linear    = quad_true_coefs[2],
        quad.true.quad      = quad_true_coefs[3]
      ))
    }
    
    ## --- Latent basis-specific summary ---
    if (model == "latent_basis") {
      summary_row <- cbind(summary_row, data.frame(
        b2.m     = mean(recall_ok$b2, na.rm = TRUE),
        b3.m     = mean(recall_ok$b3, na.rm = TRUE),
        b4.m     = mean(recall_ok$b4, na.rm = TRUE),
        b2.se.m  = mean(recall_ok$b2.SE, na.rm = TRUE),
        b3.se.m  = mean(recall_ok$b3.SE, na.rm = TRUE),
        b4.se.m  = mean(recall_ok$b4.SE, na.rm = TRUE),
        ## Bias of basis estimates against true basis
        b2.bias  = mean(recall_ok$b2 - b_true[2], na.rm = TRUE),
        b3.bias  = mean(recall_ok$b3 - b_true[3], na.rm = TRUE),
        b4.bias  = mean(recall_ok$b4 - b_true[4], na.rm = TRUE),
        ## True basis values
        b2.true  = b_true[2],
        b3.true  = b_true[3],
        b4.true  = b_true[4],
        ## Linear basis for comparison
        b2.linear = 0.25,
        b3.linear = 0.50,
        b4.linear = 0.75,
        ## Total change (latent basis k2 = E[S_T] - E[S_1])
        lb.k2.true = lb_k2_true
      ))
    }
    
    results_list[[n_idx]] <- list(
      raw     = recall_ok,
      summary = summary_row
    )
    
    cat(sprintf("  -> Complete (Converged %d/%d, Non-converged %.1f%%)\n",
                nrow(recall_ok), nrep, err.perc_ * 100))
  } # end N_vec loop
  
  # ================================================================
  # 5. Aggregate
  # ================================================================
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
    ## Composite population truth (mean trajectory)
    composite_pop = list(
      E_St            = E_St,
      E_Mt            = E_Mt,
      a_t_sum         = a_t_sum,  # sum composite loading weights
      a_t_mean        = a_t_mean,   # mean composite loading weights
      c_t_sum         = c_t_sum,  # sum composite intercept sums
      c_t_mean        = c_t_mean,   # mean composite intercept means
      b_true          = b_true,
      quad_true_coefs = quad_true_coefs
    ),
    conditions = list(
      model = model,
      mni_location = mni_location,
      mni_size = mni_size,
      noninv_items = noninv_items,
      tau_mat = tau_mat,
      lam_mat = lam_mat,
      theta_mat = theta_mat
    )
  )
}
