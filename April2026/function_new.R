# Lavaan output convergence check !
.lav_admissible <- function(fit) {
  if (is.null(fit) || !inherits(fit, "lavaan")) return(FALSE)
  if (!isTRUE(lavaan::lavInspect(fit, "converged"))) return(FALSE)
  if (!isTRUE(lavaan::lavTech(fit, "post.check"))) return(FALSE) 
  pe <- tryCatch(
    as.data.frame(lavaan::parameterEstimates(fit, standardized = TRUE)),
    error = function(e) NULL)
  if (is.null(pe) || nrow(pe) == 0L) return(FALSE) # no output
  if (any(is.na(pe$est))) return(FALSE)            # no estimates
  is_var <- pe$op == "~~" & pe$lhs == pe$rhs       # about variance...
  if (any(is_var & pe$est < 0, na.rm = TRUE)) return(FALSE) # negative variance
  is_cov_like <- pe$op == "~~"         # about correlation >1 or <-1
  if (any(is_cov_like & abs(pe$std.all) > 1.000001, na.rm = TRUE)) return(FALSE)
  TRUE
}


run <- function(
    nrep = 2000,
    N_vec = c(200, 500, 800),
    model = c("latent", "mean", "sum", "mean_scaled", "sum_scaled"),
    
    ## Design ----
    TT = 5,    # time points
    J = 5,    # number of indicators / time
    
    ## Growth Parameters ----
    mu_alpha  = 1,                                                                                      
    mu_beta   = 0.125,
    var_alpha = 1,                                                                                      
    var_beta  = 0.2025,                                                                                 
    cov_ab    = 0.18,                                        
    cor_ab    = cov_ab / sqrt(var_alpha * var_beta),                                                 
    
    ## Measurement Parameters                       
    target_loading = c(.7, .15, .3),     # Manipulate measurement quality: .7, .5, .3       
    tau_un    = 1,                                                                    
    lambda_un = 1,      
    sd.lambda = 0,                                                                                   
    rho = c(0, 0.15),                                                                                         
    psi_mat = c(1.097^2, 0.98^2, 0.95^2, 0.93^2, 0.91^2),  # latent disturbances
    
    ## Measurement Non-Invariance
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
  
  ##### 1. MNI setup #####
  mni_params <- switch(
    mni_size,
    none = list(n_noninv = 0, d_tau = 0, d_lam = 0, d_res = 0),
    weak = list(n_noninv = 1, d_tau = 0.25, d_lam = 0.15, d_res = 0.30),
    strong = list(n_noninv = 3, d_tau = 0.50, d_lam = 0.30, d_res = 0.60)
  )
  
  ## the last n Non-Invariant items
  noninv_items <- if (mni_params$n_noninv > 0L) {
    tail(seq_len(J), mni_params$n_noninv)
  } else {
    integer(0)
  }
  
  ## drift properties: non-linear change
  t_idx  <- 0:(TT - 1)
  u  <- t_idx / (TT - 1)    # intervals in 0~1 
  drift_t <- u^2
  
  ## Parameters
  tau_vec <- rep(tau_un, times=J)
  tau_mat <- matrix(tau_un, nrow = J, ncol = TT)
  
  set.seed(seed) 
  lam_vec <- if (sd.lambda == 0) {
    rep(lambda_un, J)
  } else {
    truncnorm::rtruncnorm(
      J,
      a    = lambda_un - .2,
      b    = lambda_un + .2,
      mean = lambda_un,
      sd   = sd.lambda)
  }
  lam_mat <- matrix(lam_vec, nrow = J, ncol = TT)
  
  var_eta   <- var_alpha + psi_mat[1]    # targetted at t=1
  theta_vec <- (lam_vec)^2 * var_eta * (1 / target_loading^2 - 1)
  theta_mat <- matrix(rep(theta_vec, TT), nrow = J)
  
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
    D %*% R %*% D      # cov = sd * cor * sd
  }
  
  ##### 2. Lavaan model syntax ######
  cu_items <- '
# Correlated uniqueness
i1t1 ~~ i1t2; i1t2 ~~ i1t3; i1t3 ~~ i1t4; i1t4 ~~ i1t5
i2t1 ~~ i2t2; i2t2 ~~ i2t3; i2t3 ~~ i2t4; i2t4 ~~ i2t5
i3t1 ~~ i3t2; i3t2 ~~ i3t3; i3t3 ~~ i3t4; i3t4 ~~ i3t5
i4t1 ~~ i4t2; i4t2 ~~ i4t3; i4t3 ~~ i4t4; i4t4 ~~ i4t5
i5t1 ~~ i5t2; i5t2 ~~ i5t3; i5t3 ~~ i5t4; i5t4 ~~ i5t5
'
  cu_sums <- '
# Correlated residuals 
time1 ~~ time2; time2 ~~ time3; time3 ~~ time4; time4 ~~ time5
'
  
  ### 2.1. Second-order LGCM (latent) ##### 
  if (model == "latent") {
    
    ## Structural part (2nd-order LGCM)
    lgcm_struct <- '
# Second-order LGCM
I =~ 1*f1 + 1*f2 + 1*f3 + 1*f4 + 1*f5
S =~ 0*f1 + 1*f2 + 2*f3 + 3*f4 + 4*f5

I ~~ (phi1)*I
S ~~ (phi2)*S
I ~~ (phi12)*S
I ~  (k1)*1
S ~  (k2)*1

# Factor intercepts fixed 0
f1 ~ 0*1; f2 ~ 0*1; f3 ~ 0*1; f4 ~ 0*1; f5 ~ 0*1

# Latent residual variances
f1 ~~ s1*f1; f2 ~~ s2*f2; f3 ~~ s3*f3; f4 ~~ s4*f4; f5 ~~ s5*f5
'
    
    ## Measurement parts
    # mni_size == "none" (fullly strict invariance)
    meas_inv_full <- '
# Measurement model (fixed)
f1 =~ 1*i1t1 + l2*i2t1 + l3*i3t1 + l4*i4t1 + l5*i5t1
f2 =~ 1*i1t2 + l2*i2t2 + l3*i3t2 + l4*i4t2 + l5*i5t2
f3 =~ 1*i1t3 + l2*i2t3 + l3*i3t3 + l4*i4t3 + l5*i5t3
f4 =~ 1*i1t4 + l2*i2t4 + l3*i3t4 + l4*i4t4 + l5*i5t4
f5 =~ 1*i1t5 + l2*i2t5 + l3*i3t5 + l4*i4t5 + l5*i5t5

# Intercepts (fixed)
i1t1 ~ 1*1; i1t2 ~ 1*1; i1t3 ~ 1*1; i1t4 ~ 1*1; i1t5 ~ 1*1
i2t1 ~ int2*1; i2t2 ~ int2*1; i2t3 ~ int2*1; i2t4 ~ int2*1; i2t5 ~ int2*1
i3t1 ~ int3*1; i3t2 ~ int3*1; i3t3 ~ int3*1; i3t4 ~ int3*1; i3t5 ~ int3*1
i4t1 ~ int4*1; i4t2 ~ int4*1; i4t3 ~ int4*1; i4t4 ~ int4*1; i4t5 ~ int4*1
i5t1 ~ int5*1; i5t2 ~ int5*1; i5t3 ~ int5*1; i5t4 ~ int5*1; i5t5 ~ int5*1
'
    
    # Loading MNI 
    meas_load_weak <- '
# Measurement model (weak / free)
f1 =~ 1*i1t1 + l2*i2t1 + l3*i3t1 + l4*i4t1 + l51*i5t1
f2 =~ 1*i1t2 + l2*i2t2 + l3*i3t2 + l4*i4t2 + l52*i5t2
f3 =~ 1*i1t3 + l2*i2t3 + l3*i3t3 + l4*i4t3 + l53*i5t3
f4 =~ 1*i1t4 + l2*i2t4 + l3*i3t4 + l4*i4t4 + l54*i5t4
f5 =~ 1*i1t5 + l2*i2t5 + l3*i3t5 + l4*i4t5 + l55*i5t5

# Intercepts (fixed)
i1t1 ~ 1*1; i1t2 ~ 1*1; i1t3 ~ 1*1; i1t4 ~ 1*1; i1t5 ~ 1*1
i2t1 ~ int2*1; i2t2 ~ int2*1; i2t3 ~ int2*1; i2t4 ~ int2*1; i2t5 ~ int2*1
i3t1 ~ int3*1; i3t2 ~ int3*1; i3t3 ~ int3*1; i3t4 ~ int3*1; i3t5 ~ int3*1
i4t1 ~ int4*1; i4t2 ~ int4*1; i4t3 ~ int4*1; i4t4 ~ int4*1; i4t5 ~ int4*1
i5t1 ~ int5*1; i5t2 ~ int5*1; i5t3 ~ int5*1; i5t4 ~ int5*1; i5t5 ~ int5*1
'
    meas_load_strong <- '
# Measurement model (strong / free)
f1 =~ 1*i1t1 + l2*i2t1 + l31*i3t1 + l41*i4t1 + l51*i5t1
f2 =~ 1*i1t2 + l2*i2t2 + l32*i3t2 + l42*i4t2 + l52*i5t2
f3 =~ 1*i1t3 + l2*i2t3 + l33*i3t3 + l43*i4t3 + l53*i5t3
f4 =~ 1*i1t4 + l2*i2t4 + l34*i3t4 + l44*i4t4 + l54*i5t4
f5 =~ 1*i1t5 + l2*i2t5 + l35*i3t5 + l45*i4t5 + l55*i5t5

# Intercepts (fixed)
i1t1 ~ 1*1; i1t2 ~ 1*1; i1t3 ~ 1*1; i1t4 ~ 1*1; i1t5 ~ 1*1
i2t1 ~ int2*1; i2t2 ~ int2*1; i2t3 ~ int2*1; i2t4 ~ int2*1; i2t5 ~ int2*1
i3t1 ~ int3*1; i3t2 ~ int3*1; i3t3 ~ int3*1; i3t4 ~ int3*1; i3t5 ~ int3*1
i4t1 ~ int4*1; i4t2 ~ int4*1; i4t3 ~ int4*1; i4t4 ~ int4*1; i4t5 ~ int4*1
i5t1 ~ int5*1; i5t2 ~ int5*1; i5t3 ~ int5*1; i5t4 ~ int5*1; i5t5 ~ int5*1
'
    
    # Intercept MNI 
    meas_int_weak <- '
# Measurement model (fixed)
f1 =~ 1*i1t1 + l2*i2t1 + l3*i3t1 + l4*i4t1 + l5*i5t1
f2 =~ 1*i1t2 + l2*i2t2 + l3*i3t2 + l4*i4t2 + l5*i5t2
f3 =~ 1*i1t3 + l2*i2t3 + l3*i3t3 + l4*i4t3 + l5*i5t3
f4 =~ 1*i1t4 + l2*i2t4 + l3*i3t4 + l4*i4t4 + l5*i5t4
f5 =~ 1*i1t5 + l2*i2t5 + l3*i3t5 + l4*i4t5 + l5*i5t5

# Intercepts (weak / free)
i1t1 ~ 1*1; i1t2 ~ 1*1; i1t3 ~ 1*1; i1t4 ~ 1*1; i1t5 ~ 1*1
i2t1 ~ int2*1; i2t2 ~ int2*1; i2t3 ~ int2*1; i2t4 ~ int2*1; i2t5 ~ int2*1
i3t1 ~ int3*1; i3t2 ~ int3*1; i3t3 ~ int3*1; i3t4 ~ int3*1; i3t5 ~ int3*1
i4t1 ~ int4*1; i4t2 ~ int4*1; i4t3 ~ int4*1; i4t4 ~ int4*1; i4t5 ~ int4*1
i5t1 ~ int51*1; i5t2 ~ int52*1; i5t3 ~ int53*1; i5t4 ~ int54*1; i5t5 ~ int55*1
'
    meas_int_strong <- '
# Measurement model (fixed)
f1 =~ 1*i1t1 + l2*i2t1 + l3*i3t1 + l4*i4t1 + l5*i5t1
f2 =~ 1*i1t2 + l2*i2t2 + l3*i3t2 + l4*i4t2 + l5*i5t2
f3 =~ 1*i1t3 + l2*i2t3 + l3*i3t3 + l4*i4t3 + l5*i5t3
f4 =~ 1*i1t4 + l2*i2t4 + l3*i3t4 + l4*i4t4 + l5*i5t4
f5 =~ 1*i1t5 + l2*i2t5 + l3*i3t5 + l4*i4t5 + l5*i5t5

# Intercepts (strong / free)
i1t1 ~ 1*1; i1t2 ~ 1*1; i1t3 ~ 1*1; i1t4 ~ 1*1; i1t5 ~ 1*1
i2t1 ~ int2*1; i2t2 ~ int2*1; i2t3 ~ int2*1; i2t4 ~ int2*1; i2t5 ~ int2*1
i3t1 ~ int31*1; i3t2 ~ int32*1; i3t3 ~ int33*1; i3t4 ~ int34*1; i3t5 ~ int35*1
i4t1 ~ int41*1; i4t2 ~ int42*1; i4t3 ~ int43*1; i4t4 ~ int44*1; i4t5 ~ int45*1
i5t1 ~ int51*1; i5t2 ~ int52*1; i5t3 ~ int53*1; i5t4 ~ int54*1; i5t5 ~ int55*1
'
    
    ## Residual variances
    resvar_inv_full <- '
# Residual variances (fixed)
i1t1 ~~ e1*i1t1; i1t2 ~~ e1*i1t2; i1t3 ~~ e1*i1t3; i1t4 ~~ e1*i1t4; i1t5 ~~ e1*i1t5
i2t1 ~~ e2*i2t1; i2t2 ~~ e2*i2t2; i2t3 ~~ e2*i2t3; i2t4 ~~ e2*i2t4; i2t5 ~~ e2*i2t5
i3t1 ~~ e3*i3t1; i3t2 ~~ e3*i3t2; i3t3 ~~ e3*i3t3; i3t4 ~~ e3*i3t4; i3t5 ~~ e3*i3t5
i4t1 ~~ e4*i4t1; i4t2 ~~ e4*i4t2; i4t3 ~~ e4*i4t3; i4t4 ~~ e4*i4t4; i4t5 ~~ e4*i4t5
i5t1 ~~ e5*i5t1; i5t2 ~~ e5*i5t2; i5t3 ~~ e5*i5t3; i5t4 ~~ e5*i5t4; i5t5 ~~ e5*i5t5
'
    resvar_res_weak <- '
# Residual variances (weak / free)
i1t1 ~~ e1*i1t1; i1t2 ~~ e1*i1t2; i1t3 ~~ e1*i1t3; i1t4 ~~ e1*i1t4; i1t5 ~~ e1*i1t5
i2t1 ~~ e2*i2t1; i2t2 ~~ e2*i2t2; i2t3 ~~ e2*i2t3; i2t4 ~~ e2*i2t4; i2t5 ~~ e2*i2t5
i3t1 ~~ e3*i3t1; i3t2 ~~ e3*i3t2; i3t3 ~~ e3*i3t3; i3t4 ~~ e3*i3t4; i3t5 ~~ e3*i3t5
i4t1 ~~ e4*i4t1; i4t2 ~~ e4*i4t2; i4t3 ~~ e4*i4t3; i4t4 ~~ e4*i4t4; i4t5 ~~ e4*i4t5
i5t1 ~~ e51*i5t1; i5t2 ~~ e52*i5t2; i5t3 ~~ e53*i5t3; i5t4 ~~ e54*i5t4; i5t5 ~~ e55*i5t5
'
    resvar_res_strong <- '
# Residual variances (strong / freee)
i1t1 ~~ e1*i1t1; i1t2 ~~ e1*i1t2; i1t3 ~~ e1*i1t3; i1t4 ~~ e1*i1t4; i1t5 ~~ e1*i1t5
i2t1 ~~ e2*i2t1; i2t2 ~~ e2*i2t2; i2t3 ~~ e2*i2t3; i2t4 ~~ e2*i2t4; i2t5 ~~ e2*i2t5
i3t1 ~~ e31*i3t1; i3t2 ~~ e32*i3t2; i3t3 ~~ e33*i3t3; i3t4 ~~ e34*i3t4; i3t5 ~~ e35*i3t5
i4t1 ~~ e41*i4t1; i4t2 ~~ e42*i4t2; i4t3 ~~ e43*i4t3; i4t4 ~~ e44*i4t4; i4t5 ~~ e45*i4t5
i5t1 ~~ e51*i5t1; i5t2 ~~ e52*i5t2; i5t3 ~~ e53*i5t3; i5t4 ~~ e54*i5t4; i5t5 ~~ e55*i5t5
'
    
    ## Pick blocks per (location × size)
    if (mni_size == "none" || mni_location == "none") {
      meas_part  <- meas_inv_full
      resvar_inv <- resvar_inv_full
      
    } else if (mni_location == "loading" && mni_size == "weak")   {
      meas_part  <- meas_load_weak
      resvar_inv <- resvar_inv_full
      
    } else if (mni_location == "loading" && mni_size == "strong") {
      meas_part  <- meas_load_strong
      resvar_inv <- resvar_inv_full
      
    } else if (mni_location == "intercept" && mni_size == "weak")   {
      meas_part  <- meas_int_weak
      resvar_inv <- resvar_inv_full
      
    } else if (mni_location == "intercept" && mni_size == "strong") {
      meas_part  <- meas_int_strong
      resvar_inv <- resvar_inv_full
      
    } else if (mni_location == "residual"  && mni_size == "weak")   {
      meas_part  <- meas_inv_full
      resvar_inv <- resvar_res_weak
      
    } else if (mni_location == "residual"  && mni_size == "strong") {
      meas_part  <- meas_inv_full
      resvar_inv <- resvar_res_strong
    }
    
    # Assemble
    lav_model <- paste0(meas_part, resvar_inv, lgcm_struct)
    if (rho != 0) lav_model <- paste0(lav_model, cu_items)
    
    
    ### 2.2. Sum-score / Mean-score LGCM #####
  } else {
    lav_model <- '
i =~ 1*time1 + 1*time2 + 1*time3 + 1*time4 + 1*time5
s =~ 0*time1 + 1*time2 + 2*time3 + 3*time4 + 4*time5

i ~~ (phi1)*i
s ~~ (phi2)*s
i ~~ (phi12)*s
i ~  (k1)*1
s ~  (k2)*1

time1 ~~ s1*time1; time2 ~~ s2*time2; time3 ~~ s3*time3; 
time4 ~~ s4*time4; time5 ~~ s5*time5
'
    if (rho != 0) lav_model <- paste0(lav_model, cu_sums)
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
  
   p_names <- p_base

  
  
  ##### 4. Population truth for composite models under MI #####
  # Sum 
  a_t_sum <- colSums(lam_mat)    # sum of loadings 
  c_t_sum <- colSums(tau_mat)    # sum of intercepts
  E_St <- c_t_sum + a_t_sum * (mu_alpha + t_idx * mu_beta)
  
  # Mean 
  a_t_mean <- colSums(lam_mat) / J 
  c_t_mean <- colSums(tau_mat) / J 
  E_Mt <- c_t_mean + a_t_mean * (mu_alpha + t_idx * mu_beta)
  
  ## Truth by linear projection
  B_lin <- cbind(1, t_idx)
  sum_lin_true  <- as.numeric(solve(crossprod(B_lin), t(B_lin) %*% E_St))
  mean_lin_true <- as.numeric(solve(crossprod(B_lin), t(B_lin) %*% E_Mt))
  
  ## k1/k2 truth by model
  k1_true <- switch(model,
                    "latent"       = mu_alpha,
                    "mean"         = mean_lin_true[1],
                    "sum"          = sum_lin_true[1],
                    "mean_scaled"  = mu_alpha,
                    "sum_scaled"   = mu_alpha)
  k2_true <- switch(model,
                    "latent"       = mu_beta,
                    "mean"         = mean_lin_true[2],
                    "sum"          = sum_lin_true[2],
                    "mean_scaled"  = mu_beta,
                    "sum_scaled"   = mu_beta)
  
  ## phi truth by model
  phi_scale <- switch(model, "sum" = mean(a_t_sum)^2,
                      "mean" = mean(a_t_mean)^2, 
                      1  # latent, scaled
                      )
  
  phi1_true  <- phi_scale * var_alpha
  phi2_true  <- phi_scale * var_beta
  phi12_true <- phi_scale * cov_ab
  
  
  ##### 5. Simulation loop #####
  zval <- 1.96
  SS_t <- sum((t_idx - mean(t_idx))^2)   # for slope reliability
  resid_labs <- paste0("s", seq_len(TT))
  
  ci_in <- function(est, se, true) {  # confidence interval
    as.numeric((est - zval * se) <= true & true <= (est + zval * se))
  }
  
  rbias_fun <- function(est, true) {  # relative bias
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
      ## (A.1) Growth Factor
      phi.mtx <- matrix(c(var_alpha, cov_ab,
                          cov_ab,   var_beta), 2, 2, byrow = TRUE)
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
      
      if (!.lav_admissible(fit)) {
        recall[i, "err.mark"] <- 1L
        next  # skip the left over for loop codes. 
      }
      recall[i, "err.mark"] <- 0L
      
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
      
    } # end nrep
    
    ######  6. Summary #####
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
      
      ## Scale-free quantities
      rho_IS.m = mean( recall_ok$phi12 / (sqrt(recall_ok$phi11) * sqrt(recall_ok$phi22)), na.rm = TRUE),
      d_slope.m = mean( recall_ok$kappa2 / sqrt(recall_ok$phi22), na.rm = TRUE),
      std_tot_eff.m = mean( recall_ok$kappa2 * (TT-1) / sqrt(recall_ok$phi11), na.rm = TRUE),
      slope_rel.m = mean( recall_ok$phi22 / (recall_ok$phi22  + recall_ok$sigma2_bar / SS_t) , na.rm = TRUE),                   
    
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
    
    cat(sprintf("  -> Complete (Converged %d/%d, Non-converged %.1f%%)\n",
                nrow(recall_ok), nrep, err.perc_ * 100))
  } # end N_vec loop
  
  
  ##### 5. Aggregate #####
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
    lav_out = lavaan::summary(fit, standardized=TRUE),
    ## Composite population truth (mean trajectory)
    composite_pop = list(
      E_St            = E_St,
      E_Mt            = E_Mt,
      a_t_sum         = a_t_sum,  # sum composite loading weights
      a_t_mean        = a_t_mean, # mean composite loading weights
      c_t_sum         = c_t_sum,  # sum composite intercept sums
      c_t_mean        = c_t_mean  # mean composite intercept means
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
