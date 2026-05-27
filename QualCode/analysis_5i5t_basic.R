run5i5t <- function(nrep, model = c("latent", "composite"), 
                    k1, k2, phi1, phi2, phi12,
                    psi1 = NULL, psi2 = NULL, psi3 = NULL, psi4 = NULL, psi5 = NULL,
                    tau, sigma_mu = sqrt(1 - lambda^2), lambda, sd.lambda, nsize,
                    auto_var1 = TRUE, seed = 123, rho = 0, I = 5, tlen = 5,
                    inv.type = c("strict", "scalar", "metric"),
                    nu_drift_max   = 0.3,  # intercept drift: diff first-last
                    resid_mult_max = 1.5   # residual var drift: times first/last
) {
  model   <- match.arg(model)
  inv.type <- match.arg(inv.type)
  
  if (model %in% c("latent", "composite")) {
    if (I != 5L || tlen != 5L) {
      stop("This mode is for I = 5, tlen = 5. ")
    }
  }
  
  if (!is.null(seed)) set.seed(seed)
  
  ## Common parts: lam_vec, phi.mtx, nu_mat_pop, EC-population
  lam_vec <- truncnorm::rtruncnorm(I, a = lambda-.2, b = lambda+.2, 
                                   mean = lambda, sd = sd.lambda)
  
  phi.mtx <- matrix(c(phi1^2, phi12,
                      phi12,  phi2^2),
                    nrow = 2L, byrow = TRUE)
  
  ## (MI-0) non-invariant item 
  prop_noninv_items <- 0.4                # % non-invariance (I=5→2, I=10→4)
  n_noninv_items   <- max(1L, round(I * prop_noninv_items))
  noninv_items     <- tail(seq_len(I), n_noninv_items) # The last n items
  
  ## (MI-1) intercept structure: nu_mat (I x tlen)
  nu_base      <- rep(tau, I)       # base intercept
  
  if (inv.type %in% c("strict", "scalar")) {
    ## tau invariant
    nu_mat <- matrix(nu_base, nrow = I, ncol = tlen)
    
  } else if (inv.type == "metric") {
    ## tau drift among noninv_items 
    delta_t <- seq(0, nu_drift_max, length.out = tlen)
    nu_mat  <- matrix(nu_base, nrow = I, ncol = tlen)
    
    nu_mat_noninv <- nu_mat[noninv_items, , drop = FALSE]
    nu_mat_noninv <- nu_mat_noninv + matrix(delta_t, 
                                            nrow = length(noninv_items),
                                            ncol = tlen, byrow = TRUE)
    nu_mat[noninv_items, ] <- nu_mat_noninv
    
  } else {
    stop("Unknown inv.type: ", inv.type)
  }
  
  ## EC fitted population 
  ec_pop <- theta_from_theta_gen(
    phi.mtx   = phi.mtx,
    theta_gen = c(k1, k2),
    lam_vec   = lam_vec,
    nu_mat    = nu_mat,
    b         = NULL    # By default, linear basis
  )
  k1p   <- ec_pop$means[1L]; k2p <- ec_pop$means[2L]
  phi1p <- ec_pop$covs[1L,1L]; phi12p <- ec_pop$covs[2L,1L]; phi2p <- ec_pop$covs[2L,2L]
  
  zval <- 1.96  # 95% CI
  
  ##########################################################
  ## 1) 2nd-order LGCM (effect coding, strict model로 적합)
  ##########################################################
  if (model == "latent") {
    recall.1 <- data.frame(matrix(NA_real_, nrow = nrep, ncol = 21))
    colnames(recall.1) <- c(
      "phi11","phi11.SE","phi11.p",
      "phi22","phi22.SE","phi22.p",
      "phi12","phi12.SE","phi12.p",
      "kappa1","kappa1.SE","kappa1.p",
      "kappa2","kappa2.SE","kappa2.p",
      "err.mark",
      "phi11.cov","phi12.cov","phi22.cov",
      "kappa1.cov","kappa2.cov"
    )
    
    for (i in seq_len(nrep)) {
      if (!is.null(seed)) set.seed(seed + i)
      
      ## (1) Growth Factor
      phi_draw <- MASS::mvrnorm(nsize, mu = c(0,0), Sigma = phi.mtx)
      w1 <- phi_draw[,1]; w2 <- phi_draw[,2]
      xi1 <- k1 + w1;     xi2 <- k2 + w2
      
      
      ## (2) Latent factor residuals at each time point (eta_t의 오차)
      if (auto_var1) {
        At <- rbind(c(1,0), c(1,1), c(1,2), c(1,3), c(1,4))
        v_growth <- as.numeric(diag(At %*% phi.mtx %*% t(At)))
        buf <- 1e-12
        psi1 <- sqrt(pmax(1 - v_growth[1], buf))
        psi2 <- sqrt(pmax(1 - v_growth[2], buf))
        psi3 <- sqrt(pmax(1 - v_growth[3], buf))
        psi4 <- sqrt(pmax(1 - v_growth[4], buf))
        psi5 <- sqrt(pmax(1 - v_growth[5], buf))
      } # else: psi1~psi4 그대로 사용
      
      psi.mtx <- matrix(c(psi1^2, 0, 0, 0, 0,
                          0, psi2^2, 0, 0, 0,
                          0, 0, psi3^2, 0, 0,
                          0, 0, 0, psi4^2, 0,
                          0, 0, 0, 0, psi5^2),
                        nrow = 5L, byrow = TRUE)
      psi_disturb <- MASS::mvrnorm(nsize, mu = c(0,0,0,0,0), Sigma = psi.mtx)
      zeta1 <- psi_disturb[,1]; zeta2 <- psi_disturb[,2]
      zeta3 <- psi_disturb[,3]; zeta4 <- psi_disturb[,4]
      zeta5 <- psi_disturb[,5];
      
      eta1 <- xi1 + 0*xi2 + zeta1
      eta2 <- xi1 + 1*xi2 + zeta2
      eta3 <- xi1 + 2*xi2 + zeta3
      eta4 <- xi1 + 3*xi2 + zeta4
      eta5 <- xi1 + 4*xi2 + zeta5
      
      eta_mat <- cbind(eta1, eta2, eta3, eta4, eta5)  # nsize x tlen
      
      
      ## (MI-2) residual structure: sigma.mtx
      base_sig_diag <- rep(sigma_mu^2, I * tlen)
      
      if (inv.type %in% c("strict", "metric")) {
        ## residual var invariant
        sig_diag <- base_sig_diag
        
      } else if (inv.type == "scalar") {
        ## residual var drift among noninv_items
        scale_t <- seq(1.0, resid_mult_max, length.out = tlen)
        sig_diag <- base_sig_diag 
        
        for (tt in seq_len(tlen)) {
          idx_t <- (tt - 1L) * I + noninv_items   # location of non-invariant items
          sig_diag[idx_t] <- sigma_mu^2 * scale_t[tt]
        }
        
      } else {
        stop("Unknown inv.type: ", inv.type)
      }
      
      sigma.mtx <- diag(sig_diag, nrow = I * tlen)
      
      
      ## (3-1) Correlated uniqueness (rho) across adjacent time points
      {
        rvec <- if (length(rho) == 1L) rep(rho, I) else as.numeric(rho)
        stopifnot(length(rvec) == I)
        dth <- diag(sigma.mtx)
        for (ii in seq_len(I)) {
          for (tt in 1:(tlen - 1L)) {
            id1 <- ii + (tt - 1L) * I
            id2 <- ii +  tt      * I
            cov_it <- rvec[ii] * sqrt(dth[id1] * dth[id2])
            sigma.mtx[id1, id2] <- cov_it
            sigma.mtx[id2, id1] <- cov_it
          }
        }
        evmin <- min(eigen(sigma.mtx, symmetric = TRUE,
                           only.values = TRUE)$values)
        if (evmin <= 1e-10)
          stop("Non-PD measurement error covariance; reduce |rho| or drift.")
      }
      
      ## (3-2) Measurement error 
      errors <- MASS::mvrnorm(nsize, mu = rep(0, I * tlen), Sigma = sigma.mtx)
      e <- as.data.frame(errors)
      
      ## (3-3) y
      y <- list(
        # Time 1 (eta1)
        i1t1 = nu_mat[1,1] + lam_vec[1]*eta1 + e[,1],    
        i2t1 = nu_mat[2,1] + lam_vec[2]*eta1 + e[,2],
        i3t1 = nu_mat[3,1] + lam_vec[3]*eta1 + e[,3],
        i4t1 = nu_mat[4,1] + lam_vec[4]*eta1 + e[,4],
        i5t1 = nu_mat[5,1] + lam_vec[5]*eta1 + e[,5],
        
        # Time 2 (eta2)
        i1t2 = nu_mat[1,2] + lam_vec[1]*eta2 + e[,6],    
        i2t2 = nu_mat[2,2] + lam_vec[2]*eta2 + e[,7],
        i3t2 = nu_mat[3,2] + lam_vec[3]*eta2 + e[,8],
        i4t2 = nu_mat[4,2] + lam_vec[4]*eta2 + e[,9],
        i5t2 = nu_mat[5,2] + lam_vec[5]*eta2 + e[,10],
        
        # Time 3 (eta3)
        i1t3 = nu_mat[1,3] + lam_vec[1]*eta3 + e[,11],   
        i2t3 = nu_mat[2,3] + lam_vec[2]*eta3 + e[,12],
        i3t3 = nu_mat[3,3] + lam_vec[3]*eta3 + e[,13],
        i4t3 = nu_mat[4,3] + lam_vec[4]*eta3 + e[,14],
        i5t3 = nu_mat[5,3] + lam_vec[5]*eta3 + e[,15],
        
        # Time 4 (eta4)
        i1t4 = nu_mat[1,4] + lam_vec[1]*eta4 + e[,16],   
        i2t4 = nu_mat[2,4] + lam_vec[2]*eta4 + e[,17],
        i3t4 = nu_mat[3,4] + lam_vec[3]*eta4 + e[,18],
        i4t4 = nu_mat[4,4] + lam_vec[4]*eta4 + e[,19],
        i5t4 = nu_mat[5,4] + lam_vec[5]*eta4 + e[,20],
        
        # Time 5 (eta5)
        i1t5 = nu_mat[1,5] + lam_vec[1]*eta5 + e[,21],   
        i2t5 = nu_mat[2,5] + lam_vec[2]*eta5 + e[,22],
        i3t5 = nu_mat[3,5] + lam_vec[3]*eta5 + e[,23],
        i4t5 = nu_mat[4,5] + lam_vec[4]*eta5 + e[,24],
        i5t5 = nu_mat[5,5] + lam_vec[5]*eta5 + e[,25]
      )
      
      data <- as.data.frame(y)
      
      
      ## (4) lavaan  
      if (inv.type == "strict") {
        
        ## strict: intercept invariant + residual invariant 
        L.strc.EC <- '
  eta1 =~ (lmb1)*i1t1 + (lmb2)*i2t1 + (lmb3)*i3t1 + (lmb4)*i4t1 + (lmb5)*i5t1
  eta2 =~ (lmb1)*i1t2 + (lmb2)*i2t2 + (lmb3)*i3t2 + (lmb4)*i4t2 + (lmb5)*i5t2
  eta3 =~ (lmb1)*i1t3 + (lmb2)*i2t3 + (lmb3)*i3t3 + (lmb4)*i4t3 + (lmb5)*i5t3
  eta4 =~ (lmb1)*i1t4 + (lmb2)*i2t4 + (lmb3)*i3t4 + (lmb4)*i4t4 + (lmb5)*i5t4
  eta5 =~ (lmb1)*i1t5 + (lmb2)*i2t5 + (lmb3)*i3t5 + (lmb4)*i4t5 + (lmb5)*i5t5

  # intercepts:
  i1t1 ~ (nu1)*1; i2t1 ~ (nu2)*1; i3t1 ~ (nu3)*1; i4t1 ~ (nu4)*1; i5t1 ~ (nu5)*1
  i1t2 ~ (nu1)*1; i2t2 ~ (nu2)*1; i3t2 ~ (nu3)*1; i4t2 ~ (nu4)*1; i5t2 ~ (nu5)*1
  i1t3 ~ (nu1)*1; i2t3 ~ (nu2)*1; i3t3 ~ (nu3)*1; i4t3 ~ (nu4)*1; i5t3 ~ (nu5)*1
  i1t4 ~ (nu1)*1; i2t4 ~ (nu2)*1; i3t4 ~ (nu3)*1; i4t4 ~ (nu4)*1; i5t4 ~ (nu5)*1
  i1t5 ~ (nu1)*1; i2t5 ~ (nu2)*1; i3t5 ~ (nu3)*1; i4t5 ~ (nu4)*1; i5t5 ~ (nu5)*1

  # residual variances: 
  i1t1 ~~ (sig1)*i1t1; i2t1 ~~ (sig2)*i2t1; i3t1 ~~ (sig3)*i3t1; i4t1 ~~ (sig4)*i4t1; i5t1 ~~ (sig5)*i5t1
  i1t2 ~~ (sig1)*i1t2; i2t2 ~~ (sig2)*i2t2; i3t2 ~~ (sig3)*i3t2; i4t2 ~~ (sig4)*i4t2; i5t2 ~~ (sig5)*i5t2
  i1t3 ~~ (sig1)*i1t3; i2t3 ~~ (sig2)*i2t3; i3t3 ~~ (sig3)*i3t3; i4t3 ~~ (sig4)*i4t3; i5t3 ~~ (sig5)*i5t3
  i1t4 ~~ (sig1)*i1t4; i2t4 ~~ (sig2)*i2t4; i3t4 ~~ (sig3)*i3t4; i4t4 ~~ (sig4)*i4t4; i5t4 ~~ (sig5)*i5t4
  i1t5 ~~ (sig1)*i1t5; i2t5 ~~ (sig2)*i2t5; i3t5 ~~ (sig3)*i3t5; i4t5 ~~ (sig4)*i4t5; i5t5 ~~ (sig5)*i5t5

  # 2nd-order growth
  i =~ 1*eta1 + 1*eta2 + 1*eta3 + 1*eta4 + 1*eta5
  s =~ 0*eta1 + 1*eta2 + 2*eta3 + 3*eta4 + 4*eta5
  i ~~ (phi1)*i
  s ~~ (phi2)*s
  i ~~ (phi12)*s
  i ~ (k1)*1
  s ~ (k2)*1

  eta1 ~ 0*1; eta2 ~ 0*1; eta3 ~ 0*1; eta4 ~ 0*1;  eta5 ~ 0*1
  eta1 ~~ (psi1)*eta1; eta2 ~~ (psi2)*eta2; eta3 ~~ (psi3)*eta3; 
  eta4 ~~ (psi4)*eta4; eta5 ~~ (psi5)*eta5

  # effect coding: 
  lmb1 + lmb2 + lmb3 + lmb4 + lmb5 == 1
  nu1  + nu2  + nu3  + nu4  + nu5  == 0
  '
        
      } else if (inv.type == "scalar") {
        
        L.strc.EC <- '
  eta1 =~ (lmb1)*i1t1 + (lmb2)*i2t1 + (lmb3)*i3t1 + (lmb4)*i4t1 + (lmb5)*i5t1
  eta2 =~ (lmb1)*i1t2 + (lmb2)*i2t2 + (lmb3)*i3t2 + (lmb4)*i4t2 + (lmb5)*i5t2
  eta3 =~ (lmb1)*i1t3 + (lmb2)*i2t3 + (lmb3)*i3t3 + (lmb4)*i4t3 + (lmb5)*i5t3
  eta4 =~ (lmb1)*i1t4 + (lmb2)*i2t4 + (lmb3)*i3t4 + (lmb4)*i4t4 + (lmb5)*i5t4
  eta5 =~ (lmb1)*i1t5 + (lmb2)*i2t5 + (lmb3)*i3t5 + (lmb4)*i4t5 + (lmb5)*i5t5

  # intercepts: 
  i1t1 ~ (nu1)*1; i2t1 ~ (nu2)*1; i3t1 ~ (nu3)*1; i4t1 ~ (nu4)*1; i5t1 ~ (nu5)*1
  i1t2 ~ (nu1)*1; i2t2 ~ (nu2)*1; i3t2 ~ (nu3)*1; i4t2 ~ (nu4)*1; i5t2 ~ (nu5)*1
  i1t3 ~ (nu1)*1; i2t3 ~ (nu2)*1; i3t3 ~ (nu3)*1; i4t3 ~ (nu4)*1; i5t3 ~ (nu5)*1
  i1t4 ~ (nu1)*1; i2t4 ~ (nu2)*1; i3t4 ~ (nu3)*1; i4t4 ~ (nu4)*1; i5t4 ~ (nu5)*1
  i1t5 ~ (nu1)*1; i2t5 ~ (nu2)*1; i3t5 ~ (nu3)*1; i4t5 ~ (nu4)*1; i5t5 ~ (nu5)*1

  # residual variances: NON-invariant
  i1t1 ~~ (sig1)*i1t1; i1t2 ~~ (sig1)*i1t2; i1t3 ~~ (sig1)*i1t3; i1t4 ~~ (sig1)*i1t4; i1t5 ~~ (sig1)*i1t5
  i2t1 ~~ (sig2)*i2t1; i2t2 ~~ (sig2)*i2t2; i2t3 ~~ (sig2)*i2t3; i2t4 ~~ (sig2)*i2t4; i2t5 ~~ (sig2)*i2t5
  i3t1 ~~ (sig3)*i3t1; i3t2 ~~ (sig3)*i3t2; i3t3 ~~ (sig3)*i3t3; i3t4 ~~ (sig3)*i3t4; i3t5 ~~ (sig3)*i3t5

  i4t1 ~~ (sig41)*i4t1; i4t2 ~~ (sig42)*i4t2; i4t3 ~~ (sig43)*i4t3; i4t4 ~~ (sig44)*i4t4; i4t5 ~~ (sig45)*i4t5
  i5t1 ~~ (sig51)*i5t1; i5t2 ~~ (sig52)*i5t2; i5t3 ~~ (sig53)*i5t3; i5t4 ~~ (sig54)*i5t4; i5t5 ~~ (sig55)*i5t5

  # correlated uniqueness
  i1t1 ~~ i1t2;  i1t2 ~~ i1t3;  i1t3 ~~ i1t4;  i1t4 ~~ i1t5
  i2t1 ~~ i2t2;  i2t2 ~~ i2t3;  i2t3 ~~ i2t4;  i2t4 ~~ i2t5
  i3t1 ~~ i3t2;  i3t2 ~~ i3t3;  i3t3 ~~ i3t4;  i3t4 ~~ i3t5
  i4t1 ~~ i4t2;  i4t2 ~~ i4t3;  i4t3 ~~ i4t4;  i4t4 ~~ i4t5
  i5t1 ~~ i5t2;  i5t2 ~~ i5t3;  i5t3 ~~ i5t4;  i5t4 ~~ i5t5

  # 2nd-order growth
  i =~ 1*eta1 + 1*eta2 + 1*eta3 + 1*eta4 + 1*eta5
  s =~ 0*eta1 + 1*eta2 + 2*eta3 + 3*eta4 + 4*eta5
  i ~~ (phi1)*i
  s ~~ (phi2)*s
  i ~~ (phi12)*s
  i ~ (k1)*1
  s ~ (k2)*1

  eta1 ~ 0*1; eta2 ~ 0*1; eta3 ~ 0*1; eta4 ~ 0*1; eta5 ~ 0*1
  eta1 ~~ (psi1)*eta1; eta2 ~~ (psi2)*eta2; eta3 ~~ (psi3)*eta3; eta4 ~~ (psi4)*eta4; eta5 ~~ (psi5)*eta5

  # effect coding: 
  lmb1 + lmb2 + lmb3 + lmb4 + lmb5 == 1
  nu1  + nu2  + nu3  + nu4  + nu5  == 0
  '
        
      } else if (inv.type == "metric") {

        L.strc.EC <- '
  eta1 =~ (lmb1)*i1t1 + (lmb2)*i2t1 + (lmb3)*i3t1 + (lmb4)*i4t1 + (lmb5)*i5t1
  eta2 =~ (lmb1)*i1t2 + (lmb2)*i2t2 + (lmb3)*i3t2 + (lmb4)*i4t2 + (lmb5)*i5t2
  eta3 =~ (lmb1)*i1t3 + (lmb2)*i2t3 + (lmb3)*i3t3 + (lmb4)*i4t3 + (lmb5)*i5t3
  eta4 =~ (lmb1)*i1t4 + (lmb2)*i2t4 + (lmb3)*i3t4 + (lmb4)*i4t4 + (lmb5)*i5t4
  eta5 =~ (lmb1)*i1t5 + (lmb2)*i2t5 + (lmb3)*i3t5 + (lmb4)*i4t5 + (lmb5)*i5t5

  # intercepts: NON-invariant
  i1t1 ~ (nu1)*1;  i2t1 ~ (nu2)*1;  i3t1 ~ (nu3)*1;  i4t1 ~ (nu41)*1; i5t1 ~ (nu51)*1
  i1t2 ~ (nu1)*1;  i2t2 ~ (nu2)*1;  i3t2 ~ (nu3)*1;  i4t2 ~ (nu42)*1; i5t2 ~ (nu52)*1
  i1t3 ~ (nu1)*1;  i2t3 ~ (nu2)*1;  i3t3 ~ (nu3)*1;  i4t3 ~ (nu43)*1; i5t3 ~ (nu53)*1
  i1t4 ~ (nu1)*1;  i2t4 ~ (nu2)*1;  i3t4 ~ (nu3)*1;  i4t4 ~ (nu44)*1; i5t4 ~ (nu54)*1
  i1t5 ~ (nu1)*1;  i2t5 ~ (nu2)*1;  i3t5 ~ (nu3)*1;  i4t5 ~ (nu45)*1; i5t5 ~ (nu55)*1

  # residual variances: 
  i1t1 ~~ (sig1)*i1t1; i2t1 ~~ (sig2)*i2t1; i3t1 ~~ (sig3)*i3t1; i4t1 ~~ (sig4)*i4t1; i5t1 ~~ (sig5)*i5t1
  i1t2 ~~ (sig1)*i1t2; i2t2 ~~ (sig2)*i2t2; i3t2 ~~ (sig3)*i3t2; i4t2 ~~ (sig4)*i4t2; i5t2 ~~ (sig5)*i5t2
  i1t3 ~~ (sig1)*i1t3; i2t3 ~~ (sig2)*i2t3; i3t3 ~~ (sig3)*i3t3; i4t3 ~~ (sig4)*i4t3; i5t3 ~~ (sig5)*i5t3
  i1t4 ~~ (sig1)*i1t4; i2t4 ~~ (sig2)*i2t4; i3t4 ~~ (sig3)*i3t4; i4t4 ~~ (sig4)*i4t4; i5t4 ~~ (sig5)*i5t4
  i1t5 ~~ (sig1)*i1t5; i2t5 ~~ (sig2)*i2t5; i3t5 ~~ (sig3)*i3t5; i4t5 ~~ (sig4)*i4t5; i5t5 ~~ (sig5)*i5t5

  # correlated uniqueness 
  i1t1 ~~ i1t2;  i1t2 ~~ i1t3;  i1t3 ~~ i1t4;  i1t4 ~~ i1t5
  i2t1 ~~ i2t2;  i2t2 ~~ i2t3;  i2t3 ~~ i2t4;  i2t4 ~~ i2t5
  i3t1 ~~ i3t2;  i3t2 ~~ i3t3;  i3t3 ~~ i3t4;  i3t4 ~~ i3t5
  i4t1 ~~ i4t2;  i4t2 ~~ i4t3;  i4t3 ~~ i4t4;  i4t4 ~~ i4t5
  i5t1 ~~ i5t2;  i5t2 ~~ i5t3;  i5t3 ~~ i5t4;  i5t4 ~~ i5t5

  # 2nd-order growth
  i =~ 1*eta1 + 1*eta2 + 1*eta3 + 1*eta4 + 1*eta5
  s =~ 0*eta1 + 1*eta2 + 2*eta3 + 3*eta4 + 4*eta5
  i ~~ (phi1)*i
  s ~~ (phi2)*s
  i ~~ (phi12)*s
  i ~ (k1)*1
  s ~ (k2)*1

  eta1 ~ 0*1; eta2 ~ 0*1; eta3 ~ 0*1; eta4 ~ 0*1; eta5 ~ 0*1
  eta1 ~~ (psi1)*eta1; eta2 ~~ (psi2)*eta2; eta3 ~~ (psi3)*eta3; eta4 ~~ (psi4)*eta4; eta5 ~~ (psi5)*eta5

  # effect coding:
  lmb1 + lmb2 + lmb3 + lmb4 + lmb5 == 1

  nu1 + nu2 + nu3 + nu41 + nu51 == 0
  nu1 + nu2 + nu3 + nu42 + nu52 == 0
  nu1 + nu2 + nu3 + nu43 + nu53 == 0
  nu1 + nu2 + nu3 + nu44 + nu54 == 0
  nu1 + nu2 + nu3 + nu45 + nu55 == 0
  '
        
      } else {
        stop("Unknown inv.type: ", inv.type)
      }
      
      
      fit <- lavaan::sem(
        L.strc.EC, data = data,
        auto.fix.first  = FALSE,
        auto.fix.single = FALSE,
        effect.coding   = TRUE,
        control = list(reltol = 1e-12, xtol_rel = 1e-12)
      )
      out_summary <- summary(fit, standardized=TRUE, rsq=TRUE)
      pe  <- as.data.frame(lavaan::parameterEstimates(fit, standardized = TRUE))
      pe  <- pe[pe$label != "", c("op","label","est","se","pvalue","std.all")]
      
      ## non-convergence
      pe$err.mark <- dplyr::case_when(
        Reduce(`|`, lapply(pe, is.na))            ~ 1,
        pe$op == "~~" & pe$std.all >=  1.000001   ~ 1,
        pe$op == "~~" & pe$std.all <= -1.000001   ~ 1,
        pe$op == "~~" & pe$label != "phi12" & pe$est < 0 ~ 1,
        TRUE                                      ~ NA_real_
      )
      
      gv <- function(lbl, col) pe[pe$label == lbl, col, drop = TRUE]
      recall.1[i,1]  <- gv("phi1","est")
      recall.1[i,2]  <- gv("phi1","se")
      recall.1[i,3]  <- gv("phi1","pvalue")
      recall.1[i,4]  <- gv("phi2","est")
      recall.1[i,5]  <- gv("phi2","se")
      recall.1[i,6]  <- gv("phi2","pvalue")
      recall.1[i,7]  <- gv("phi12","est")
      recall.1[i,8]  <- gv("phi12","se")
      recall.1[i,9]  <- gv("phi12","pvalue")
      recall.1[i,10] <- gv("k1","est")
      recall.1[i,11] <- gv("k1","se")
      recall.1[i,12] <- gv("k1","pvalue")
      recall.1[i,13] <- gv("k2","est")
      recall.1[i,14] <- gv("k2","se")
      recall.1[i,15] <- gv("k2","pvalue")
      recall.1[i,16] <- ifelse(sum(pe$err.mark, na.rm = TRUE) >= 1, 1, NA)
      
      ## (5) 95% coverage
      ci_in <- function(est, se, true) {
        lo <- est - zval * se; hi <- est + zval * se
        as.numeric(true >= lo & true <= hi)
      }
      recall.1[i,17] <- ci_in(recall.1[i,1],  recall.1[i,2],  phi1p)  # phi11.cov
      recall.1[i,18] <- ci_in(recall.1[i,7],  recall.1[i,8],  phi12p) # phi12.cov
      recall.1[i,19] <- ci_in(recall.1[i,4],  recall.1[i,5],  phi2p)  # phi22.cov
      recall.1[i,20] <- ci_in(recall.1[i,10], recall.1[i,11], k1p)    # kappa1.cov
      recall.1[i,21] <- ci_in(recall.1[i,13], recall.1[i,14], k2p)    # kappa2.cov
    } # for i
    
    err.perc_ <- mean(!is.na(recall.1$err.mark))
    na.perc_  <- mean( is.na(recall.1$err.mark))
    
    recall.1 <- dplyr::filter(recall.1, is.na(err.mark))
    recall.1 <- dplyr::mutate(
      recall.1,
      phi11.mark   = (phi11.p <= .05),
      phi12.mark   = (phi12.p <= .05),
      phi22.mark   = (phi22.p <= .05),
      kappa1.mark  = (kappa1.p <= .05),
      kappa2.mark  = (kappa2.p <= .05),
      
      phi11.m      = mean(phi11, na.rm = TRUE),
      phi12.m      = mean(phi12, na.rm = TRUE),
      phi22.m      = mean(phi22, na.rm = TRUE),
      kappa1.m     = mean(kappa1, na.rm = TRUE),
      kappa2.m     = mean(kappa2, na.rm = TRUE),
      
      phi11.SE.m   = mean(phi11.SE, na.rm = TRUE),
      phi12.SE.m   = mean(phi12.SE, na.rm = TRUE),
      phi22.SE.m   = mean(phi22.SE, na.rm = TRUE),
      kappa1.SE.m  = mean(kappa1.SE, na.rm = TRUE),
      kappa2.SE.m  = mean(kappa2.SE, na.rm = TRUE),
      
      phi11.power  = mean(phi11.mark, na.rm = TRUE),
      phi12.power  = mean(phi12.mark, na.rm = TRUE),
      phi22.power  = mean(phi22.mark, na.rm = TRUE),
      kappa1.power = mean(kappa1.mark, na.rm = TRUE),
      kappa2.power = mean(kappa2.mark, na.rm = TRUE),
      
      phi11.cov.rate = mean(phi11.cov, na.rm = TRUE),
      phi12.cov.rate = mean(phi12.cov, na.rm = TRUE),
      phi22.cov.rate = mean(phi22.cov, na.rm = TRUE),
      kappa1.cov.rate= mean(kappa1.cov, na.rm = TRUE),
      kappa2.cov.rate= mean(kappa2.cov, na.rm = TRUE),
      
      corr         = phi12 / (sqrt(phi11) * sqrt(phi22)),
      
      phi11.bias   = mean(phi11 - phi1p),
      phi12.bias   = mean(phi12 - phi12p),
      phi22.bias   = mean(phi22 - phi2p),
      kappa1.bias  = mean(kappa1 - k1p),
      kappa2.bias  = mean(kappa2 - k2p),
      
      phi11.R.bias = mean((phi11 - phi1p)/phi1p),
      phi12.R.bias = mean((phi12 - phi12p)/phi12p),
      phi22.R.bias = mean((phi22 - phi2p)/phi2p),
      kappa1.R.bias= mean((kappa1 - k1p)/k1p),
      kappa2.R.bias= mean((kappa2 - k2p)/k2p),
      
      phi11.rmse   = sqrt(mean((phi11 - phi1p)^2)),
      phi12.rmse   = sqrt(mean((phi12 - phi12p)^2)),
      phi22.rmse   = sqrt(mean((phi22 - phi2p)^2)),
      kappa1.rmse  = sqrt(mean((kappa1 - k1p)^2)),
      kappa2.rmse  = sqrt(mean((kappa2 - k2p)^2)),
      
      phi11.rel.rmse = sqrt(mean((phi11 - phi1p)^2)) / phi1p,
      phi12.rel.rmse = sqrt(mean((phi12 - phi12p)^2)) / phi12p,
      phi22.rel.rmse = sqrt(mean((phi22 - phi2p)^2)) / phi2p,
      kappa1.rel.rmse = sqrt(mean((kappa1 - k1p)^2)) / k1p,
      kappa2.rel.rmse = sqrt(mean((kappa2 - k2p)^2)) / k2p,
      
      phi11.coverage = mean((phi11 - 1.96 * phi11.SE <= phi1p) & (phi11 + 1.96 * phi11.SE >= phi1p), na.rm = TRUE),
      phi12.coverage = mean((phi12 - 1.96 * phi12.SE <= phi12p) & (phi12 + 1.96 * phi12.SE >= phi12p), na.rm = TRUE),
      phi22.coverage = mean((phi22 - 1.96 * phi22.SE <= phi2p) & (phi22 + 1.96 * phi22.SE >= phi2p), na.rm = TRUE),
      kappa1.coverage = mean((kappa1 - 1.96 * kappa1.SE <= k1p) & (kappa1 + 1.96 * kappa1.SE >= k1p), na.rm = TRUE),
      kappa2.coverage = mean((kappa2 - 1.96 * kappa2.SE <= k2p) & (kappa2 + 1.96 * kappa2.SE >= k2p), na.rm = TRUE),
      
      phi11.ARB = mean(abs(phi11 - phi1p))/phi1p,
      phi12.ARB = mean(abs(phi12 - phi12p))/phi12p,
      phi22.ARB = mean(abs(phi22 - phi2p))/phi2p,
      kappa1.ARB = mean(abs(kappa1 - k1p))/k1p,
      kappa2.ARB = mean(abs(kappa2 - k2p))/k2p,
      
      err.perc = err.perc_,
      admissible = na.perc_
    )
    
    return(list(
      data    = data,
      output  = out_summary,
      summary = recall.1,
      ec_pop  = ec_pop
    ))
  } # end latent
  
  ##########################################################
  ## 2) Composite (sum score 1-LGCM)
  ##########################################################
  if (model == "composite") {
    recall.c <- data.frame(matrix(NA_real_, nrow = nrep, ncol = 21))
    colnames(recall.c) <- c(
      "phi11","phi11.SE","phi11.p",
      "phi22","phi22.SE","phi22.p",
      "phi12","phi12.SE","phi12.p",
      "alpha1","alpha1.SE","alpha1.p",
      "alpha2","alpha2.SE","alpha2.p",
      "err.mark",
      "phi11.cov","phi12.cov","phi22.cov",
      "alpha1.cov","alpha2.cov"
    )
    
    for (i in seq_len(nrep)) {
      if (!is.null(seed)) set.seed(seed + i)
      
      ## (1) Growth Factor
      phi_draw <- MASS::mvrnorm(nsize, mu = c(0,0), Sigma = phi.mtx)
      w1 <- phi_draw[,1]; w2 <- phi_draw[,2]
      xi1 <- k1 + w1;     xi2 <- k2 + w2
      
      
      ## (2) Latent factor residuals at each time point (eta_t의 오차)
      if (auto_var1) {
        At <- rbind(c(1,0), c(1,1), c(1,2), c(1,3), c(1,4))
        v_growth <- as.numeric(diag(At %*% phi.mtx %*% t(At)))
        buf <- 1e-12
        psi1 <- sqrt(pmax(1 - v_growth[1], buf))
        psi2 <- sqrt(pmax(1 - v_growth[2], buf))
        psi3 <- sqrt(pmax(1 - v_growth[3], buf))
        psi4 <- sqrt(pmax(1 - v_growth[4], buf))
        psi5 <- sqrt(pmax(1 - v_growth[5], buf))
      } # else: psi1~psi4 그대로 사용
      
      psi.mtx <- matrix(c(psi1^2, 0, 0, 0, 0,
                          0, psi2^2, 0, 0, 0,
                          0, 0, psi3^2, 0, 0,
                          0, 0, 0, psi4^2, 0,
                          0, 0, 0, 0, psi5^2),
                        nrow = 5L, byrow = TRUE)
      psi_disturb <- MASS::mvrnorm(nsize, mu = c(0,0,0,0,0), Sigma = psi.mtx)
      zeta1 <- psi_disturb[,1]; zeta2 <- psi_disturb[,2]
      zeta3 <- psi_disturb[,3]; zeta4 <- psi_disturb[,4]
      zeta5 <- psi_disturb[,5];
      
      eta1 <- xi1 + 0*xi2 + zeta1
      eta2 <- xi1 + 1*xi2 + zeta2
      eta3 <- xi1 + 2*xi2 + zeta3
      eta4 <- xi1 + 3*xi2 + zeta4
      eta5 <- xi1 + 4*xi2 + zeta5
      
      eta_mat <- cbind(eta1, eta2, eta3, eta4, eta5)  # nsize x tlen
      
      
      ## (MI-2) residual structure: sigma.mtx
      base_sig_diag <- rep(sigma_mu^2, I * tlen)
      
      if (inv.type %in% c("strict", "metric")) {
        ## residual var invariant
        sig_diag <- base_sig_diag
        
      } else if (inv.type == "scalar") {
        ## residual var drift among noninv_items
        scale_t <- seq(1.0, resid_mult_max, length.out = tlen)
        sig_diag <- base_sig_diag 
        
        for (tt in seq_len(tlen)) {
          idx_t <- (tt - 1L) * I + noninv_items   # location of non-invariant items
          sig_diag[idx_t] <- sigma_mu^2 * scale_t[tt]
        }
        
      } else {
        stop("Unknown inv.type: ", inv.type)
      }
      
      sigma.mtx <- diag(sig_diag, nrow = I * tlen)
      
      
      ## (3-1) Correlated uniqueness (rho) across adjacent time points
      {
        rvec <- if (length(rho) == 1L) rep(rho, I) else as.numeric(rho)
        stopifnot(length(rvec) == I)
        dth <- diag(sigma.mtx)
        for (ii in seq_len(I)) {
          for (tt in 1:(tlen - 1L)) {
            id1 <- ii + (tt - 1L) * I
            id2 <- ii +  tt      * I
            cov_it <- rvec[ii] * sqrt(dth[id1] * dth[id2])
            sigma.mtx[id1, id2] <- cov_it
            sigma.mtx[id2, id1] <- cov_it
          }
        }
        evmin <- min(eigen(sigma.mtx, symmetric = TRUE,
                           only.values = TRUE)$values)
        if (evmin <= 1e-10)
          stop("Non-PD measurement error covariance; reduce |rho| or drift.")
      }
      
      ## (3-2) Measurement error 
      errors <- MASS::mvrnorm(nsize, mu = rep(0, I * tlen), Sigma = sigma.mtx)
      e <- as.data.frame(errors)
      
      ## (3-3) y
      y <- list(
        # Time 1 (eta1)
        i1t1 = nu_mat[1,1] + lam_vec[1]*eta1 + e[,1],    
        i2t1 = nu_mat[2,1] + lam_vec[2]*eta1 + e[,2],
        i3t1 = nu_mat[3,1] + lam_vec[3]*eta1 + e[,3],
        i4t1 = nu_mat[4,1] + lam_vec[4]*eta1 + e[,4],
        i5t1 = nu_mat[5,1] + lam_vec[5]*eta1 + e[,5],
        
        # Time 2 (eta2)
        i1t2 = nu_mat[1,2] + lam_vec[1]*eta2 + e[,6],    
        i2t2 = nu_mat[2,2] + lam_vec[2]*eta2 + e[,7],
        i3t2 = nu_mat[3,2] + lam_vec[3]*eta2 + e[,8],
        i4t2 = nu_mat[4,2] + lam_vec[4]*eta2 + e[,9],
        i5t2 = nu_mat[5,2] + lam_vec[5]*eta2 + e[,10],
        
        # Time 3 (eta3)
        i1t3 = nu_mat[1,3] + lam_vec[1]*eta3 + e[,11],   
        i2t3 = nu_mat[2,3] + lam_vec[2]*eta3 + e[,12],
        i3t3 = nu_mat[3,3] + lam_vec[3]*eta3 + e[,13],
        i4t3 = nu_mat[4,3] + lam_vec[4]*eta3 + e[,14],
        i5t3 = nu_mat[5,3] + lam_vec[5]*eta3 + e[,15],
        
        # Time 4 (eta4)
        i1t4 = nu_mat[1,4] + lam_vec[1]*eta4 + e[,16],   
        i2t4 = nu_mat[2,4] + lam_vec[2]*eta4 + e[,17],
        i3t4 = nu_mat[3,4] + lam_vec[3]*eta4 + e[,18],
        i4t4 = nu_mat[4,4] + lam_vec[4]*eta4 + e[,19],
        i5t4 = nu_mat[5,4] + lam_vec[5]*eta4 + e[,20],
        
        # Time 5 (eta5)
        i1t5 = nu_mat[1,5] + lam_vec[1]*eta5 + e[,21],   
        i2t5 = nu_mat[2,5] + lam_vec[2]*eta5 + e[,22],
        i3t5 = nu_mat[3,5] + lam_vec[3]*eta5 + e[,23],
        i4t5 = nu_mat[4,5] + lam_vec[4]*eta5 + e[,24],
        i5t5 = nu_mat[5,5] + lam_vec[5]*eta5 + e[,25]
      )
      
      data <- as.data.frame(y)
      
      time1 <- rowSums(data[,  1: 5])
      time2 <- rowSums(data[,  6:10])
      time3 <- rowSums(data[, 11:15])
      time4 <- rowSums(data[, 16:20])
      time5 <- rowSums(data[, 21:25])
      data.comp <- data.frame(time1,time2,time3,time4,time5)
      
      
      ## (4) lavaan composite growth 
      C.strc.EC <- '
        i =~ 1*time1 + 1*time2 + 1*time3 + 1*time4 + 1*time5
        s =~ 0*time1 + 1*time2 + 2*time3 + 3*time4 + 4*time5
        i ~~ (phi1)*i
        s ~~ (phi2)*s
        i ~~ (phi12)*s
        i ~  (a1)*1
        s ~  (a2)*1
      '
      fit <- lavaan::growth(C.strc.EC, data = data.comp)
      out_summary <- summary(fit, standardized=TRUE, rsq=TRUE)
      pe  <- as.data.frame(lavaan::parameterEstimates(fit, standardized = TRUE))
      pe  <- pe[pe$label != "", c("op","label","est","se","pvalue","std.all")]
      pe$err.mark <- dplyr::case_when(
        Reduce(`|`, lapply(pe, is.na)) ~ 1,
        pe$op=="~~" & pe$std.all >=  1.000001 ~ 1,
        pe$op=="~~" & pe$std.all <= -1.000001 ~ 1,
        pe$op=="~~" & pe$label!="phi12" & pe$est < 0 ~ 1,
        TRUE ~ NA_real_
      )
      gv <- function(lbl, col) pe[pe$label==lbl, col, drop = TRUE]
      
      recall.c[i,1]  <- gv("phi1","est")
      recall.c[i,2]  <- gv("phi1","se")
      recall.c[i,3]  <- gv("phi1","pvalue")
      recall.c[i,4]  <- gv("phi2","est")
      recall.c[i,5]  <- gv("phi2","se")
      recall.c[i,6]  <- gv("phi2","pvalue")
      recall.c[i,7]  <- gv("phi12","est")
      recall.c[i,8]  <- gv("phi12","se")
      recall.c[i,9]  <- gv("phi12","pvalue")
      recall.c[i,10] <- gv("a1","est")
      recall.c[i,11] <- gv("a1","se")
      recall.c[i,12] <- gv("a1","pvalue")
      recall.c[i,13] <- gv("a2","est")
      recall.c[i,14] <- gv("a2","se")
      recall.c[i,15] <- gv("a2","pvalue")
      recall.c[i,16] <- ifelse(sum(pe$err.mark, na.rm=TRUE) >= 1, 1, NA)
      
      ## 95% coverage
      ci_in <- function(est, se, true) {
        lo <- est - zval*se; hi <- est + zval*se
        as.numeric(true >= lo & true <= hi)
      }
      recall.c[i,17] <- ci_in(recall.c[i,1],  recall.c[i,2],  phi1p)  # phi11.cov
      recall.c[i,18] <- ci_in(recall.c[i,7],  recall.c[i,8],  phi12p) # phi12.cov
      recall.c[i,19] <- ci_in(recall.c[i,4],  recall.c[i,5],  phi2p)  # phi22.cov
      recall.c[i,20] <- ci_in(recall.c[i,10], recall.c[i,11], k1p)    # alpha1.cov
      recall.c[i,21] <- ci_in(recall.c[i,13], recall.c[i,14], k2p)    # alpha2.cov
    } # for i
    
    err.perc_ <- mean(!is.na(recall.c$err.mark))
    na.perc_  <- mean( is.na(recall.c$err.mark))
    
    recall.c <- dplyr::filter(recall.c, is.na(err.mark))
    recall.c <- dplyr::mutate(
      recall.c,
      phi11.mark   = (phi11.p <= .05),
      phi12.mark   = (phi12.p <= .05),
      phi22.mark   = (phi22.p <= .05),
      alpha1.mark  = (alpha1.p <= .05),
      alpha2.mark  = (alpha2.p <= .05),
      
      phi11.m      = mean(phi11, na.rm = TRUE),
      phi12.m      = mean(phi12, na.rm = TRUE),
      phi22.m      = mean(phi22, na.rm = TRUE),
      alpha1.m     = mean(alpha1, na.rm = TRUE),
      alpha2.m     = mean(alpha2, na.rm = TRUE),
      
      phi11.SE.m   = mean(phi11.SE, na.rm = TRUE),
      phi12.SE.m   = mean(phi12.SE, na.rm = TRUE),
      phi22.SE.m   = mean(phi22.SE, na.rm = TRUE),
      alpha1.SE.m  = mean(alpha1.SE, na.rm = TRUE),
      alpha2.SE.m  = mean(alpha2.SE, na.rm = TRUE),
      
      phi11.power  = mean(phi11.mark, na.rm = TRUE),
      phi12.power  = mean(phi12.mark, na.rm = TRUE),
      phi22.power  = mean(phi22.mark, na.rm = TRUE),
      kappa1.power = mean(alpha1.mark, na.rm = TRUE),
      kappa2.power = mean(alpha2.mark, na.rm = TRUE),
      
      phi11.cov.rate = mean(phi11.cov, na.rm = TRUE),
      phi12.cov.rate = mean(phi12.cov, na.rm = TRUE),
      phi22.cov.rate = mean(phi22.cov, na.rm = TRUE),
      kappa1.cov.rate= mean(alpha1.cov, na.rm = TRUE),
      kappa2.cov.rate= mean(alpha2.cov, na.rm = TRUE),
      
      corr         = phi12 / (sqrt(phi11) * sqrt(phi22)),
      
      phi11.bias   = mean(phi11 - phi1p),
      phi12.bias   = mean(phi12 - phi12p),
      phi22.bias   = mean(phi22 - phi2p),
      kappa1.bias  = mean(alpha1 - k1p),
      kappa2.bias  = mean(alpha2 - k2p),
      
      phi11.R.bias = mean((phi11 - phi1p)/phi1p),
      phi12.R.bias = mean((phi12 - phi12p)/phi12p),
      phi22.R.bias = mean((phi22 - phi2p)/phi2p),
      kappa1.R.bias= mean((alpha1 - k1p)/k1p),
      kappa2.R.bias= mean((alpha2 - k2p)/k2p),
      
      phi11.rmse   = sqrt(mean((phi11 - phi1p)^2)),
      phi12.rmse   = sqrt(mean((phi12 - phi12p)^2)),
      phi22.rmse   = sqrt(mean((phi22 - phi2p)^2)),
      kappa1.rmse  = sqrt(mean((alpha1 - k1p)^2)),
      kappa2.rmse  = sqrt(mean((alpha2 - k2p)^2)),
      
      phi11.rel.rmse = sqrt(mean((phi11 - phi1p)^2)) / phi1p,
      phi12.rel.rmse = sqrt(mean((phi12 - phi12p)^2)) / phi12p,
      phi22.rel.rmse = sqrt(mean((phi22 - phi2p)^2)) / phi2p,
      kappa1.rel.rmse = sqrt(mean((alpha1 - k1p)^2)) / k1p,
      kappa2.rel.rmse = sqrt(mean((alpha2 - k2p)^2)) / k2p,
      
      phi11.coverage = mean((phi11 - 1.96 * phi11.SE <= phi1p) & (phi11 + 1.96 * phi11.SE >= phi1p), na.rm = TRUE),
      phi12.coverage = mean((phi12 - 1.96 * phi12.SE <= phi12p) & (phi12 + 1.96 * phi12.SE >= phi12p), na.rm = TRUE),
      phi22.coverage = mean((phi22 - 1.96 * phi22.SE <= phi2p) & (phi22 + 1.96 * phi22.SE >= phi2p), na.rm = TRUE),
      kappa1.coverage = mean((alpha1 - 1.96 * alpha1.SE <= k1p) & (alpha1 + 1.96 * alpha1.SE >= k1p), na.rm = TRUE),
      kappa2.coverage = mean((alpha2 - 1.96 * alpha2.SE <= k2p) & (alpha2 + 1.96 * alpha2.SE >= k2p), na.rm = TRUE),
      
      phi11.ARB = mean(abs(phi11 - phi1p))/phi1p,
      phi12.ARB = mean(abs(phi12 - phi12p))/phi12p,
      phi22.ARB = mean(abs(phi22 - phi2p))/phi2p,
      kappa1.ARB = mean(abs(alpha1 - k1p))/k1p,
      kappa2.ARB = mean(abs(alpha2 - k2p))/k2p,
      
      err.perc = err.perc_,
      admissible  = na.perc_
    )
    
    return(list(
      data    = data.comp,
      output  = out_summary,
      summary = recall.c,
      ec_pop  = ec_pop
    ))
  } # end composite.
  
  stop("model should be either 'latnet' or 'composite'")
}


