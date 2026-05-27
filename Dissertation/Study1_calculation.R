### Growth Factor Covariances =================

## 1) Growth basis (B) 
make_B <- function(T, b = NULL) {
  if (is.null(b)) b <- 0:(T-1)    # default slope basis = t-1
  stopifnot(length(b) == T)
  cbind(Intercept = 1, Slope = b)
}

## 2) Sum of factor loadings 
##    Lambda_mat:  (JJ x T) mat 
sum_loadings_by_time <- function(Lambda_mat = NULL) {
  colSums(as.matrix(Lambda_mat))
}


## 3) Make M 
make_M <- function(B, s) {
  # Diagonal matrix of sum of factor loadings
  D <- diag(as.numeric(s), nrow = length(s))
  solve(crossprod(B), t(B) %*% D %*% B)  # (B'B) * X = B'D B
}                                        # X = (B'B)^(-1) B'D B


## 4) Covariance by EC:  Σ_xi^(ec) = M Σ_xi^(gen) M'
transform_Sigma_xi <- function(Phi_gen, B, s) {
  D <- diag(as.numeric(s), nrow = length(s))
  M <- make_M(B, s)
  
  Phi_ec <- M %*% Phi_gen %*% t(M)
  
  Phi_ec
}


## 5) Re-parameterization to EC fitted pop
theta_from_theta_gen <- function(phi.mtx,
                                 theta_gen,
                                 lam_vec,
                                 tau_mat,
                                 b = NULL) {
  I <- nrow(tau_mat)     # number of items
  tlen <- ncol(tau_mat)    # time points
  
  ## (1) growth basis: B
  B <- make_B(tlen, b = b)  # linear change 0:(tlen-1)
  
  ## (2) loading matrix 
  Lambda_mat <- matrix(lam_vec, nrow = I, ncol = tlen)
  
  ## (3) Sum of loadings: s_t
  s <- sum_loadings_by_time(Lambda_mat)
  
  ## (4) make M (the complex one)
  M <- make_M(B, s)
  
  ## (5) Sum of intercepts c_t 
  cvec <- colSums(tau_mat) 
  
  ## (6) phi_ec = M * phi_gen * M'
  ec_cov <- transform_Sigma_xi(phi.mtx, B, s)
  
  ## (7) θ_ec = M θ_gen + (B'B)^(-1) B' c
  part_M <- as.numeric(M %*% matrix(theta_gen, ncol = 1L))
  part_c <- as.numeric(solve(crossprod(B), t(B) %*% cvec))
  ec_means <- part_M + part_c
  
  list(means = ec_means, covs = ec_cov)
}


### Population moments under the generating model =================
## Used for analytic plim (pseudo-true values) of each fitted model.
## Mirrors the data-generating mechanism in study1_function_new.R (lines 229-250).

## Residual covariance with correlated uniqueness (AR Lag-1), item j over time.
## (Identical to the local helper inside run(); lifted to top level for reuse.)
make_resid_cov <- function(var_e, TT, rho) {
  R <- diag(TT)
  for (t in seq_len(TT - 1L)) {
    R[t, t + 1L] <- rho     # AR Lag-1 only
    R[t + 1L, t] <- rho
  }
  D <- diag(sqrt(var_e), nrow = TT)
  D %*% R %*% D             # cov = sd * cor * sd
}

## Exact population covariance/mean of the JT observed items, and of the
## T sum scores. Ordering matches col_names: variable index = (t-1)*J + j.
##   Cov(eta) = B Phi_gen B' + diag(psi)
##   Sigma_item = Lambda Cov(eta) Lambda' + Theta
##   E[y_jt]    = tau_jt + lambda_jt E[eta_t],   E[eta] = B mu
##   sum score: aggregate with A (A[t,(t-1)J+j] = 1)
pop_moments <- function(phi.mtx, mu_vec, psi_mat, lam_mat, tau_mat, theta_mat,
                        rho, TT, J, b = NULL) {
  B       <- make_B(TT, b)                          # TT x 2 growth basis
  Cov_eta <- B %*% phi.mtx %*% t(B) + diag(psi_mat[seq_len(TT)], TT)
  E_eta   <- as.numeric(B %*% matrix(mu_vec, ncol = 1L))   # length TT

  ## Loading matrix Lambda: (J*TT) x TT, time-major ordering
  Lambda <- matrix(0, J * TT, TT)
  for (t in seq_len(TT))
    for (j in seq_len(J))
      Lambda[(t - 1L) * J + j, t] <- lam_mat[j, t]

  ## Residual covariance Theta (within-item AR lag-1; cross-item = 0)
  Theta <- matrix(0, J * TT, J * TT)
  for (j in seq_len(J)) {
    Sig_ej <- make_resid_cov(theta_mat[j, ], TT, rho)
    for (t in seq_len(TT))
      for (u in seq_len(TT))
        Theta[(t - 1L) * J + j, (u - 1L) * J + j] <- Sig_ej[t, u]
  }

  ## Item-level population moments
  Sigma_item <- Lambda %*% Cov_eta %*% t(Lambda) + Theta
  mu_item    <- numeric(J * TT)
  for (t in seq_len(TT))
    for (j in seq_len(J))
      mu_item[(t - 1L) * J + j] <- tau_mat[j, t] + lam_mat[j, t] * E_eta[t]

  ## Sum-score aggregation: A (TT x J*TT)
  A <- matrix(0, TT, J * TT)
  for (t in seq_len(TT)) A[t, ((t - 1L) * J + 1L):(t * J)] <- 1

  list(Sigma_item = Sigma_item,
       mu_item    = mu_item,
       Sigma_sum  = A %*% Sigma_item %*% t(A),
       mu_sum     = as.numeric(A %*% mu_item),
       Cov_eta    = Cov_eta)
}

