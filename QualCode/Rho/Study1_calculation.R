### Growth Factor Covariances =================

## 1) Growth basis (B) 
make_B <- function(T, b = NULL) {
  if (is.null(b)) b <- 0:(T-1)  # default slope basis = t-1
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
  I <- length(lam_vec)                                 # number of items
  if (!is.matrix(tau_mat)) tau_mat <- as.matrix(tau_mat)  # intercepts vector
  if (nrow(tau_mat) != I) stop("nrow(tau_mat) must equal length(lam_vec)")
  
  tlen <- ncol(tau_mat)    # time points
  
  ## (1) growth basis B
  B <- make_B(tlen, b = b)  # linear change 0:(tlen-1)
  
  ## (2) loading matrix (lam_vec_t)
  Lambda_mat <- matrix(lam_vec, nrow = I, ncol = tlen)
  
  ## (3) Sum of loadings s_t
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

