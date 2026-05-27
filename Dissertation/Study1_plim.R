## =====================================================================
## Study 1 - Analytic plim (pseudo-true values) of the fitted models
## ---------------------------------------------------------------------
## Fits the SAME lavaan models used in study1_function_new.R directly to
## the exact population moments of the generating model. The resulting
## point estimates are the probability limits (plim) of the ML estimators,
## i.e. the values the finite-N estimates converge to as N -> Inf.
##
## Purpose (reviewer issue #2: large-N consistency):
##   * correctly specified cells  -> plim(EC) = phi_ec  (consistency proof)
##   * TT=3, rho=0                 -> plim = phi_ec for BOTH  => the large
##                                    finite-N bias there is purely
##                                    finite-sample (vanishes as N grows)
##   * rho=.15 with CU omitted     -> plim(sum) - phi_ec != 0  => genuine
##                                    asymptotic (projection) bias
##
## Point estimates are invariant to sample.nobs; nobs only scales SE/chisq.
## sample.cov.rescale = FALSE keeps the fed moments exact.
## =====================================================================

suppressMessages({ library(lavaan) })
source("Study1_calculation.R")   # make_B, theta_from_theta_gen, pop_moments, make_resid_cov

## ---- fixed growth / measurement parameters (study1_function_new.R defaults) ----
MU_ALPHA  <- 1
MU_BETA   <- 0.125
VAR_ALPHA <- 1
VAR_BETA  <- 0.2025
COV_AB    <- 0.18
COR_AB    <- COV_AB / sqrt(VAR_ALPHA * VAR_BETA)
TAU_UN    <- 0.1
LAMBDA_UN <- 0.2
PSI_MAT   <- c(1.097^2, 0.98^2, 0.95^2, 0.93^2, 0.91^2, 0.89^2, 0.87^2, 0.85^2)

## ---- build cell-specific parameter set (mirrors run() lines 61-91) ----
build_params <- function(TT, J, sd.lambda, target_loading) {
  phi.mtx <- matrix(c(VAR_ALPHA, COV_AB,
                      COV_AB,    VAR_BETA), 2, 2, byrow = TRUE)
  tau_mat <- matrix(TAU_UN, nrow = J, ncol = TT)

  lam_vec <- if (sd.lambda == 0) {
    rep(LAMBDA_UN, J)
  } else {
    intervals <- seq(-(J - 1) / 2, (J - 1) / 2, by = 1)
    ratio     <- sd.lambda / sd(intervals)
    LAMBDA_UN + intervals * ratio
  }
  lam_mat <- matrix(lam_vec, nrow = J, ncol = TT)

  var_eta   <- VAR_ALPHA + PSI_MAT[1]
  theta_vec <- (lam_vec)^2 * var_eta * (1 / target_loading^2 - 1)
  theta_mat <- matrix(rep(theta_vec, TT), nrow = J)

  list(phi.mtx = phi.mtx, tau_mat = tau_mat, lam_vec = lam_vec,
       lam_mat = lam_mat, theta_mat = theta_mat)
}

## ---- build the exact lavaan syntax used in run() (lines 93-173) ----
build_syntax <- function(model, TT, J, rho) {
  load_t <- function(t, J) paste0("f", t, " =~ ", paste(sprintf("lmb%d*i%dt%d", 1:J, 1:J, t), collapse = " + "))
  int_t  <- function(t, J) paste(sprintf("i%dt%d ~ nu%d*1", 1:J, t, 1:J), collapse = "; ")
  res_j  <- function(j, TT) paste(sprintf("i%dt%d ~~ e%d*i%dt%d", j, 1:TT, j, j, 1:TT), collapse = "; ")
  cu_j   <- function(j, TT) paste(sprintf("i%dt%d ~~ i%dt%d", j, 1:(TT-1), j, 2:TT), collapse = "; ")

  if (model == "EC_latent") {
    I_eq  <- paste0("I =~ ", paste(sprintf("1*f%d", 1:TT), collapse = " + "))
    S_eq  <- paste0("S =~ ", paste(sprintf("%d*f%d", 0:(TT-1), 1:TT), collapse = " + "))
    f_int <- paste(sprintf("f%d ~ 0*1", 1:TT), collapse = "; ")
    f_var <- paste(sprintf("f%d ~~ s%d*f%d", 1:TT, 1:TT, 1:TT), collapse = "; ")
    meas  <- paste(sapply(1:TT, load_t, J = J),   collapse = "\n")
    ints  <- paste(sapply(1:TT, int_t,  J = J),   collapse = "\n")
    rvars <- paste(sapply(1:J,  res_j,  TT = TT), collapse = "\n")
    ec_con <- paste0(paste(sprintf("lmb%d", 1:J), collapse = " + "), " == 1\n",
                     paste(sprintf("nu%d",  1:J), collapse = " + "), " == 0")
    lav_model <- paste("# Second-order LGCM (effect coding)",
                       I_eq, S_eq,
                       "I ~~ (phi1)*I", "S ~~ (phi2)*S", "I ~~ (phi12)*S",
                       "I ~ (k1)*1", "S ~ (k2)*1",
                       f_int, f_var,
                       "# Measurement", meas, "# Intercepts", ints,
                       "# Residual variances", rvars, "# EC constraints", ec_con,
                       sep = "\n")
    if (rho != 0 && TT != 3) {
      cu_items <- paste(sapply(1:J, cu_j, TT = TT), collapse = "\n")
      lav_model <- paste(lav_model, cu_items, sep = "\n")
    }
  } else {
    i_eq  <- paste0("i =~ ", paste(sprintf("1*time%d", 1:TT), collapse = " + "))
    s_eq  <- paste0("s =~ ", paste(sprintf("%d*time%d", 0:(TT-1), 1:TT), collapse = " + "))
    t_var <- paste(sprintf("time%d ~~ s%d*time%d", 1:TT, 1:TT, 1:TT), collapse = "; ")
    lav_model <- paste(i_eq, s_eq,
                       "i ~~ (phi1)*i", "s ~~ (phi2)*s", "i ~~ (phi12)*s",
                       "i ~ (k1)*1", "s ~ (k2)*1", t_var, sep = "\n")
    if (rho != 0 && TT != 3) {
      cu_sums <- paste(sprintf("time%d ~~ time%d", 1:(TT-1), 2:TT), collapse = "; ")
      lav_model <- paste(lav_model, cu_sums, sep = "\n")
    }
  }
  lav_model
}

## ---- fit a model to population moments; return plim point estimates ----
fit_plim <- function(model, lav_model, mom, TT, J, NOBS = 1e5) {
  if (model == "EC_latent") {
    col_names <- character(J * TT)
    for (t in seq_len(TT)) for (j in seq_len(J))
      col_names[(t - 1L) * J + j] <- paste0("i", j, "t", t)
    Sigma <- mom$Sigma_item; dimnames(Sigma) <- list(col_names, col_names)
    mu    <- mom$mu_item;    names(mu) <- col_names
    fit <- lavaan::sem(lav_model, sample.cov = Sigma, sample.mean = mu,
                       sample.nobs = NOBS, sample.cov.rescale = FALSE,
                       auto.fix.first = FALSE, auto.fix.single = FALSE,
                       effect.coding = TRUE, meanstructure = TRUE)
  } else {
    tn    <- paste0("time", seq_len(TT))
    Sigma <- mom$Sigma_sum; dimnames(Sigma) <- list(tn, tn)
    mu    <- mom$mu_sum;     names(mu) <- tn
    fit <- lavaan::growth(lav_model, sample.cov = Sigma, sample.mean = mu,
                          sample.nobs = NOBS, sample.cov.rescale = FALSE)
  }
  pe <- lavaan::parameterEstimates(fit)
  gv <- function(l) { v <- pe$est[pe$label == l]; if (length(v)) v[1] else NA_real_ }
  ss <- lavaan::standardizedSolution(fit)
  cor12 <- { v <- ss$est.std[ss$label == "phi12"]; if (length(v)) v[1] else NA_real_ }
  fm <- tryCatch(lavaan::fitMeasures(fit, c("chisq", "df", "rmsea", "cfi", "srmr")),
                 error = function(e) rep(NA_real_, 5))
  list(phi1 = gv("phi1"), phi2 = gv("phi2"), phi12 = gv("phi12"),
       cor12 = cor12, k1 = gv("k1"), k2 = gv("k2"),
       fmin = lavaan::fitMeasures(fit, "fmin"),   # 0.5 * F_ml at optimum
       rmsea = unname(fm["rmsea"]), cfi = unname(fm["cfi"]), srmr = unname(fm["srmr"]),
       df = unname(fm["df"]))
}

## ---- design grid (32 unique parameter cells; plim is N-invariant) ----
grid <- expand.grid(
  rho            = c(0, 0.15),
  TT             = c(3, 5),
  J              = c(3, 5),
  loading_type   = c("tau_eq", "congeneric"),
  loading_quality= c("vgood", "poor"),
  stringsAsFactors = FALSE
)
grid$sd.lambda      <- ifelse(grid$loading_type == "tau_eq", 0, 0.05)
grid$target_loading <- ifelse(grid$loading_quality == "vgood", 0.6708, 0.4)
grid$cu_modeled     <- with(grid, rho != 0 & TT != 3)
grid$spec_status    <- with(grid, ifelse(
  rho == 0, "correct (no CU)",
  ifelse(TT != 3, "correct (CU modeled)", "MISSPEC (CU omitted)")))

## ---- run ----
out <- vector("list", nrow(grid))
for (r in seq_len(nrow(grid))) {
  g  <- grid[r, ]
  pr <- build_params(g$TT, g$J, g$sd.lambda, g$target_loading)
  mom <- pop_moments(pr$phi.mtx, c(MU_ALPHA, MU_BETA), PSI_MAT,
                     pr$lam_mat, pr$tau_mat, pr$theta_mat,
                     rho = g$rho, TT = g$TT, J = g$J, b = NULL)

  ## analytic EC-fitted target (phi_ec) for this cell
  ec_pop <- theta_from_theta_gen(pr$phi.mtx, c(MU_ALPHA, MU_BETA),
                                 pr$lam_vec, pr$tau_mat, b = NULL)
  tgt <- list(phi1 = ec_pop$covs[1, 1], phi2 = ec_pop$covs[2, 2],
              phi12 = ec_pop$covs[2, 1],
              cor12 = ec_pop$covs[2, 1] / sqrt(ec_pop$covs[1, 1] * ec_pop$covs[2, 2]),
              k1 = ec_pop$means[1], k2 = ec_pop$means[2])

  ec  <- fit_plim("EC_latent", build_syntax("EC_latent", g$TT, g$J, g$rho), mom, g$TT, g$J)
  sm  <- fit_plim("sum",       build_syntax("sum",       g$TT, g$J, g$rho), mom, g$TT, g$J)

  mk <- function(model, p) data.frame(
    model = model, rho = g$rho, TT = g$TT, J = g$J,
    loading_type = g$loading_type, loading_quality = g$loading_quality,
    cu_modeled = g$cu_modeled, spec_status = g$spec_status,
    # plim point estimates
    phi1.plim = p$phi1, phi2.plim = p$phi2, phi12.plim = p$phi12,
    cor12.plim = p$cor12, k1.plim = p$k1, k2.plim = p$k2,
    # analytic target (phi_ec)
    phi1.tgt = tgt$phi1, phi2.tgt = tgt$phi2, phi12.tgt = tgt$phi12,
    cor12.tgt = tgt$cor12, k1.tgt = tgt$k1, k2.tgt = tgt$k2,
    # asymptotic relative bias = (plim - target)/target
    phi1.arbias  = (p$phi1  - tgt$phi1 ) / tgt$phi1,
    phi2.arbias  = (p$phi2  - tgt$phi2 ) / tgt$phi2,
    phi12.arbias = (p$phi12 - tgt$phi12) / tgt$phi12,
    cor12.arbias = (p$cor12 - tgt$cor12) / tgt$cor12,
    k2.arbias    = (p$k2    - tgt$k2   ) / tgt$k2,
    # population misfit of the fitted model
    rmsea.pop = p$rmsea, cfi.pop = p$cfi, srmr.pop = p$srmr, df.pop = p$df,
    stringsAsFactors = FALSE)

  out[[r]] <- rbind(mk("EC_latent", ec), mk("sum", sm))
}

plim_table <- do.call(rbind, out)
rownames(plim_table) <- NULL
write.csv(plim_table, "study1_plim_T35.csv", row.names = FALSE)

## ---- console summary: asymptotic relative bias by model x TT x rho ----
cat("\n=== Asymptotic relative bias (plim - phi_ec)/phi_ec, averaged over J/quality/type ===\n")
agg <- aggregate(cbind(phi2.arbias, phi12.arbias, cor12.arbias) ~ model + TT + rho + spec_status,
                 data = plim_table, FUN = function(x) round(mean(x), 4))
print(agg)

cat("\n=== Max |asymptotic rel. bias| of phi2 in CORRECTLY specified cells (should be ~0) ===\n")
ok <- plim_table[plim_table$spec_status != "MISSPEC (CU omitted)", ]
print(tapply(abs(ok$phi2.arbias), ok$model, max))

cat("\nWrote study1_plim_T35.csv (", nrow(plim_table), "rows )\n")
