
sample_size <- data.frame(
  N     = c("small", "medium", "large"),
  N_size = c(200, 500, 800),
  stringsAsFactors = FALSE
)


rho_levels <- data.frame(
  residual_cov  = c("no rho",  "small rho"),
  rho_val       = c(0, 0.1),
  stringsAsFactors = FALSE
)


conditions <- expand.grid(
  model    = c("latent", "composite"),
  rho_val  = c(0, 0.1),
  N_size = c(200, 500, 800),
  stringsAsFactors = FALSE
)

conditions <- merge(conditions, sample_size,  by = "N_size")
conditions <- merge(conditions, rho_levels,  by = "rho_val")

# Fixed parameters
conditions$var_alpha <- 1
conditions$var_beta <- 0.2
conditions$cov_ab    <- 0.1
conditions$mu_alpha  <- 0
conditions$mu_beta   <- 0.5
conditions$psi       <- 0.2
conditions$mni_location <- "none"
conditions$mni_size <- "none"
conditions$nrep      <- 2000
conditions$seed      <- 218

# Sort and assign job_id
conditions <- conditions[order(conditions$model, conditions$N_size,
                                conditions$rho_val), ]
conditions$job_id <- seq_len(nrow(conditions))

write.csv(conditions, "conditions.csv", row.names = FALSE)
cat("Created conditions.csv with", nrow(conditions), "rows.\n")
