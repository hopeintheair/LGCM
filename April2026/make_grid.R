# Grid generator for function_mine.R::run()
# Varies: model x (mni_location x mni_size) x N
# Other arguments use function defaults (set in run_one).

sample_size <- data.frame(
  N      = c("small", "medium", "large"),
  N_size = c(200, 500, 800),
  stringsAsFactors = FALSE
)

# Valid (location, size) combinations:
# "none" only pairs with "none"; otherwise weak/strong on each location.
mni_grid <- rbind(
  data.frame(mni_location = "none",      mni_size = "none"),
  expand.grid(
    mni_location = c("intercept", "loading", "residual"),
    mni_size     = c("weak", "strong"),
    stringsAsFactors = FALSE
  )
)
mni_grid$mni_location <- as.character(mni_grid$mni_location)
mni_grid$mni_size     <- as.character(mni_grid$mni_size)

models <- c("latent", "mean", "sum", "mean_scaled", "sum_scaled")

loading_grid <- data.frame(
  loading_quality = c("vgood", "good", "poor"),
  target_loading  = c(0.7, 0.5, 0.3),
  stringsAsFactors = FALSE
)

conditions <- expand.grid(
  model        = models,
  mni_idx      = seq_len(nrow(mni_grid)),
  N_size       = sample_size$N_size,
  rho          = c(0, 0.15),
  loading_quality = loading_grid$loading_quality,
  stringsAsFactors = FALSE
)
conditions$mni_location <- mni_grid$mni_location[conditions$mni_idx]
conditions$mni_size     <- mni_grid$mni_size[conditions$mni_idx]
conditions$mni_idx      <- NULL

conditions <- merge(conditions, sample_size, by = "N_size", sort = FALSE)
conditions <- merge(conditions, loading_grid, by = "loading_quality", sort = FALSE)

# Fixed parameters (function_mine.R defaults; override here as needed)
conditions$TT             <- 5
conditions$J              <- 5
conditions$mu_alpha       <- 1
conditions$mu_beta        <- 0.125
conditions$var_alpha      <- 1
conditions$var_beta       <- 0.2025
conditions$cov_ab         <- 0.18
conditions$tau_un         <- 1
conditions$lambda_un      <- 1
conditions$sd.lambda      <- 0
conditions$nrep           <- 2000
conditions$seed           <- 218

# Sort and assign job_id
conditions$loading_quality <- factor(conditions$loading_quality,
                                     levels = c("vgood", "good", "poor"))
conditions <- conditions[order(conditions$model,
                               conditions$mni_location,
                               conditions$mni_size,
                               conditions$rho,
                               conditions$loading_quality,
                               conditions$N_size), ]
conditions$loading_quality <- as.character(conditions$loading_quality)
conditions$job_id <- seq_len(nrow(conditions))

# Reorder columns: job_id first
front <- c("job_id", "model", "N", "N_size",
           "mni_location", "mni_size", "rho",
           "loading_quality", "target_loading")
conditions <- conditions[, c(front, setdiff(names(conditions), front))]

write.csv(conditions, "conditions_new.csv", row.names = FALSE)
cat("Created conditions_new.csv with", nrow(conditions), "rows.\n")
