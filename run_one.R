
args   <- commandArgs(trailingOnly = TRUE)
job_id <- as.integer(args[1])
if (is.na(job_id)) stop("Usage: Rscript run_one.R <job_id>")

# ── Project path (set in submit.sh via PROJECT_DIR env var) ─────────────────
project_dir <- Sys.getenv("PROJECT_DIR", unset = ".")
setwd(project_dir)

source("mean_function.R")
library(lavaan)
library(MASS)
library(dplyr)

# ── Load grid row ────────────────────────────────────────────────────────────
cond <- read.csv("conditions.csv", stringsAsFactors = FALSE)
row  <- cond[cond$job_id == job_id, ]
if (nrow(row) == 0) stop("job_id not found: ", job_id)

cat(sprintf(
  "[job %02d] model=%-9s N=%-3d rho=%-4.2f",
  row$job_id,
  row$model,
  row$N_size,
  row$rho_val
))

# ── Run simulation ───────────────────────────────────────────────────────────
result <- run_vary_N(
  nrep         = row$nrep,
  N_vec        = row$N_size,     
  model        = row$model,
  var_alpha    = row$var_alpha,
  var_beta     = row$var_beta,
  cov_ab       = row$cov_ab,
  mu_alpha     = row$mu_alpha,
  mu_beta      = row$mu_beta,
  rho          = row$rho_val,
  psi          = row$psi,
  mni_location = row$mni_location,
  mni_size     = row$mni_size,
  seed         = row$seed
)

# ── Save output ──────────────────────────────────────────────────────────────
dir.create("results", showWarnings = FALSE)

out_file <- sprintf(
  "results/job%02d_%s_N%d_rho.RData",
  row$job_id,
  row$model,
  row$N_size,
  gsub("\\.", "p", as.character(row$rho_val))
)

save(result, row, file = out_file)
cat(sprintf("[job %02d] saved: %s\n", row$job_id, out_file))