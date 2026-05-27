############################################################
## I = 5, tlen = 5  (run5i5t) ###############################
############################################################

source("calculation_qual.R")
source("analysis_5i5t_basic.R")
library(lavaan)

############################################################
## small variance (phi2 = 0.32, phi12 = 0.13)
############################################################

## ---- small / poor (sigma_mu = 1.57) ----
latent_5i5_small_poor_ST_200_basic <- run5i5t(
  nrep = 1000, model = "latent",
  k1 = 1, k2 = 0.125,
  phi1 = 1, phi2 = 0.32, phi12 = 0.13,
  sigma_mu = 1.57,
  psi1 = 1.097, psi2 = .98, psi3 = .95, psi4 = .93, psi5 = .91,
  inv.type = c("strict"),
  tau = .10, lambda = .40, sd.lambda = 0.0,
  nsize = 200, auto_var1 = FALSE, rho = 0, seed = 3201
)
composite_5i5_small_poor_ST_200_basic <- run5i5t(
  nrep = 1000, model = "composite",
  k1 = 1, k2 = 0.125,
  phi1 = 1, phi2 = 0.32, phi12 = 0.13,
  sigma_mu = 1.57,
  psi1 = 1.097, psi2 = .98, psi3 = .95, psi4 = .93, psi5 = .91,
  inv.type = c("strict"),
  tau = .10, lambda = .40, sd.lambda = 0.0,
  nsize = 200, auto_var1 = FALSE, rho = 0, seed = 3201
)

latent_5i5_small_poor_ST_500_basic <- run5i5t(
  nrep = 1000, model = "latent",
  k1 = 1, k2 = 0.125,
  phi1 = 1, phi2 = 0.32, phi12 = 0.13,
  sigma_mu = 1.57,
  psi1 = 1.097, psi2 = .98, psi3 = .95, psi4 = .93, psi5 = .91,
  inv.type = c("strict"),
  tau = .10, lambda = .40, sd.lambda = 0.0,
  nsize = 500, auto_var1 = FALSE, rho = 0, seed = 3201
)
composite_5i5_small_poor_ST_500_basic <- run5i5t(
  nrep = 1000, model = "composite",
  k1 = 1, k2 = 0.125,
  phi1 = 1, phi2 = 0.32, phi12 = 0.13,
  sigma_mu = 1.57,
  psi1 = 1.097, psi2 = .98, psi3 = .95, psi4 = .93, psi5 = .91,
  inv.type = c("strict"),
  tau = .10, lambda = .40, sd.lambda = 0.0,
  nsize = 500, auto_var1 = FALSE, rho = 0, seed = 3201
)

latent_5i5_small_poor_ST_800_basic <- run5i5t(
  nrep = 1000, model = "latent",
  k1 = 1, k2 = 0.125,
  phi1 = 1, phi2 = 0.32, phi12 = 0.13,
  sigma_mu = 1.57,
  psi1 = 1.097, psi2 = .98, psi3 = .95, psi4 = .93, psi5 = .91,
  inv.type = c("strict"),
  tau = .10, lambda = .40, sd.lambda = 0.0,
  nsize = 800, auto_var1 = FALSE, rho = 0, seed = 3201
)
composite_5i5_small_poor_ST_800_basic <- run5i5t(
  nrep = 1000, model = "composite",
  k1 = 1, k2 = 0.125,
  phi1 = 1, phi2 = 0.32, phi12 = 0.13,
  sigma_mu = 1.57,
  psi1 = 1.097, psi2 = .98, psi3 = .95, psi4 = .93, psi5 = .91,
  inv.type = c("strict"),
  tau = .10, lambda = .40, sd.lambda = 0.0,
  nsize = 800, auto_var1 = FALSE, rho = 0, seed = 3201
)



## ---- small / vgood (sigma_mu = 0.36) ----
latent_5i5_small_vgood_ST_200_basic <- run5i5t(
  nrep = 1000, model = "latent",
  k1 = 1, k2 = 0.125,
  phi1 = 1, phi2 = 0.32, phi12 = 0.13,
  sigma_mu = 0.36,
  psi1 = 1.097, psi2 = .98, psi3 = .95, psi4 = .93, psi5 = .91,
  inv.type = c("strict"),
  tau = .10, lambda = .40, sd.lambda = 0.0,
  nsize = 200, auto_var1 = FALSE, rho = 0, seed = 3201
)
composite_5i5_small_vgood_ST_200_basic <- run5i5t(
  nrep = 1000, model = "composite",
  k1 = 1, k2 = 0.125,
  phi1 = 1, phi2 = 0.32, phi12 = 0.13,
  sigma_mu = 0.36,
  psi1 = 1.097, psi2 = .98, psi3 = .95, psi4 = .93, psi5 = .91,
  inv.type = c("strict"),
  tau = .10, lambda = .40, sd.lambda = 0.0,
  nsize = 200, auto_var1 = FALSE, rho = 0, seed = 3201
)

latent_5i5_small_vgood_ST_500_basic <- run5i5t(
  nrep = 1000, model = "latent",
  k1 = 1, k2 = 0.125,
  phi1 = 1, phi2 = 0.32, phi12 = 0.13,
  sigma_mu = 0.36,
  psi1 = 1.097, psi2 = .98, psi3 = .95, psi4 = .93, psi5 = .91,
  inv.type = c("strict"),
  tau = .10, lambda = .40, sd.lambda = 0.0,
  nsize = 500, auto_var1 = FALSE, rho = 0, seed = 3201
)
composite_5i5_small_vgood_ST_500_basic <- run5i5t(
  nrep = 1000, model = "composite",
  k1 = 1, k2 = 0.125,
  phi1 = 1, phi2 = 0.32, phi12 = 0.13,
  sigma_mu = 0.36,
  psi1 = 1.097, psi2 = .98, psi3 = .95, psi4 = .93, psi5 = .91,
  inv.type = c("strict"),
  tau = .10, lambda = .40, sd.lambda = 0.0,
  nsize = 500, auto_var1 = FALSE, rho = 0, seed = 3201
)

latent_5i5_small_vgood_ST_800_basic <- run5i5t(
  nrep = 1000, model = "latent",
  k1 = 1, k2 = 0.125,
  phi1 = 1, phi2 = 0.32, phi12 = 0.13,
  sigma_mu = 0.36,
  psi1 = 1.097, psi2 = .98, psi3 = .95, psi4 = .93, psi5 = .91,
  inv.type = c("strict"),
  tau = .10, lambda = .40, sd.lambda = 0.0,
  nsize = 800, auto_var1 = FALSE, rho = 0, seed = 3201
)
composite_5i5_small_vgood_ST_800_basic <- run5i5t(
  nrep = 1000, model = "composite",
  k1 = 1, k2 = 0.125,
  phi1 = 1, phi2 = 0.32, phi12 = 0.13,
  sigma_mu = 0.36,
  psi1 = 1.097, psi2 = .98, psi3 = .95, psi4 = .93, psi5 = .91,
  inv.type = c("strict"),
  tau = .10, lambda = .40, sd.lambda = 0.0,
  nsize = 800, auto_var1 = FALSE, rho = 0, seed = 3201
)



############################################################
## large variance (phi2 = 0.55, phi12 = 0.22)
############################################################

## ---- large / poor (sigma_mu = 1.57) ----
latent_5i5_large_poor_ST_200_basic <- run5i5t(
  nrep = 1000, model = "latent",
  k1 = 1, k2 = 0.125,
  phi1 = 1, phi2 = 0.55, phi12 = 0.22,
  sigma_mu = 1.57,
  psi1 = 1.097, psi2 = .98, psi3 = .95, psi4 = .93, psi5 = .91,
  inv.type = c("strict"),
  tau = .10, lambda = .40, sd.lambda = 0.0,
  nsize = 200, auto_var1 = FALSE, rho = 0, seed = 3201
)
composite_5i5_large_poor_ST_200_basic <- run5i5t(
  nrep = 1000, model = "composite",
  k1 = 1, k2 = 0.125,
  phi1 = 1, phi2 = 0.55, phi12 = 0.22,
  sigma_mu = 1.57,
  psi1 = 1.097, psi2 = .98, psi3 = .95, psi4 = .93, psi5 = .91,
  inv.type = c("strict"),
  tau = .10, lambda = .40, sd.lambda = 0.0,
  nsize = 200, auto_var1 = FALSE, rho = 0, seed = 3201
)

latent_5i5_large_poor_ST_500_basic <- run5i5t(
  nrep = 1000, model = "latent",
  k1 = 1, k2 = 0.125,
  phi1 = 1, phi2 = 0.55, phi12 = 0.22,
  sigma_mu = 1.57,
  psi1 = 1.097, psi2 = .98, psi3 = .95, psi4 = .93, psi5 = .91,
  inv.type = c("strict"),
  tau = .10, lambda = .40, sd.lambda = 0.0,
  nsize = 500, auto_var1 = FALSE, rho = 0, seed = 3201
)
composite_5i5_large_poor_ST_500_basic <- run5i5t(
  nrep = 1000, model = "composite",
  k1 = 1, k2 = 0.125,
  phi1 = 1, phi2 = 0.55, phi12 = 0.22,
  sigma_mu = 1.57,
  psi1 = 1.097, psi2 = .98, psi3 = .95, psi4 = .93, psi5 = .91,
  inv.type = c("strict"),
  tau = .10, lambda = .40, sd.lambda = 0.0,
  nsize = 500, auto_var1 = FALSE, rho = 0, seed = 3201
)

latent_5i5_large_poor_ST_800_basic <- run5i5t(
  nrep = 1000, model = "latent",
  k1 = 1, k2 = 0.125,
  phi1 = 1, phi2 = 0.55, phi12 = 0.22,
  sigma_mu = 1.57,
  psi1 = 1.097, psi2 = .98, psi3 = .95, psi4 = .93, psi5 = .91,
  inv.type = c("strict"),
  tau = .10, lambda = .40, sd.lambda = 0.0,
  nsize = 800, auto_var1 = FALSE, rho = 0, seed = 3201
)
composite_5i5_large_poor_ST_800_basic <- run5i5t(
  nrep = 1000, model = "composite",
  k1 = 1, k2 = 0.125,
  phi1 = 1, phi2 = 0.55, phi12 = 0.22,
  sigma_mu = 1.57,
  psi1 = 1.097, psi2 = .98, psi3 = .95, psi4 = .93, psi5 = .91,
  inv.type = c("strict"),
  tau = .10, lambda = .40, sd.lambda = 0.0,
  nsize = 800, auto_var1 = FALSE, rho = 0, seed = 3201
)



## ---- large / vgood (sigma_mu = 0.36) ----
latent_5i5_large_vgood_ST_200_basic <- run5i5t(
  nrep = 1000, model = "latent",
  k1 = 1, k2 = 0.125,
  phi1 = 1, phi2 = 0.55, phi12 = 0.22,
  sigma_mu = 0.36,
  psi1 = 1.097, psi2 = .98, psi3 = .95, psi4 = .93, psi5 = .91,
  inv.type = c("strict"),
  tau = .10, lambda = .40, sd.lambda = 0.0,
  nsize = 200, auto_var1 = FALSE, rho = 0, seed = 3201
)
composite_5i5_large_vgood_ST_200_basic <- run5i5t(
  nrep = 1000, model = "composite",
  k1 = 1, k2 = 0.125,
  phi1 = 1, phi2 = 0.55, phi12 = 0.22,
  sigma_mu = 0.36,
  psi1 = 1.097, psi2 = .98, psi3 = .95, psi4 = .93, psi5 = .91,
  inv.type = c("strict"),
  tau = .10, lambda = .40, sd.lambda = 0.0,
  nsize = 200, auto_var1 = FALSE, rho = 0, seed = 3201
)

latent_5i5_large_vgood_ST_500_basic <- run5i5t(
  nrep = 1000, model = "latent",
  k1 = 1, k2 = 0.125,
  phi1 = 1, phi2 = 0.55, phi12 = 0.22,
  sigma_mu = 0.36,
  psi1 = 1.097, psi2 = .98, psi3 = .95, psi4 = .93, psi5 = .91,
  inv.type = c("strict"),
  tau = .10, lambda = .40, sd.lambda = 0.0,
  nsize = 500, auto_var1 = FALSE, rho = 0, seed = 3201
)
composite_5i5_large_vgood_ST_500_basic <- run5i5t(
  nrep = 1000, model = "composite",
  k1 = 1, k2 = 0.125,
  phi1 = 1, phi2 = 0.55, phi12 = 0.22,
  sigma_mu = 0.36,
  psi1 = 1.097, psi2 = .98, psi3 = .95, psi4 = .93, psi5 = .91,
  inv.type = c("strict"),
  tau = .10, lambda = .40, sd.lambda = 0.0,
  nsize = 500, auto_var1 = FALSE, rho = 0, seed = 3201
)

latent_5i5_large_vgood_ST_800_basic <- run5i5t(
  nrep = 1000, model = "latent",
  k1 = 1, k2 = 0.125,
  phi1 = 1, phi2 = 0.55, phi12 = 0.22,
  sigma_mu = 0.36,
  psi1 = 1.097, psi2 = .98, psi3 = .95, psi4 = .93, psi5 = .91,
  inv.type = c("strict"),
  tau = .10, lambda = .40, sd.lambda = 0.0,
  nsize = 800, auto_var1 = FALSE, rho = 0, seed = 3201
)
composite_5i5_large_vgood_ST_800_basic <- run5i5t(
  nrep = 1000, model = "composite",
  k1 = 1, k2 = 0.125,
  phi1 = 1, phi2 = 0.55, phi12 = 0.22,
  sigma_mu = 0.36,
  psi1 = 1.097, psi2 = .98, psi3 = .95, psi4 = .93, psi5 = .91,
  inv.type = c("strict"),
  tau = .10, lambda = .40, sd.lambda = 0.0,
  nsize = 800, auto_var1 = FALSE, rho = 0, seed = 3201
)

save.image(file = "5i5t_ST_basic3_2026.RData")
