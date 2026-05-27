library(dplyr)
library(tidyr)
library(purrr)
library(stringr)
library(ggplot2)
library(forcats)

load("~/Documents/GitHub/SOG_SIM/Qual_Code/1217/5i5t_ST_basic3.RData")
load("~/Documents/GitHub/SOG_SIM/Qual_Code/1217/rhos/5i5t_ST_rho3.RData")

want_cols <- c(
  "phi11.power","phi12.power","phi22.power","kappa1.power","kappa2.power",
  "phi11.R.bias","phi12.R.bias","phi22.R.bias","kappa1.R.bias","kappa2.R.bias",
  "phi11.rel.rmse","phi12.rel.rmse","phi22.rel.rmse","kappa1.rel.rmse","kappa2.rel.rmse",
  "phi11.coverage","phi12.coverage","phi22.coverage","kappa1.coverage","kappa2.coverage",
  "na.perc","admissible"  # convergence rate
)

obj_names <- ls(pattern = "^(latent|composite)_.+_ST_\\d+_(basic|rho|sd)$")


one_to_long <- function(nm) {
  x <- get(nm, envir = .GlobalEnv)
  sm <- x[["summary"]]
  
  if (is.null(sm)) {
    warning("No [[summary]] in: ", nm)
    return(NULL)
  }
  
  sm <- as.data.frame(sm)
  if (nrow(sm) > 1) sm <- sm[1, , drop = FALSE]  
  
  # --- composite: alpha -> kappa ---
  cn <- colnames(sm)
  cn <- str_replace(cn, "^alpha1\\b", "kappa1")
  cn <- str_replace(cn, "^alpha2\\b", "kappa2")
  colnames(sm) <- cn
  
  # --- keep relevant ones only ---
  keep <- intersect(want_cols, colnames(sm))
  sm <- sm[, keep, drop = FALSE]
  
  # --- name parcing ---
  m <- str_match(
    nm,
    "^(latent|composite)_([^_]+)_([^_]+)_([^_]+)_ST_(\\d+)_(basic|rho|sd)$"
  )
  if (anyNA(m)) {
    warning("Name does not match expected pattern: ", nm)
    return(NULL)
  }
  
  meta <- tibble(
    object_name = nm,
    model      = m[,2],
    design     = m[,3],
    slope_var  = m[,4],
    quality    = m[,5],
    nsize      = as.integer(m[,6]),
    type       = m[,7]
  )
  
  # --- into long ---
  special_keys <- c("na.perc", "admissible")
  
  out <- sm %>%
    tibble::as_tibble() %>%
    pivot_longer(cols = everything(), names_to = "key", values_to = "value") %>%
    mutate(
      key = if_else(key %in% special_keys, "na.perc", key),   # <- 표준화
      param  = if_else(key == "na.perc", "convergence", str_extract(key, "^[^.]+")),
      metric = if_else(key == "na.perc", "cov.rate",   str_extract(key, "(?<=\\.).+$"))
    ) %>%
    select(-key) %>%
    bind_cols(meta, .)
  
  out
}

# final summary
sim_long_orig <- map_dfr(obj_names, one_to_long) %>%
  select(model, type, design, slope_var, quality, nsize, param, metric, value, object_name) %>%
  arrange(model, type, slope_var, quality, nsize, param, metric)



##### Plotting ======

sim_long <- sim_long_orig %>%
  mutate(quality = recode(quality, "poor" = "Poor", "good" = "Good", "vgood" = "Very Good"),
         metric= recode(metric, "rel.rmse" = "R.RMSE", "coverage"= "95% Coverage", "power" = "Power",
                        "cov.rate" ="Convergence"),
         slope_var = recode(slope_var, "small" = "Small", "medium" = "Medium","large" = "Large"),
         type = recode(type, "basic" = "Basic", "rho" = "Cor.Unique"))

sim_long$quality <- factor(sim_long$quality, levels = c("Poor", "Good", "Very Good"))
sim_long$metric <- factor(sim_long$metric, levels = c("R.RMSE", "R.bias", "95% Coverage", 
                                                      "Power", "Convergence"))
sim_long$slope_var <- factor(sim_long$slope_var, levels = c("Small", "Medium", "Large"))
sim_long$type <- factor(sim_long$type, levels = c("Basic", "Cor.Unique"))


hline_data_fixed <- tibble(
  metric = c("Power", "R.bias", "R.bias", "95% Coverage"),
  yint   = c(0.8,     0.1,      -0.1,     0.8)
) %>%
  crossing(
    slope_var = levels(sim_long$slope_var),
    type      = unique(sim_long$type),
    quality   = levels(sim_long$quality)
  ) %>%
  mutate(
    quality   = factor(quality,   levels = levels(sim_long$quality)),
    type      = factor(type,      levels = levels(sim_long$type)),
    slope_var = factor(slope_var, levels = levels(sim_long$slope_var)),
    metric    = factor(metric,    levels = levels(sim_long$metric))
  )




##### phi11 #####
pdat1 <- sim_long %>%
  filter(metric %in% c("Power", "95% Coverage"),
         param  %in% c("phi11")) %>%
  mutate(
    quality = fct_relevel(as.factor(quality), "Poor", "Good", "Very Good"),
    type    = fct_relevel(as.factor(type), "Basic", "Cor.Unique")   # 필요하면
  )

ggplot(pdat1, aes(x
                  = factor(nsize), y = value, color = model, linetype = model,
                  group = interaction(model, quality))) +
  geom_point(position = position_dodge(width = .3)) +
  geom_line(position = position_dodge(width = .3)) +
  geom_hline(
    data =  hline_data_fixed %>% filter(metric %in% c("Power", "95% Coverage")),
    aes(yintercept = yint),
    inherit.aes = FALSE,
    linetype = "dashed", color = "gray"
  ) +
  facet_grid(metric + slope_var ~ type + quality, scales = "free_y") +
  scale_linetype_manual(values = c("composite" = "twodash", "latent" = "solid")) +
  ggtitle("Variance of Growth Intercept") +
  labs(x = "Sample size (n)", y = "Value", color = "Model", linetype = "Model") + 
  theme_bw(base_family = "Times New Roman") +
  theme(text = element_text(family = "Times New Roman"))


pdat2 <- sim_long %>%
  filter(metric %in% c("R.RMSE", "R.bias"),
         param  %in% c("phi11")) %>%
  mutate(
    quality = fct_relevel(as.factor(quality), "Poor", "Good", "Very Good"),
    type    = fct_relevel(as.factor(type), "Basic", "Cor.Unique")   # 필요하면
  )

ggplot(pdat2, aes(x
                  = factor(nsize), y = value, color = model, linetype = model,
                  group = interaction(model, quality))) +
  geom_point(position = position_dodge(width = .3)) +
  geom_line(position = position_dodge(width = .3)) +
  geom_hline(
    data =  hline_data_fixed %>% filter(metric %in% c("R.RMSE", "R.bias")),
    aes(yintercept = yint),
    inherit.aes = FALSE,
    linetype = "dashed", color = "gray"
  ) +
  facet_grid(metric + slope_var ~ type + quality, scales = "free_y") +
  scale_linetype_manual(values = c("composite" = "twodash", "latent" = "solid")) +
  ggtitle("Variance of Growth Intercept") +
  labs(x = "Sample size (n)", y = "Value", color = "Model", linetype = "Model") +
  theme_bw(base_family = "Times New Roman") +
  theme(text = element_text(family = "Times New Roman"))


##### phi12 #####
phi12_power <- sim_long %>%
  filter(metric %in% c("Power", "95% Coverage"),
         param  %in% c("phi12")) %>%
  mutate(
    quality = fct_relevel(as.factor(quality), "Poor", "Good", "Very Good"),
    type    = fct_relevel(as.factor(type), "Basic", "Cor.Unique")   # 필요하면
  )

ggplot(phi12_power, aes(x
                  = factor(nsize), y = value, color = model, linetype = model,
                  group = interaction(model, quality))) +
  geom_point(position = position_dodge(width = .3)) +
  geom_line(position = position_dodge(width = .3)) +
  geom_hline(
    data =  hline_data_fixed %>% filter(metric %in% c("Power", "95% Coverage")),
    aes(yintercept = yint),
    inherit.aes = FALSE,
    linetype = "dashed", color = "gray"
  ) +
  facet_grid(metric + slope_var ~ type + quality, scales = "free_y") +
  scale_linetype_manual(values = c("composite" = "twodash", "latent" = "solid")) +
  ggtitle("Covariance between Growth Intercept & Slope") +
  labs(x = "Sample size (n)", y = "Value", color = "Model", linetype = "Model") +
  theme_bw(base_family = "Times New Roman") +
  theme(text = element_text(family = "Times New Roman"))


phi12_accu <- sim_long %>%
  filter(metric %in% c("R.RMSE", "R.bias"),
         param  %in% c("phi12")) %>%
  mutate(
    quality = fct_relevel(as.factor(quality), "Poor", "Good", "Very Good"),
    type    = fct_relevel(as.factor(type), "Basic", "Cor.Unique")   # 필요하면
  )

ggplot(phi12_accu, aes(x
                  = factor(nsize), y = value, color = model, linetype = model,
                  group = interaction(model, quality))) +
  geom_point(position = position_dodge(width = .3)) +
  geom_line(position = position_dodge(width = .3)) +
  geom_hline(
    data =  hline_data_fixed %>% filter(metric %in% c("R.RMSE", "R.bias")),
    aes(yintercept = yint),
    inherit.aes = FALSE,
    linetype = "dashed", color = "gray"
  ) +
  facet_grid(metric + slope_var ~ type + quality, scales = "free_y") +
  scale_linetype_manual(values = c("composite" = "twodash", "latent" = "solid")) +
  ggtitle("Covariance between Growth Intercept & Slope") +
  labs(x = "Sample size (n)", y = "Value", color = "Model", linetype = "Model") +
  theme_bw(base_family = "Times New Roman") +
  theme(text = element_text(family = "Times New Roman"))


##### phi22 #####
phi22_power <- sim_long %>%
  filter(metric %in% c("Power", "95% Coverage"),
         param  %in% c("phi22")) %>%
  mutate(
    quality = fct_relevel(as.factor(quality), "Poor", "Good", "Very Good"),
    type    = fct_relevel(as.factor(type), "Basic", "Cor.Unique")   # 필요하면
  )

ggplot(phi22_power, aes(x
                        = factor(nsize), y = value, color = model, linetype = model,
                        group = interaction(model, quality))) +
  geom_point(position = position_dodge(width = .3)) +
  geom_line(position = position_dodge(width = .3)) +
  geom_hline(
    data =  hline_data_fixed %>% filter(metric %in% c("Power", "95% Coverage")),
    aes(yintercept = yint),
    inherit.aes = FALSE,
    linetype = "dashed", color = "gray"
  ) +
  facet_grid(metric + slope_var ~ type + quality, scales = "free_y") +
  scale_linetype_manual(values = c("composite" = "twodash", "latent" = "solid")) + 
  ggtitle("Variance of Growth Slope") +
  labs(x = "Sample size (n)", y = "Value", color = "Model", linetype = "Model") +
  theme_bw(base_family = "Times New Roman") +
  theme(text = element_text(family = "Times New Roman"))


phi22_accu <- sim_long %>%
  filter(metric %in% c("R.RMSE", "R.bias"),
         param  %in% c("phi22")) %>%
  mutate(
    quality = fct_relevel(as.factor(quality), "Poor", "Good", "Very Good"),
    type    = fct_relevel(as.factor(type), "Basic", "Cor.Unique")   # 필요하면
  )

ggplot(phi22_accu, aes(x
                       = factor(nsize), y = value, color = model, linetype = model,
                       group = interaction(model, quality))) +
  geom_point(position = position_dodge(width = .3)) +
  geom_line(position = position_dodge(width = .3)) +
  geom_hline(
    data =  hline_data_fixed %>% filter(metric %in% c("R.RMSE", "R.bias")),
    aes(yintercept = yint),
    inherit.aes = FALSE,
    linetype = "dashed", color = "gray"
  ) +
  facet_grid(metric + slope_var ~ type + quality, scales = "free_y") +
  scale_linetype_manual(values = c("composite" = "twodash", "latent" = "solid")) +
  ggtitle("Variance of Growth Slope") +
  labs(x = "Sample size (n)", y = "Value", color = "Model", linetype = "Model") +
  theme_bw(base_family = "Times New Roman") +
  theme(text = element_text(family = "Times New Roman"))





##### kappa1 #####
k_power <- sim_long %>%
  filter(metric %in% c("Power", "95% Coverage"),
         param  %in% c("kappa1","kappa2")) %>%
  mutate(
    quality = fct_relevel(as.factor(quality), "Poor", "Good", "Very Good"),
    type    = fct_relevel(as.factor(type), "Basic", "Cor.Unique")   # 필요하면
  )

ggplot(k_power, aes(x
                        = factor(nsize), y = value, color = model, linetype = model,
                        group = interaction(model, quality))) +
  geom_point(position = position_dodge(width = .3)) +
  geom_line(position = position_dodge(width = .3)) +
  geom_hline(
    data =  hline_data_fixed %>% filter(metric %in% c("Power", "95% Coverage")),
    aes(yintercept = yint),
    inherit.aes = FALSE,
    linetype = "dashed", color = "gray"
  ) +
  facet_grid(metric + slope_var ~ param + type + quality, scales = "free_y") +
  scale_linetype_manual(values = c("composite" = "twodash", "latent" = "solid")) +
  ggtitle("Mean of Growth Intercept & Slope") +
  labs(x = "Sample size (n)", y = "Value", color = "Model", linetype = "Model") +
  theme_bw(base_family = "Times New Roman") +
  theme(text = element_text(family = "Times New Roman"))


k_accu <- sim_long %>%
  filter(metric %in% c("R.RMSE", "R.bias"),
         param  %in% c("kappa1","kappa2")) %>%
  mutate(
    quality = fct_relevel(as.factor(quality), "Poor", "Good", "Very Good"),
    type    = fct_relevel(as.factor(type), "Basic", "Cor.Unique")   # 필요하면
  )

ggplot(k_accu, aes(x
                       = factor(nsize), y = value, color = model, linetype = model,
                       group = interaction(model, quality))) +
  geom_point(position = position_dodge(width = .3)) +
  geom_line(position = position_dodge(width = .3)) +
  geom_hline(
    data =  hline_data_fixed %>% filter(metric %in% c("R.RMSE", "R.bias")),
    aes(yintercept = yint),
    inherit.aes = FALSE,
    linetype = "dashed", color = "gray"
  ) +
  facet_grid(metric + slope_var ~ param + type + quality, scales = "free_y") +
  scale_linetype_manual(values = c("composite" = "twodash", "latent" = "solid")) +
  ggtitle("Mean of Growth Intercept & Slope") +
  labs(x = "Sample size (n)", y = "Value", color = "Model", linetype = "Model") +
  theme_bw(base_family = "Times New Roman") +
  theme(text = element_text(family = "Times New Roman"))


# convergence
conv <- sim_long %>%
  filter(metric %in% c("Convergence"),
         param  %in% c("convergence")) %>%
  mutate(
    quality = fct_relevel(as.factor(quality), "Poor", "Good", "Very Good"),
    type    = fct_relevel(as.factor(type), "Basic", "Cor.Unique")   # 필요하면
  )

ggplot(conv, aes(x
                   = factor(nsize), y = value, color = model, linetype = model,
                   group = interaction(model, quality))) +
  geom_point(position = position_dodge(width = .3)) +
  geom_line(position = position_dodge(width = .3)) +
  geom_hline(
    data =  hline_data_fixed %>% filter(metric %in% c("Convergence")),
    aes(yintercept = yint),
    inherit.aes = FALSE,
    linetype = "dashed", color = "gray"
  ) +
  facet_grid(slope_var ~ type + quality, scales = "free_y") +
  scale_linetype_manual(values = c("composite" = "twodash", "latent" = "solid")) +
  ggtitle("Admissible Output Proportion") +
  labs(x = "Sample size (n)", y = "Value", color = "Model", linetype = "Model") +
  theme_bw(base_family = "Times New Roman") +
  theme(text = element_text(family = "Times New Roman"))

