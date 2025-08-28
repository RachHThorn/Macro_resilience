# R Thornley
# 27/08/2025
# Project: P1_COMPADRE_DRAGNET
# get predictions from the GLMMs

# install the mixed up package directly from Github
# remotes::install_github('m-clark/mixedup')

# load packages
packages <- c("tidyverse", "glmmTMB", "rstan", "assertthat", "purrr", "remotes", "mixedup")

for (p in packages) {
  if (!requireNamespace(p, quietly = TRUE)) {
    install.packages(p)
  }
  library(p, character.only = TRUE)
}

################################################################################

data <- read_csv("results/DRAGNet_T0_T1_all.csv")
experiments <- c("Control", "Disturbance")
data <- data %>% filter(trt %in% experiments)
names(data)
betareg_mod <- glmmTMB(
  new_max_cover ~ trt * year_trt + (1|New_taxon),
  family = ordbeta(link = "logit"),
  data = data)

data$predicted <- predict(betareg_mod)

# Calculate mean x per group
group_means <- data %>%
  group_by(New_taxon, trt, year_trt) %>%
  summarise(new_max_cover = mean(new_max_cover), .groups = "drop")
group_means$predicted <- predict(betareg_mod, newdata = group_means, re.form = NULL)

# random effect deviations from the population mean
# the use these to preditc for each level of the REs
fixef(betareg_mod)["(Intercept)"] + ranef(betareg_mod)$New_taxon


################################################################################

# Load and filter
data <- read_csv("results/DRAGNet_T0_T1_all.csv")
experiments <- c("Control", "Disturbance")
data <- data %>% filter(trt %in% experiments)

# Fit model
betareg_mod <- glmmTMB(
  new_max_cover ~ trt * year_trt + (1 | New_taxon),
  family = ordbeta(link = "logit"),
  data = data
)

# check if the random intercepts vary between models
VarCorr(betareg_mod)

# Predictions for each observation (fixed + random effects)
data$predicted <- predict(betareg_mod)

# Group-level means and predictions
group_means <- data %>%
  group_by(New_taxon, trt, year_trt) %>%
  summarise(new_max_cover = mean(new_max_cover), .groups = "drop")

group_means$predicted <- predict(betareg_mod, newdata = group_means, re.form = NULL)

# Extract random effect deviations
random_effects <- ranef(betareg_mod)$New_taxon

# Combine fixed intercept with random deviations
intercepts <- fixef(betareg_mod)$cond["(Intercept)"] + random_effects[, "(Intercept)"]
intercepts

coef(betareg_mod)$New_taxon


# Fixed intercept
fixed_int <- fixef(betareg_mod)$cond["(Intercept)"]

# Random intercept deviations
re <- ranef(betareg_mod)$New_taxon[["(Intercept)"]]

# Combine fixed + random to get group-specific intercepts
taxon_intercepts <- fixed_int + re

head(taxon_intercepts)


ranef(betareg_mod)$New_taxon
VarCorr(betareg_mod)
