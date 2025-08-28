# R Thornley
# 27/08/2025
# Project: P1_COMPADRE_DRAGNET
# get predictions from the GLMMs

# install the mixed up package directly from Github
# remotes::install_github('m-clark/mixedup')

# load packages
packages <- c("tidyverse", "glmmTMB", "rstan", "assertthat", "purrr", "remotes", 
              "mixedup", "easystats", "purrr")

for (p in packages) {
  if (!requireNamespace(p, quietly = TRUE)) {
    install.packages(p)
  }
  library(p, character.only = TRUE)
}

################################################################################

# test the different RE outputs as a function of 
data <- read_csv("results/DRAGNet_T0_T1_all.csv")
experiments <- c("Control", "Disturbance")
data <- data %>% filter(trt %in% experiments)
names(data)
data <- data %>% mutate(site_block = paste0(site_name, "_", block))
data <- data %>% mutate(site_block_taxon = paste0(site_name, "_", block, "_", New_taxon))
n_distinct(data$site_block) # there are 216 unique combinations of site and block
n_distinct(data$site_block_taxon) # there are 4002 unique combinations of site, block, taxon

data %>% group_by(site_block) %>% count()
data %>% group_by(site_block_taxon) %>% count()

# all these models converge 

# model 1: the model we have been using to date that just takes taxon as a RE (fixed effect is the BACI)
betareg_mod_1 <- glmmTMB(
  new_max_cover ~ trt * year_trt + (1|New_taxon),
  family = ordbeta(link = "logit"),
  data = data)

# model 2: addition of site_block as the RE
betareg_mod_2 <- glmmTMB(
  new_max_cover ~ trt * year_trt + (1|New_taxon) + (1|site_block),
  family = ordbeta(link = "logit"),
  data = data)

# model 3: addition of site_block_taxon as the RE (I think this is the right specification for pairing - repeat measures)
betareg_mod_3 <- glmmTMB(
  new_max_cover ~ trt * year_trt + (1|New_taxon) + (1|site_block_taxon),
  family = ordbeta(link = "logit"),
  data = data)

# model 4: specification of taxa as a fixed effect (wiht the repeat measures random effect)
betareg_mod_4 <- glmmTMB(
  new_max_cover ~ trt * year_trt + New_taxon + (1|site_block_taxon),
  family = ordbeta(link = "logit"),
  data = data)

# model 4: specification of taxa as a fixed effect (wiht the repeat measures random effect)
betareg_mod_5 <- glmmTMB(
  new_max_cover ~ trt * year_trt + New_taxon + (1|site_block),
  family = ordbeta(link = "logit"),
  data = data)

################################################################################

# check if the random intercepts vary 
VarCorr(betareg_mod_1)
VarCorr(betareg_mod_2)
VarCorr(betareg_mod_3)
VarCorr(betareg_mod_4) # here taxa is fixed - the repeat measures captures more of the variance
VarCorr(betareg_mod_5) # here taxa is fixed - site block captures very little of the variance

# all these values are moderate to high which suggests that species identity and the spatial
# grouping variable are larger drivers of cover values than time or treatment
# this is to be expected

################################################################################

summary(betareg_mod_4)
summary(betareg_mod_5)

################################################################################

create_model_summary <- function(model){
  d1 <- r2_nakagawa(model)
  d1a <- (d1[[1]])
  d1b <- (d1[[2]])
  d2 <-icc(model)
  d3 <- get_variance_residual(model)
  d4 <- get_variance_fixed(model)
  d5 <- get_variance_random(model)
  dat <- c(d1a, d1b, d2, d3, d4, d5)
  dat <- as.data.frame(dat)
  dat$var_total <- dat$var.residual + dat$var.fixed + dat$var.random
  dat$prop_fixed <- dat$var.fixed/dat$var_total *100
  dat$prop_random <- dat$var.random/dat$var_total *100
  dat$prop_residual <- dat$var.residual/dat$var_total*100
  return(dat)
}

mod_1 <- betareg_mod_1 %>%
  create_model_summary() %>%
  mutate(mod_nos = "mod_1") %>%
  mutate(random = "taxa") %>% 
  mutate(experiment = "DIST") %>%
  mutate(model = "Ordbeta")

mod_2 <- betareg_mod_2 %>%
  create_model_summary() %>%
  mutate(mod_nos = "mod_2") %>%
  mutate(random = "taxa") %>% 
  mutate(experiment = "DIST") %>%
  mutate(model = "Ordbeta")

mod_3 <- betareg_mod_3 %>%
  create_model_summary() %>%
  mutate(mod_nos = "mod_3") %>%
  mutate(random = "taxa") %>% 
  mutate(experiment = "DIST") %>%
  mutate(model = "Ordbeta")

mod_4 <- betareg_mod_4 %>%
  create_model_summary() %>%
  mutate(mod_nos = "mod_4") %>%
  mutate(random = "taxa") %>% 
  mutate(experiment = "DIST") %>%
  mutate(model = "Ordbeta")

mod_5 <- betareg_mod_5 %>%
  create_model_summary() %>%
  mutate(mod_nos = "mod_5") %>%
  mutate(random = "taxa") %>% 
  mutate(experiment = "DIST") %>%
  mutate(model = "Ordbeta")

mods <- bind_rows(mod_1, mod_2, mod_3, mod_4, mod_5)


################################################################################

# only for dist T0-T1

#function to extract model R2 and variances
create_model_summary <- function(model){
  
  r2_vals <- r2_nakagawa(model)
  dat <- tibble(
    r2_marginal    = unlist(r2_vals[1]),
    r2_conditional = unlist(r2_vals[2]),
    icc_adjusted   = unlist(icc(model)[[1]]),
    icc_conditional = unlist(icc(model)[[2]]),
    icc_unadjusted = unlist(icc(model)[[3]]),
    var.residual   = unlist(get_variance_residual(model)),
    var.fixed      = unlist(get_variance_fixed(model)),
    var.random     = unlist(get_variance_random(model))
  ) %>%
    mutate(
      var_total     = var.residual + var.fixed + var.random,
      prop_fixed    = 100 * var.fixed / var_total,
      prop_random   = 100 * var.random / var_total,
      prop_residual = 100 * var.residual / var_total
    )
  
  dat
}

# Define models and their labels
model_list <- list(
  mod_1 = betareg_mod_1,
  mod_2 = betareg_mod_2,
  mod_3 = betareg_mod_3,
  mod_4 = betareg_mod_4,
  mod_5 = betareg_mod_5
)

# Apply function to each model and bind results
mods_dist_T1 <- imap_dfr(model_list, ~ create_model_summary(.x) %>%
                   mutate(
                     mod_nos = .y,
                     experiment = "DIST",
                     model = "Ordbeta",
                     time_period = "T1"
                   ))


random_coeffs <- imap_dfr(model_list, ~ extract_random_coefs(.x))


# what is the difference between a random coefficient and a random effect here?
RCo_mod_1 <- extract_random_coefs(betareg_mod_1)
RCo_mod_1$type <- "coef"
RCo_mod_1$group
RCo_mod_2 <- extract_random_coefs(betareg_mod_2)
RCo_mod_3 <- extract_random_coefs(betareg_mod_3)

RE_mod_1 <- extract_random_effects(betareg_mod_1)
RE_mod_1$type <- "eff"
RE_mod_1$group
RE_mod_2 <- extract_random_effects(betareg_mod_2)
RE_mod_3 <- extract_random_effects(betareg_mod_3)

fixed_mod_4 <- extract_fixed_effects(betareg_mod_4)
fixed_mod_4 <- extract_fixed_effects(betareg_mod_4)

mod_1_effects <- bind_rows(RCo_mod_1, RE_mod_1)

# use extract_random_effects() when you want to know how individuals deviate from the whole model
# use extract_random_coeff() when you want to examine the variability in the relationships between
# predictors and the outcomes across different groups

################################################################################

# visualise these model results
# are there any differences between using the random coeffs and the random effects???
names(mod_1_effects)
mod_1_effects$group
mod_1_effects %>%
  ggplot(aes(value, reorder(group, value), colour = type))+
  theme_classic()+
  geom_point(size = 1)+
  geom_errorbar(aes(xmin = value - se, xmax = value + se), size = 0.3)+
  ylab("Taxon") +
  xlab("Standardised effects with SE")

# These values are both given on the model scale

# check to see if they are actually different from each other
RCo_mod_1 <- RCo_mod_1 %>% ungroup() %>% arrange(as.character(group))
RE_mod_1 <- RE_mod_1 %>% ungroup() %>% arrange(as.character(group))

# perfect positive correlation here - so no difference
cor.test(RCo_mod_1$value, RE_mod_1$value, method = "kendall")

################################################################################

# from now on lets stick to random coeffs to compare the diff models
nrow(RCo_mod_1) # 987
nrow(RCo_mod_2) # 1203 - more taxa???

n_distinct(RCo_mod_2)

# lets look at the difference between the random effetcs taken from the different models
RCo_mod_2 <- RCo_mod_2 %>% ungroup() %>% arrange(as.character(group))
cor.test(RCo_mod_1$value, RCo_mod_2$value, method = "kendall")





################################################################################

perf_mods <- compare_performance(betareg_mod_1, betareg_mod_2, betareg_mod_3, betareg_mod_4, betareg_mod_5)


see

# we need to extract these taxon level effects for each of the models and 

# visualise the variance decomposition between fixed, random and residual
plot_3  <-
  all %>%
  dplyr::filter(random == "both") %>%
  dplyr::select(prop_fixed:time_period) %>%
  pivot_longer(cols = prop_fixed:prop_residual, names_to = "part", values_to = "value") %>%
  ggplot(aes(x = experiment, y = value, fill= part, label = value)) +
  theme_classic()+
  theme(legend.title = element_blank())+
  geom_bar(position = "stack", stat = "identity") +
  facet_wrap(~time_period)+
  scale_fill_manual(values = c("black", "darkolivegreen2", "grey80"), 
                    labels = c("Fixed", "Random", "Residual"))+
  ylab("Variance Explained")+
  xlab("Experiment")+
  theme(legend.position = "bottom")+
  theme(text = element_text(size = 16))+
  theme(legend.text=element_text(size=rel(1)))

plot_3






library(marginaleffects)
avg_slopes(betareg_mod_3)




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
summary(betareg_mod)

re <- ranef(betareg_mod)[["New_taxon"]][["(Intercept)"]]
head(re)


library(modelsummary)
modelsummary(betareg_mod)
get_estimates(betareg_mod)

ranef(betareg_mod)$New_taxon[["(Intercept)"]]
re <- ranef(betareg_mod)[["New_taxon"]][["(Intercept)"]]
head(re)



# Fixed intercept
fixed_int <- fixef(betareg_mod)$cond["(Intercept)"]
fixed_int # -2.477423

# Random deviations
re <- ranef(betareg_mod)$New_taxon[, "(Intercept)"]
re

# random intercept deviations
str(ranef(betareg_mod))
betareg_mod$obj$
re <- ranef(betareg_mod)$New_taxon[["(Intercept)"]]
re <- ranef(betareg_mod)[["New_taxon"]][["(Intercept)"]]

coef(betareg_mod) # actual intercepts per 

# Combine
taxon_intercepts <- fixed_int + re
taxon_intercepts


# Predict for all levels (using a data frame with all levels)
# predictions_all_levels <- predict(mod, newdata = all_possible_newdata, re.form = NULL, type = "response")
                                                                       rdrr.io
                                                                                             How to calculate predicted values using glmmTMB?
                                                                                               13 Jul 2023
                                                                                             
                                                                                             Stack Overflow
                                                                                             
                                                                                             


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
