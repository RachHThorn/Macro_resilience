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

# model 2: addition of site_block as the RE alongside the taxon
betareg_mod_2 <- glmmTMB(
  new_max_cover ~ trt * year_trt + (1|New_taxon) + (1|site_block),
  family = ordbeta(link = "logit"),
  data = data)

# model 3: addition of site_block_taxon as the RE (I think this is the right specification for pairing - repeat measures)
betareg_mod_3 <- glmmTMB(
  new_max_cover ~ trt * year_trt + (1|New_taxon) + (1|site_block_taxon),
  family = ordbeta(link = "logit"),
  data = data)

# model 4: specification of taxa as a fixed effect (with the repeat measures random effect)
betareg_mod_4 <- glmmTMB(
  new_max_cover ~ trt * year_trt + New_taxon + (1|site_block_taxon),
  family = ordbeta(link = "logit"),
  data = data)

# model 5: specification of taxa as a fixed effect (with the site_block random effect)
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

mod_4_summary <- summary(betareg_mod_4)
mod_5_summary <- summary(betareg_mod_5)

################################################################################

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
# this is only for the DIST and T1 data - needs re-writing for the other data sets too
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


# from models 4 and 5 we wnaat to extract the fixed effects only
fixed_mod_4 <- extract_fixed_effects(betareg_mod_4)
fixed_mod_4 <- 
  fixed_mod_4 %>% 
  filter(str_detect(term, "New_taxon")) %>%
  mutate(new_term = str_remove(term, "New_taxon"))

fixed_mod_5 <- extract_fixed_effects(betareg_mod_5)
fixed_mod_5 <- 
  fixed_mod_5 %>% 
  filter(str_detect(term, "New_taxon")) %>%
  mutate(new_term = str_remove(term, "New_taxon"))

# use extract_random_effects() when you want to know how individuals deviate from the whole model
# use extract_random_coeff() when you want to examine the variability in the relationships between
# predictors and the outcomes across different groups

################################################################################

# visualise these model results
# are there any differences between using the random coeffs and the random effects???
# we can bind the two types of effects and see if they are any different from each other
mod_1_effects <- bind_rows(RCo_mod_1, RE_mod_1)
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
nrow(RCo_mod_2) # 1203 - more rows than we need - 
# the effects for the other grouping variable are also included in this table

# filter the larger data frame so both only have the taxa effects listed
names(RCo_mod_2)
RCo_mod_2 <- RCo_mod_2 %>% filter(group_var == "New_taxon")

# lets look at the difference between the random effects taken from the different models
RCo_mod_2 <- RCo_mod_2 %>% ungroup() %>% arrange(as.character(group))
cor.test(RCo_mod_1$value, RCo_mod_2_small$value, method = "kendall")
# a very strong positive correlation - there is not much difference between the effects

################################################################################

nrow(RCo_mod_3) 
RCo_mod_3 <- RCo_mod_3 %>% filter(group_var == "New_taxon")
RCo_mod_3 <- RCo_mod_3 %>% ungroup() %>% arrange(as.character(group))
cor.test(RCo_mod_2$value, RCo_mod_3$value, method = "kendall")
cor.test(RCo_mod_1$value, RCo_mod_3$value, method = "kendall")

plot(RCo_mod_2$value, RCo_mod_3$value)
plot(RCo_mod_1$value, RCo_mod_2$value)
plot(RCo_mod_1$value, RCo_mod_3$value)

################################################################################

# lets compare the random effects from these models to the BACI simple response
simple_BACI <- read_csv("results/BACI_T0_T1_results.csv") %>% 
  filter(trt == "DIST", time == "T0-T1")
simple_BACI

BACI_all <- simple_BACI %>%
  group_by(trt, time, Taxon) %>%
  summarise(
    mean = mean(BACI, na.rm = TRUE),
    se = sd(BACI, na.rm = TRUE) / sqrt(n()),
    .groups = "drop"
  ) 

n_distinct(BACI_all$Taxon) # 987 taxa that we have estimates for
BACI_all <- BACI_all %>% ungroup() %>% arrange(as.character(Taxon))
names(BACI_all)

# compare mods 1-3 wiht BACI response
cor.test(RCo_mod_1$value, BACI_all$mean, method = "kendall")
plot(RCo_mod_1$value, BACI_all$mean)

cor.test(RCo_mod_2$value, BACI_all$mean, method = "kendall")
plot(RCo_mod_2$value, BACI_all$mean)

cor.test(RCo_mod_3$value, BACI_all$mean, method = "kendall")
plot(RCo_mod_3$value, BACI_all$mean)

################################################################################

# lets look at the other models where the fixed effect is taxon

fixed_mod_4 <- 
  fixed_mod_4 %>% 
  ungroup() %>% 
  arrange(as.character(new_term))

fixed_mod_5 <- 
  fixed_mod_5 %>% 
  ungroup() %>% 
  arrange(as.character(new_term))

# so the correlation between the two fixed effects modes are very similar
cor.test(fixed_mod_4$value, fixed_mod_5$value, method = "kendall")
plot(fixed_mod_4$value, fixed_mod_5$value)


# lets look at the correlation between effects from the RE models (1-3) and the fixed effects (4-5)
n_distinct(RCo_mod_1$group) # 987
n_distinct(fixed_mod_4$new_term) # 986

# difference between the vectors
setdiff(RCo_mod_1$group, fixed_mod_4$new_term)
# fixed mod 4 contains the taxon "Abutilon_fruticosum"
RCo_mod_1 <- RCo_mod_1 %>% filter(!group == "Abutilon_fruticosum")

# now test for the correlation
cor.test(fixed_mod_4$value, RCo_mod_1$value, method = "kendall")
plot(fixed_mod_4$value, RCo_mod_1$value)

# how about the correlation between the BACI and the fixed model
n_distinct(BACI_all$Taxon) # 987
n_distinct(fixed_mod_4$new_term) # 986

# difference between the vectors
setdiff(BACI_all$Taxon, fixed_mod_4$new_term)
# fixed mod 4 contains the taxon "Abutilon_fruticosum"
BACI_all <- BACI_all %>% filter(!Taxon == "Abutilon_fruticosum")
cor.test(fixed_mod_4$value, BACI_all$mean, method = "kendall")
plot(fixed_mod_4$value, BACI_all$mean)

################################################################################

# see if the interaction term is significant in the fixed effetc model
fixed_mod_4 <- extract_fixed_effects(betareg_mod_4)
fixed_results <- fixed_mod_4 %>% filter(!str_detect(term, "New_taxon"))
fixed_results

fixed_mod_5 <- extract_fixed_effects(betareg_mod_5)
fixed_results <- fixed_mod_5 %>% filter(!str_detect(term, "New_taxon"))
fixed_results

################################################################################

# compare the performance of the 5 models
perf_mods <- compare_performance(betareg_mod_1, betareg_mod_2, betareg_mod_3, betareg_mod_4, betareg_mod_5)

################################################################################

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



################################################################################

# try to predict from the models onto the scale of the response
# this is all very rough and doesn't work....


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
