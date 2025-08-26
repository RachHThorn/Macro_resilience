# R Thornley
# 10/02/2025
# Project: P1_COMPADRE_DRAGNET
# Script: S3_GLMM_extract_random_effects_all_DRAG_data

# find the effect sizes for all of the DRAGNet species to contextualise the overlap

# 1) For T0-T1 find effects for each species
# 2) For T0-T1 find effects for each species within site
# 3) For T0-T2 find effects for each species
# 4) For T0-T2 find effects for each species within site
# 5) For T0-T3 find effects for each species
# 6) For T0-T3 find effects for each species within site

# for each section - effects are found per experiment
# Disturbance
# NPK
# NPK + Disturbance

library(tidyverse)
# use the glmmTMB package for frequentist hurdle and ordered beta reg models
library(glmmTMB)
# use the mixed up package for extraction of the random effects
library(mixedup)

# these functions run three model types on the DRAGNet data and extract effects with Standard Errors
# 1) Beta zero-inflated hurdle model; 2) Frequentist ordered beta model;
# We extract taxon random effects in the first function and taxon within site random effects in the second function

get_RE_taxon <- function(data, experiments) {
  data <- data %>% filter(trt %in% experiments)
  hurdle_mod <- glmmTMB(new_max_cover ~ trt * year_trt + (1|New_taxon),
                        family = beta_family(link = "logit"), ziformula = ~ 1 + trt * year_trt,
                        data = data)
  betareg_mod <- glmmTMB(new_max_cover ~ trt * year_trt + 
                           (1|New_taxon),
                         family = ordbeta(link = "logit"), 
                         data = data)
  Hurdle_cond <- hurdle_mod %>% extract_random_coefs(component = "cond") %>%
    mutate(model = "Hurdle_cond")
  Hurdle_zi <- hurdle_mod %>% extract_random_coefs(component = "zi") %>%
    mutate(model = "Hurdle_zi")
  Ordbeta <- betareg_mod %>% extract_random_coefs() %>%
    mutate(model = "Ordbeta")
  results <- rbind(Hurdle_cond, Hurdle_zi, Ordbeta)
  return(results)
}

get_RE_taxon_site <- function(data, experiments) {
  data <- data %>% filter(trt %in% experiments)
  hurdle_mod <- glmmTMB(new_max_cover ~ trt * year_trt + (1|taxon_site),
                        family = beta_family(link = "logit"), ziformula = ~ 1 + trt * year_trt,
                        data = data)
  betareg_mod <- glmmTMB(new_max_cover ~ trt * year_trt + 
                           (1|taxon_site),
                         family = ordbeta(link = "logit"), 
                         data = data)
  Hurdle_cond <- hurdle_mod %>% extract_random_coefs(component = "cond") %>%
    mutate(model = "Hurdle_cond")
  Hurdle_zi <- hurdle_mod %>% extract_random_coefs(component = "zi") %>%
    mutate(model = "Hurdle_zi")
  Ordbeta <- betareg_mod %>% extract_random_coefs() %>%
    mutate(model = "Ordbeta")
  results <- rbind(Hurdle_cond, Hurdle_zi, Ordbeta)
  return(results)
}

################################################################################

# Load the data for modelling

dat_1 <- read_csv("results/July_2025/DRAGNet_all_clean_data_with_zeros_T0_T1.csv")
# check the data has no NAs first
dat_1 %>%
  filter(if_any(everything(), is.na))

# run the get_RE_taxon on each of the three experiments
experiments <- c("Control", "Disturbance")
DIST <- get_RE_taxon(dat_1, experiments)
DIST$time_period <- "T0-T1"
DIST$experiment <- "DIST"

experiments <- c("Control", "NPK")
NPK <- get_RE_taxon(dat_1, experiments)
NPK$time_period <- "T0-T1"
NPK$experiment <- "NPK"

experiments <- c("Control", "NPK+Disturbance")
INTER <- get_RE_taxon(dat_1, experiments)
INTER$time_period <- "T0-T1"
INTER$experiment <- "INTER"

# join up the results for the first time period
T1_taxon <- rbind(DIST, NPK, INTER)

# run the get_RE_taxon_site on each of the three experiments
experiments <- c("Control", "Disturbance")
DIST <- get_RE_taxon_site(dat_1, experiments)
DIST$time_period <- "T0-T1"
DIST$experiment <- "DIST"

experiments <- c("Control", "NPK")
NPK <- get_RE_taxon_site(dat_1, experiments)
NPK$time_period <- "T0-T1"
NPK$experiment <- "NPK"

experiments <- c("Control", "NPK+Disturbance")
INTER <- get_RE_taxon_site(dat_1, experiments)
INTER$time_period <- "T0-T1"
INTER$experiment <- "INTER"

T1_taxon_site <- rbind(DIST, NPK, INTER)

################################################################################

# Time period 2 #

# read in data and clean up
dat_2 <- read_csv("results/July_2025/DRAGNet_all_clean_data_with_zeros_T0_T2.csv")
dat_2 %>%
  filter(if_any(everything(), is.na))

range(dat_2$new_max_cover)

experiments <- c("Control", "Disturbance")
DIST <- get_RE_taxon(dat_2, experiments)
DIST$time_period <- "T0-T2"
DIST$experiment <- "DIST"

experiments <- c("Control", "NPK")
NPK <- get_RE_taxon(dat_2, experiments)
NPK$time_period <- "T0-T2"
NPK$experiment <- "NPK"

experiments <- c("Control", "NPK+Disturbance")
INTER <- get_RE_taxon(dat_2, experiments)
INTER$time_period <- "T0-T2"
INTER$experiment <- "INTER"

T2_taxon <- rbind(DIST, NPK, INTER)

# run the get_RE_taxon_site on each of the three experiments
experiments <- c("Control", "Disturbance")
DIST <- get_RE_taxon_site(dat_2, experiments)
DIST$time_period <- "T0-T2"
DIST$experiment <- "DIST"

experiments <- c("Control", "NPK")
NPK <- get_RE_taxon_site(dat_2, experiments)
NPK$time_period <- "T0-T2"
NPK$experiment <- "NPK"

experiments <- c("Control", "NPK+Disturbance")
INTER <- get_RE_taxon_site(dat_2, experiments)
INTER$time_period <- "T0-T2"
INTER$experiment <- "INTER"

T2_taxon_site <- rbind(DIST, NPK, INTER)

################################################################################

# Time period 3 #

# read in data and clean up
dat_3 <- read_csv("results/July_2025/DRAGNet_all_clean_data_with_zeros_T0_T3.csv")
dat_3 %>%
  filter(if_any(everything(), is.na))

range(dat_3$new_max_cover)

experiments <- c("Control", "Disturbance")
DIST <- get_RE_taxon(dat_3, experiments)
DIST$time_period <- "T0-T3"
DIST$experiment <- "DIST"

experiments <- c("Control", "NPK")
NPK <- get_RE_taxon(dat_3, experiments)
NPK$time_period <- "T0-T3"
NPK$experiment <- "NPK"

experiments <- c("Control", "NPK+Disturbance")
INTER <- get_RE_taxon(dat_3, experiments)
INTER$time_period <- "T0-T3"
INTER$experiment <- "INTER"

T3_taxon <- rbind(DIST, NPK, INTER)

# run the get_RE_taxon_site on each of the three experiments
experiments <- c("Control", "Disturbance")
DIST <- get_RE_taxon_site(dat_3, experiments)
DIST$time_period <- "T0-T3"
DIST$experiment <- "DIST"

experiments <- c("Control", "NPK")
NPK <- get_RE_taxon_site(dat_3, experiments)
NPK$time_period <- "T0-T3"
NPK$experiment <- "NPK"

experiments <- c("Control", "NPK+Disturbance")
INTER <- get_RE_taxon_site(dat_3, experiments)
INTER$time_period <- "T0-T3"
INTER$experiment <- "INTER"

T3_taxon_site <- rbind(DIST, NPK, INTER)


# join and export data for modelling
# NOTE: these effects need filtering for the COMPADRE taxa overlap

all_taxon <- rbind(T1_taxon, T2_taxon, T3_taxon)
write_csv(all_taxon, "results/July_2025/RE_SE_Taxon_all_DRAGNet.csv")

all_taxon_site <- rbind(T1_taxon_site, T2_taxon_site, T3_taxon_site)
write_csv(all_taxon_site, "results/July_2025/RE_SE_Taxon_site_all_DRAGNet.csv")

################################################################################