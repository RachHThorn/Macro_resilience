# R Thornley
# 29/08/2025
# Project: P1_COMPADRE_DRAGNET
# Script: S3_GLMM_extract_random_effects_all_DRAG_data

#-----------------------------------------------------------
# Instructions
#-----------------------------------------------------------

# Read in the tidy DRAGNet data for each time period we are interested in
# ie T0-T1, T0-T2, T0-T3
# Function Get_RE_taxon - runs two types of GLMMs for the each of the treatments
# in order to get the REs for only the taxon
# Function Get_RE_taxon_site - runs two types of GLMMs for the each of the treatments
# in order to get the REs for the taxon within site
# wrapper function enables these models to be run for each of the different experiments 
# code is then run for each separate data set for each time period
# then joined and exported

#-----------------------------------------------------------
# Packages
#-----------------------------------------------------------

# install the mixed up package directly from Github
remotes::install_github('m-clark/mixedup')

# load packages
packages <- c("tidyverse", "glmmTMB", "rstan", "assertthat", "purrr", "remotes", "mixedup")

for (p in packages) {
  if (!requireNamespace(p, quietly = TRUE)) {
    install.packages(p)
  }
  library(p, character.only = TRUE)
}

#-----------------------------------------------------------
# Functions
#-----------------------------------------------------------

get_RE_taxon <- function(data, experiments) {
  
  data <- data %>% filter(trt %in% experiments)
  data <- data %>% mutate(site_block_taxon = paste0(site_name, "_", block, "_", New_taxon))
  
  hurdle_mod <- glmmTMB(
    new_max_cover ~ trt * year_trt + (1|New_taxon) + (1|site_block_taxon),
    family = beta_family(link = "logit"),
    ziformula = ~ 1 + trt * year_trt,
    data = data
  )
  
  betareg_mod <- glmmTMB(
    new_max_cover ~ trt * year_trt + (1|New_taxon) + (1|site_block_taxon),
    family = ordbeta(link = "logit"),
    data = data
  )
  
  bind_rows(
    extract_random_coefs(hurdle_mod, component = "cond") %>% mutate(model = "Hurdle_cond"),
    extract_random_coefs(hurdle_mod, component = "zi")   %>% mutate(model = "Hurdle_zi"),
    extract_random_coefs(betareg_mod)                    %>% mutate(model = "Ordbeta")
  )
}

get_RE_taxon_site <- function(data, experiments) {
  
  data <- data %>% filter(trt %in% experiments)
  data <- data %>% mutate(site_block_taxon = paste0(site_name, "_", block, "_", New_taxon))
  
  hurdle_mod <- glmmTMB(
    new_max_cover ~ trt * year_trt + (1|taxon_site) + (1|site_block_taxon),
    family = beta_family(link = "logit"),
    ziformula = ~ 1 + trt * year_trt,
    data = data
  )
  
  betareg_mod <- glmmTMB(
    new_max_cover ~ trt * year_trt + (1|taxon_site) + (1|site_block_taxon),
    family = ordbeta(link = "logit"),
    data = data
  )
  
  bind_rows(
    extract_random_coefs(hurdle_mod, component = "cond") %>% mutate(model = "Hurdle_cond"),
    extract_random_coefs(hurdle_mod, component = "zi")   %>% mutate(model = "Hurdle_zi"),
    extract_random_coefs(betareg_mod)                    %>% mutate(model = "Ordbeta")
  )
}

#-----------------------------------------------------------
# Wrapper to process each experiment (treatment)
#-----------------------------------------------------------

process_time_period <- function(file, time_label) {
  dat <- read_csv(file)
  
  exp_list <- list(
    DIST  = c("Control", "Disturbance"),
    NPK   = c("Control", "NPK"),
    INTER = c("Control", "NPK+Disturbance")
  )
  
  taxon <- map_dfr(names(exp_list), function(e) {
    get_RE_taxon(dat, exp_list[[e]]) %>%
      mutate(time_period = time_label, experiment = e)
  })
  
  taxon_site <- map_dfr(names(exp_list), function(e) {
    get_RE_taxon_site(dat, exp_list[[e]]) %>%
      mutate(time_period = time_label, experiment = e)
  })
  
  list(taxon = taxon, taxon_site = taxon_site)
}

#-----------------------------------------------------------
# Run for all time periods
#-----------------------------------------------------------

files <- c(
  "results/DRAGNet_T0_T1_all.csv",
  "results/DRAGNet_T0_T2_all.csv",
  "results/DRAGNet_T0_T3_all.csv")

results <- imap(files, ~ process_time_period(.x, .y))

all_taxon <- map(results, "taxon") %>% bind_rows()
all_taxon_site <- map(results, "taxon_site") %>% bind_rows()

#-----------------------------------------------------------
# Export
#-----------------------------------------------------------

write_csv(all_taxon,      "results/RE_SE_Taxon_all_DRAGNet.csv")
write_csv(all_taxon_site, "results/RE_SE_Taxon_site_all_DRAGNet.csv")
