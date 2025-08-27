# R Thornley
# 10/02/2025
# Project: P1_COMPADRE_DRAGNET
# Script: S3_GLMM_extract_random_effects_all_DRAG_data (tidier code)

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
  
  hurdle_mod <- glmmTMB(
    new_max_cover ~ trt * year_trt + (1|New_taxon),
    family = beta_family(link = "logit"),
    ziformula = ~ 1 + trt * year_trt,
    data = data
  )
  
  betareg_mod <- glmmTMB(
    new_max_cover ~ trt * year_trt + (1|New_taxon),
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
  
  hurdle_mod <- glmmTMB(
    new_max_cover ~ trt * year_trt + (1|taxon_site),
    family = beta_family(link = "logit"),
    ziformula = ~ 1 + trt * year_trt,
    data = data
  )
  
  betareg_mod <- glmmTMB(
    new_max_cover ~ trt * year_trt + (1|taxon_site),
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
# Wrapper to process each time period
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
  "results/DRAGNet_T0_T1_overlap.csv",
  "results/DRAGNet_T0_T2_overlap.csv",
  "results/DRAGNet_T0_T3_overlap.csv")

results <- imap(files, ~ process_time_period(.x, .y))

all_taxon <- map(results, "taxon") %>% bind_rows()
all_taxon_site <- map(results, "taxon_site") %>% bind_rows()

#-----------------------------------------------------------
# Export
#-----------------------------------------------------------

write_csv(all_taxon,      "results/RE_SE_Taxon_all_DRAGNet.csv")
write_csv(all_taxon_site, "results/RE_SE_Taxon_site_all_DRAGNet.csv")
