# R Thornley
# 27/08/2025
# Project: P1_COMPADRE_DRAGNET
# Script: S4_simple_BACI_calculations

# calculate simple BACI equation for each of the experiments 
# for both the time points

library(tidyverse)

################################################################################

# wrapper function that creates the BACI equation for a specified time-period and treatment
compute_BACI <- function(data, control, treatment, before, after, time_label, trt_label) {
  
  wide <- data %>%
    filter(trt %in% c(control, treatment)) %>%
    group_by(site_name, New_taxon, year_trt, trt, block) %>%
    summarise(mean_cover = mean(new_max_cover), .groups = "drop") %>%
    mutate(trt_time = paste(trt, year_trt, sep = "_")) %>%
    select(site_name, New_taxon, block, trt_time, mean_cover) %>%
    pivot_wider(values_from = mean_cover, names_from = trt_time)
  
  # the four columns we need for BACI
  need_cols <- c(
    paste0(control,    "_", before),
    paste0(treatment,  "_", before),
    paste0(control,    "_", after),
    paste0(treatment,  "_", after)
  )
  
  # ensure they exist (fill missing with 0)
  for (nm in need_cols) {
    if (!nm %in% names(wide)) wide[[nm]] <- 0
  }
  
  wide %>%
    mutate(across(all_of(need_cols), ~ replace_na(.x, 0))) %>%
    # drop rows where all four inputs are 0
    filter(if_any(all_of(need_cols), ~ . != 0)) %>%
    transmute(
      Taxon = New_taxon,
      site_name,
      block,
      BACI =
        (.data[[paste0(treatment, "_", after)]] - .data[[paste0(control, "_", after)]]) -
        (.data[[paste0(treatment, "_", before)]] - .data[[paste0(control, "_", before)]]),
      group_var = "quadrat",
      time = time_label,
      trt  = trt_label
    )
}

################################################################################


# -------------------------------
# T0–T1
# -------------------------------

T1 <- read_csv("results/DRAGNet_T0_T1_all.csv")
 
Quadrat_T1_DIST  <- compute_BACI(data = T1, control = "Control", treatment = "Disturbance", 
                                 before = "T0", after = "T1", time_label = "T0-T1", trt_label = "DIST")
Quadrat_T1_NPK   <- compute_BACI(data = T1, control = "Control", treatment = "NPK", 
                                 before = "T0", after = "T1", time_label = "T0-T1", trt_label = "NPK")
Quadrat_T1_INTER <- compute_BACI(data = T1, control = "Control", treatment = "NPK+Disturbance", 
                                 before = "T0", after = "T1", time_label = "T0-T1", trt_label = "INTER")

T1_result <- bind_rows(Quadrat_T1_DIST, Quadrat_T1_NPK, Quadrat_T1_INTER)
write_csv(T1_result, "results/BACI_T0_T1_results.csv")

# -------------------------------
# T0–T2
# -------------------------------

T2 <- read_csv("results/DRAGNet_T0_T2_all.csv")

Quadrat_T2_DIST  <- compute_BACI(data = T2, control = "Control", treatment = "Disturbance", 
                                 before = "T0", after = "T2", time_label = "T0-T2", trt_label = "DIST")
Quadrat_T2_NPK   <- compute_BACI(data = T2, control = "Control", treatment = "NPK", 
                                 before = "T0", after = "T2", time_label = "T0-T2", trt_label = "NPK")
Quadrat_T2_INTER <- compute_BACI(data = T2, control = "Control", treatment = "NPK+Disturbance", 
                                 before = "T0", after = "T2", time_label = "T0-T2", trt_label = "INTER")

T2_result <- bind_rows(Quadrat_T2_DIST, Quadrat_T2_NPK, Quadrat_T2_INTER)
write_csv(T2_result, "results/BACI_T0_T2_results.csv")

# -------------------------------
# T0–T3
# -------------------------------

T3 <- read_csv("results/DRAGNet_T0_T3_all.csv")

Quadrat_T3_DIST  <- compute_BACI(data = T3, control = "Control", treatment = "Disturbance", 
                                 before = "T0", after = "T3", time_label = "T0-T3", trt_label = "DIST")
Quadrat_T3_NPK   <- compute_BACI(data = T3, control = "Control", treatment = "NPK", 
                                 before = "T0", after = "T3", time_label = "T0-T3", trt_label = "NPK")
Quadrat_T3_INTER <- compute_BACI(data = T3, control = "Control", treatment = "NPK+Disturbance", 
                                 before = "T0", after = "T3", time_label = "T0-T3", trt_label = "INTER")

T3_result <- bind_rows(Quadrat_T3_DIST, Quadrat_T3_NPK, Quadrat_T3_INTER)
write_csv(T3_result, "results/BACI_T0_T3_results.csv")

# -------------------------------
# Combine all
# -------------------------------
all <- bind_rows(T1_result, T2_result, T3_result)
write_csv(all, "results/BACI_responses_simple_quadrat_all_time_periods.csv")
