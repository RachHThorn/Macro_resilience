# R Thornley
# 27/08/2025
# Compare the results from the simple BACIs with the modelled GLMM outputs
# Project: P1_COMPADRE_DRAGNET
# Script: A3_compare_response_metrics

library(tidyverse)
library(ggpubr)

################################################################################

# read in the BACI and the RE form the GLMMs
RE_taxa <- read_csv("results/RE_SE_Taxon_all_DRAGNet.csv")
# these effects are at the quadrat level
BACI_quadrat <- read_csv("results/BACI_responses_simple_quadrat_all_time_periods.csv")
names(RE_taxa)
# group is the taxa variable
n_distinct(RE_taxa$group) # there are 1224 unique species in this data set
# Taxon is the taxa variable
n_distinct(BACI_quadrat$Taxon) # there are also 1224 unique species in this data set

# so we have effect sizes for 1224 species altogether
# these also need to be filtered for the DRAGNet species later on
# the list of species will vary depending on the time and the treatment we choose

# NOTE
# the values for BACI quadrat are used for the simple BACI calculations because 
# from them can be calculated a mean and SE around the mean 
# these can then be passed into the weighted regression models (meta-analytical model)
# this makes the model outputs more comparable with the GLMM work flow

################################################################################

# calculate the mean per taxa and the associated SE
BACI_all <- BACI_quadrat %>%
  group_by(trt, time, Taxon) %>%
  summarise(
    mean = mean(BACI, na.rm = TRUE),
    se = sd(BACI, na.rm = TRUE) / sqrt(n()),
    .groups = "drop"
  ) 

# if we drop_na - it takes out all the effects for which we have only one data point
BACI_all %>% drop_na() %>% pull(Taxon) %>% n_distinct() # if we do this there are only 873 taxa 


# find out how many taxa are in each BACI experiment and at each of the time periods

# unique taxa per time point and treatment
BACI_all %>%
  count(trt, time) # total number of taxa in each group

# number of unique taxa per time point and treatment
names(RE_taxa)
RE_taxa %>%
  filter(model == "Ordbeta") %>%
  count(experiment, time_period, model)

################################################################################

# compare the two sets of outputs to see how different they are

# so to start with we have the problem that we have more estimates in the modelling 
# scenario than in the simple BACI scenario
# we need to find out how different they are

# break the dat down into the diff. experiments and times

DIST_T1_BACI <- BACI_all %>% filter(trt == "DIST") %>% filter(time == "T0-T1") %>% arrange(Taxon)
DIST_T1_RE <- RE_taxa %>% 
  filter(model == "Ordbeta", experiment == "DIST", time_period == 1) %>% 
  rename(Taxon = group) %>%
  arrange(Taxon)

################################################################################

# Kendall's rank correlation
cor.test(DIST_T1_BACI$mean, DIST_T1_RE$value, method = "kendall")
plot(DIST_T1_BACI$mean, DIST_T1_RE$value)

# there is a significant negative correlation between the two data sets but it is very weak
# and probably not that important

# just to check this has worked how we wanted
# Join datasets by species
combined <- inner_join(DIST_T1_BACI, DIST_T1_RE, by = "Taxon")
# Compute Kendall's tau on the paired measures
cor.test(combined$mean, combined$value, method = "kendall")

################################################################################

# What about if we transform the data 

