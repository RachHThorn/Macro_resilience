# R Thornley
# 18/08/2025
# Compare the results from the simple BACIs with the modelled GLMM outputs
# Project: P1_COMPADRE_DRAGNET
# Script: A3_compare_response_metrics

library(tidyverse)
library(ggpubr)

################################################################################

# read in the BACI and the RE form the GLMMs
RE_taxa <- read_csv("results/July_2025/RE_SE_Taxon_all_DRAGNet.csv")
# these effects are at the quadrat level
BACI_quadrat <- read_csv("results/July_2025/BACI_responses_simple_quadrat_all_DRAGNet.csv")
names(RE_taxa)
# group is the taxa variable
unique(RE_taxa$group) # there are 1086 unique species in this data set
names(BACI_quadrat)
# Taxon is the taxa variable
unique(BACI_quadrat$Taxon) # there are also 1086 unique species in this data set

# so we have effect sizes for 1086 species altogether
# these also need to be filtered for the DRAGNet species later on
# the list of species will vary depending on the time and the treatment we choose

# NOTE
# the values for BACI quadrat are used for the simple BACI calculations because 
# from them can be calculated a mean and SE around the mean 
# these can then be passed into the weighted regression models (meta-analytical model)
# this makes the model outputs more comparable with the GLMM work flow

################################################################################

# shorter code (output chat GPT)

BACI_all <- BACI_quadrat %>%
  group_by(trt, time, Taxon) %>%
  summarise(
    mean = mean(BACI, na.rm = TRUE),
    se = sd(BACI, na.rm = TRUE) / sqrt(n()),
    .groups = "drop"
  ) %>%
  drop_na()

# number of unique taxa per time point
BACI_all %>%
  count(trt, time) # total number of taxa in each group

################################################################################


# compare the two sets of outputs to see how different they are

# so to start with we have the problem that we have more estimates in the modelling 
# scenario than in the simple BACI scenario
# we need to find out how different they are

# setdiff finds elements in first vector not in second vector
setdiff(RE_dist_T1$group, BACI_dist_T1$Taxon) # there are  species here that have modelled effects for which we don't have simple BACI effects
setdiff(BACI_dist_T1$Taxon, RE_dist_T1$group) # there are no species that we have BACI effects for that are not captured in the GLMMs

# filter the RE list of taxa by the BACI list of taxa to compare the ordered factors
baci_taxa <- BACI$Taxon
baci_taxa # 61 species
RE <- RE %>% filter(group %in% baci_taxa)
RE$group # 52 species
setdiff(RE$group, BACI$Taxon) # now they are the same length and can be compared

str(RE)
str(BACI)
# the taxa are both still character vectors here and we need to change them to
# factors and order them by the value of the effect

################################################################################

# one approach is to take ranks and then work out if these ranks are different

# create an ordered factor with ranks based on the filtered dfs
# for both datasets - then join them
new_RE <- RE %>% ungroup %>% arrange(group) %>% 
  mutate(RE_rank = factor(rank(value), ordered = TRUE)) %>% 
  select(group, RE_rank) %>%
  rename(Taxon = group)

new_BACI <- BACI %>% ungroup() %>% arrange(Taxon) %>% 
  mutate(BACI_rank = factor(rank(mean_dist), ordered = TRUE)) %>%
  select(Taxon, BACI_rank)

setdiff(new_RE$Taxon, new_BACI$Taxon) # the Taxon names are still the same

ranks <- new_RE %>% left_join(new_BACI)

###############################################################################

# perform chi-squared test
table_data <- table(ranks$RE_rank, ranks$BACI_rank)
table_data
# Chi-square test
chi_test <- chisq.test(table_data)
print(chi_test)
# p - value = 0.25 - so the variables are not strongly related

###############################################################################

# perform wilcox test
BACI_rank <- as.numeric(ranks$BACI_rank)
BACI_rank
RE_rank <- as.numeric(ranks$RE_rank)
RE_rank
wilcox <- wilcox.test(BACI_rank, RE_rank, paired = TRUE)
# this is in conflict with the results from the chi-squared test
wilcox
# p value = 0.6 / no significant shifts in rankings?
plot(BACI_rank, RE_rank)

# Pearson's product moment correlation
cor.test(BACI_rank, RE_rank, method = "spearman") # significantly diff

# Cohen's Kappa
library(irr)
kappa2(data.frame(RE_rank, BACI_rank)) # with numeric data
kappa2(data.frame(ranks$RE_rank, ranks$BACI_rank)) # with ranked ordered data
# very small k value - total disagreement between ranks
# p = 0.31

################################################################################

# lets look at this for ALL the species in DRAGNet


################################################################################

# Get effects for all species across the whole of DRAGnet 
# then compare them to see how drastic these differences are 

# read in the BACI and the RE form the GLMMs
BACI_quadrat <- read_csv("results/March_2025/BACI_responses_simple_quadrat_all_DRAGNet.csv")
RE_taxa <- read_csv("results/March_2025/RE_SE_Taxon_all_DRAGNet.csv")

# get the taxon level mean BACIs for one time period / trt
BACI <-
  BACI_quadrat %>%  
  select(Taxon, site_name, group_var, time, BACI_dist) %>% 
  filter(time == "T0-T1") %>%
  group_by(time, Taxon) %>%
  mutate(mean_dist = mean(BACI_dist), se = sd(BACI_dist) / sqrt(n())) %>%
  distinct(mean_dist, .keep_all= TRUE) %>%
  drop_na() %>%
  ungroup()
n_distinct(BACI$Taxon) # 592 taxon that we were able to get effects for

RE <-
  RE_taxa %>%
  filter(model == "Ordbeta") %>%
  filter(experiment == "DIST") %>%
  filter(time_period == "T0-T1") %>%
  ungroup() %>%
  rename(Taxon = group) 
n_distinct(RE$Taxon) # 716 species that we were able to get effects for

################################################################################

# compare the two sets of outputs to see how different they are

# so to start with we have the problem that we have more estimates in the modelling 
# scenario than in the simple BACI scenario
# we need to find out how different they are

# setdiff finds elements in first vector not in second vector
setdiff(RE$Taxon, BACI$Taxon) # there are 147 species here that have modelled effects for which we don't have simple BACI effects
setdiff(BACI$Taxon, RE$Taxon) # there are 23 species here that have BACI effects for which we don't have modelled effects

# filter the RE list of taxa by the BACI list of taxa to compare the ordered factors
baci_taxa <- BACI$Taxon
baci_taxa # 52 species
RE <- RE %>% filter(Taxon %in% baci_taxa)
RE$Taxon # 567 species
setdiff(RE$Taxon, BACI$Taxon) # now they are the same length and can be compared

# one approach is to take ranks and then work out if these ranks are different

# create an ordered factor with ranks based on the filtered dfs
# for both datasets - then join them
new_RE <- RE %>% ungroup %>% arrange(Taxon) %>% 
  mutate(RE_rank = factor(rank(value), ordered = TRUE)) %>% 
  select(Taxon, RE_rank) 

new_BACI <- BACI %>% ungroup() %>% arrange(Taxon) %>% 
  mutate(BACI_rank = factor(rank(mean_dist), ordered = TRUE)) %>%
  select(Taxon, BACI_rank)

setdiff(new_RE$Taxon, new_BACI$Taxon) # the Taxon names are still the same

ranks <- new_RE %>% left_join(new_BACI)

###############################################################################

# perform chi-squared test
table_data <- table(ranks$RE_rank, ranks$BACI_rank)
table_data
# Chi-square test
chi_test <- chisq.test(table_data)
print(chi_test)
# p value is very low: there is a likely association between the values
# but this is definitely not linear

# perform wilcox test
BACI_rank <- as.numeric(ranks$BACI_rank)
BACI_rank
RE_rank <- as.numeric(ranks$RE_rank)
RE_rank
wilcox <- wilcox.test(BACI_rank, RE_rank, paired = TRUE)
# this is in conflict with the results from the chi-squared test
wilcox
# p value is very small - 
plot(BACI_rank, RE_rank)

# Pearson's product moment correlation
cor.test(BACI_rank, RE_rank, method = "spearman") # significantly diff
# again a very small P value

################################################################################
