# R Thornley
# Macro-ecology DRAGNet and COMPADRE paper
# 2024-2025

################################################################################

# Required Packages #

library(tidyverse) # general data wrangling
library(Rcompadre) # package to interact with COMPADRE database
library(glmmTMB) # GLMMs package - runs the ordered beta-reg and hurdle models
library(lme4) # provides functions to interrogate GLMMs
library(mixedup) #  provides functions to interrogate GLMMs
library(popbio) # basic MPM analyses
library(popdemo) # advanced MPM transient measures
library(Rage) # life history metrics from MPMs
library(DHARMa) # residual simulation and model checking for GLMMs

################################################################################

# Data Processing Scripts #

# script 1: Get a list of species that overlap in DRAGNet and COMPADRE etc.
source("code/S1_summary_information_DRAGNet_COMPADRE.R")

# script 2: Get_master_data_DRAGNet_modelling
source("code/S2_get_master_data_DRAGNet_modelling.R")

# script 3: Run GLMMs and extract Random Effect estimates at species and species within site level
source("code/S3_GLMM_extract_random_effects_all_DRAGNet_cover_data_RE.R")

# script 4: Run the simple BACI calculation across the DRAGNet data base
source("code/S4_simple_BACI_calculations.R")

# script 5 : Get demographic variables from MPMs found in COMPADRE 
source("code/S5_get_demographic_metrics.R")

# script 6 : Get community context variables from DRAGNet plots
source("code/S6_get_community_metrics.R")

# script 7 : Get community context variables from functional databases
source("code/S7_species_strategies_calculate.R")

# script 8 : Test for phylogenetic signal in demographic traits
source("code/S8_test_for_phylogenetic_signal.R")

# script 9 : Model DRAGNet Taxon effects from GLMMs with COMPADRE variables (Meta-regression)
source("code/Jan_2025/S9_meta_regression_taxon_demographic_modelled_BACI.R")

# script 10 : Model DRAGNet Taxon effects from simple BACIs with COMPADRE variables (Meta-regression)
source("code/Jan_2025/S9_meta_regression_taxon_demographic_modelled_BACI.R")

# script 10 : Model DRAGNet Site / Taxon effects with COMPADRE and community vars (Meta-regression)
source("code/Jan_2025/S10_meta_regression_taxon_site_community.R")

################################################################################

# Extra information # 

# script E1: 
source("code/")

# script E2: Variance partitioning random effects models
source("code/Jan_2025/S11_variance_partitioning.R")

################################################################################

# Appendices / supporting information #

# script A1: get sample sizes with zeros (DRAGNet)
source("A1_get_sample_sizes_w_zeros.R")

# script A2: effect size contextualisation from simple BACI and GLMM models: Appendix 2
source("A2_random_effects_visualise.R")

# script A3: test if the different methods of 
source("A3_random_effects_visualise.R")

################################################################################

# Paper Figures #



# map out sites
# 
