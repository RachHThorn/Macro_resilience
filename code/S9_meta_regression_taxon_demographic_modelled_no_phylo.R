# R Thornley
# 22/07/2025
# meta-regression for BACI results: modelled and simple with demographic metrics
# no phylogenetic dependency

library(tidyverse)
library(meta)
library(metafor)
library(orchaRd)

################################################################################

# join the effect sizes 

# Load demography variables
demo <- read_csv("results/July_2025/all_COMPADRE_metrics.csv") %>% 
  filter(DRAGNet == TRUE) %>% 
  rename("Taxon" = "SpeciesAccepted", "Demo_trait" = "demo_var", "Demo_value" = "value") %>%
  mutate(Taxon = str_replace_all(Taxon, " ", "_")) %>%
  select(Taxon, Demo_trait, Demo_value) %>% 
  group_by(Taxon, Demo_trait) %>%
  summarise(Demo_value = mean(Demo_value, na.rm = TRUE)) %>%
  drop_na(Demo_value)

# read in the shared list of species found in both compadre and dragnet
Drag_taxa <- read_csv("results/July_2025/common_species_drag_comp.csv")
Drag_taxa <- Drag_taxa %>% pull(x)
# Load the REs from the first round of modelling / tidy names / filter for the dragnet taxa
effects <- read_csv("results/July_2025/RE_SE_Taxon_all_DRAGNet.csv") %>%
  dplyr::rename("Taxon" = "group") %>%
  filter(Taxon %in% Drag_taxa)

# join both data sets
both <- demo %>% left_join(effects)
names(both)
head(both)

# save this data set
write_csv(both, "results/July_2025/effects_demo_all_times_joined.csv")

################################################################################

# now apply the meta-regression models

# create a nested data frame
nested_df <- 
  both %>% 
  drop_na() %>% 
  rename(taxon_estimate = value) %>%
  group_by(model, time_period, experiment, Demo_trait) %>% 
  nest() 

# if you want to test function on all of the data
# data <- nested_df %>% pull(data) %>% pluck(1)

################################################################################
# function that performs a meta-regression for each of the data frames contained 
# in a nested df object with time_period, experiment, model, demographic variable
# as grouping variables

perform_meta <- function(data) {
  
  meta_mod <- rma(yi = taxon_estimate,
                  sei = se,
                  data = data,
                  method = "ML",
                  mods = ~ Demo_value,
                  test = "knha")
  result <- 
    meta_mod %>%
    permutest() %>% 
    coef(summary()) 
  result$test <- rownames(result)
  result <- result %>% filter(test == "Demo_value")
  return(result)
  
} 

################################################################################

# run the perform_meta function across each of the nested dfs
result_modelled <- 
  nested_df %>%
  mutate(result = map(data, perform_meta)) 

# tidy into suitable data frame
tidy_df_modelled <- result_modelled %>% 
  dplyr::select(Demo_trait, model, time_period, experiment, result) %>%
  unnest(cols = c(result)) 

# write_csv(tidy_df_modelled, "results/Jan_2025/meta_regression_results_Taxon_all.csv")

###############################################################################

# visualise the meta regression results for each experiment / time / demo var
# for each model type

# effect size plus CI for each demo var
# Hurdle conditional model
tidy_df %>% 
  dplyr::filter(model == "Hurdle_cond") %>%
  ggplot(aes(estimate, Demo_trait))+
  theme_bw()+
  geom_point()+
  geom_errorbar(aes(xmax = ci.lb, xmin = ci.ub))+
  geom_vline(xintercept = 0, colour = "red")+
  facet_grid(time_period ~ experiment)+
  ggtitle("Hurdle conditional model meta regression results")

# Hurdle zero-inflated
tidy_df %>% 
  dplyr::filter(model == "Hurdle_zi") %>%
  ggplot(aes(estimate, Demo_trait))+
  theme_bw()+
  geom_point()+
  geom_errorbar(aes(xmax = ci.lb, xmin = ci.ub))+
  geom_vline(xintercept = 0, colour = "red")+
  facet_grid(time_period ~ experiment)+
  ggtitle("Hurdle zi model meta regression results")

tidy_df %>% 
  dplyr::filter(model == "Ordbeta") %>%
  ggplot(aes(estimate, Demo_trait))+
  theme_bw()+
  geom_point()+
  geom_errorbar(aes(xmax = ci.lb, xmin = ci.ub))+
  geom_vline(xintercept = 0, colour = "red")+
  facet_grid(time_period ~ experiment)+
  ggtitle("Frequentist Ordered Beta model meta regression results")

################################################################################

# now do the following for the simple BACIs
# Run these models on the simple BACIs to see if the results are different
BACI <- read_csv("results/March_2025/BACI_responses_simple_quadrat.csv")
# get taxon level means and SEs across all sites
head(BACI)
unique(BACI$time)
BACI <- 
  BACI %>%
  pivot_longer(cols = c(BACI_dist:BACI_inter), names_to = "experiment", values_to = "simple_BACI") %>%
  mutate(experiment = case_when(experiment == "BACI_dist" ~ "DIST",
                                experiment == "BACI_npk" ~ "NPK",
                                experiment == "BACI_inter" ~ "INTER")) %>%
  rename(time_period = time)

unique(BACI$Taxon) # 66 taxa that we have simple BACI effect sizes for

# create means and SE at the taxon level to model in the meta-regressions
# we need to filter the data where we only have one data point or else we can't find an error term
BACI <- 
  BACI %>% 
  group_by(Taxon, experiment, time_period) %>% 
  summarise(
    count = n(),  # Number of observations
    taxon_estimate = mean(simple_BACI, na.rm = TRUE), 
    se = ifelse(count > 1, sd(simple_BACI, na.rm = TRUE) / sqrt(count), NA),  # SE only if more than 1 observation
    .groups = "drop"
  )
# drop the taxon with no SE
BACI <- BACI %>% drop_na() 
BACI <- BACI %>% select(!count)
unique(BACI$Taxon) # 59 taxa that we can get effect sizes and SE for

# join both data sets
both <- demo %>% left_join(BACI)
names(both)

# create a nested data frame
nested_df <- 
  both %>% 
  drop_na() %>%
  group_by(time_period, experiment, Demo_trait) %>% 
  nest() 

# run the perform_meta function across each of the nested dfs
result_BACI <- 
  nested_df %>%
  mutate(result = map(data, perform_meta)) 

# tidy into suitable data frame
tidy_df_BACI <- result_BACI %>% 
  dplyr::select(Demo_trait, time_period, experiment, result) %>%
  unnest(cols = c(result)) 

tidy_df_BACI %>% 
  ggplot(aes(estimate, Demo_trait))+
  theme_bw()+
  geom_point()+
  geom_errorbar(aes(xmax = ci.lb, xmin = ci.ub)) +
  facet_grid(time_period ~ experiment)+
  geom_vline(xintercept = 0, colour = "red", linetype = "dashed")+
  ggtitle("Simple BACI Effects")

###############################################################################

# join all results from simple BACI and from the modelled BACIs
all <- tidy_df_BACI %>% mutate(model = "simple_BACI") %>%
  rbind(tidy_df_modelled)
write_csv(all, "results/March_2025/meta_regressions_no_phylo_results.csv")

