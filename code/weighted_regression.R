# R Thornley
# 29/08/2025
# Regressions for the Life History Traits and RE outcomes

library(tidyverse)
library(purrr)
library(broom)
library(brms)

################################################################################
# load and tidy data
################################################################################

# Load demography variables and find the mean and SE of each Taxon and trait

# we see from plotting out these values that there are still a lot of outliers 
# we are going to keep these at the moment
# but some of the variables have very skewed distributions 
# these should be corrected by logging
# see script A4_spread_compadre_data for visualisations

needs_logging <- c("R0", "percapita_repro", "age_repro", "Lmean", "Lmax")

demo <- read_csv("results/all_COMPADRE_metrics.csv") %>% 
  filter(DRAGNet == TRUE) %>% 
  rename("Taxon" = "SpeciesAccepted", "Demo_trait" = "demo_var", "Demo_value" = "value") %>%
  mutate(Taxon = str_replace_all(Taxon, " ", "_")) %>%
  select(Taxon, Demo_trait, Demo_value) %>%
  group_by(Taxon, Demo_trait) %>%
  mutate(Sample_size = n(),
         s = sd(Demo_value, na.rm = TRUE),
         se = s/sqrt(Sample_size),
         Demo_value_mean = mean(Demo_value, na.rm = TRUE)) %>%
  mutate(Taxon_label = str_replace_all(Taxon, "_", " ")) %>%
  mutate(Taxon_label = paste0(Taxon_label, " n = ", Sample_size)) %>%
  mutate(New_demo_value = if_else(Demo_trait %in% needs_logging, log(Demo_value), Demo_value)) %>%
  mutate(Log_flag = if_else(Demo_trait %in% needs_logging, "logged", "not_logged"))

# n shows us the number of matrices we are using for each data point
demo %>% group_by(Taxon, Demo_trait) %>% count()

# save a vector of taxa names to use in the next stage of modelling
comp_taxa <- unique(demo$Taxon) # 41 taxa

# read in effect sizes from the whole of DRAGNet 
taxa <- read_csv("results/RE_SE_Taxon_all_DRAGNet.csv")
unique(taxa$time_period)
# read in the taxon RE from the compadre and dragnet overlap
# comp_taxa <- read_csv("results/common_species_drag_comp.csv") %>% pull(x)

# select the random effects for the var of interest only (New_taxon)
# create a new ID variable the makes a note if data are in the compadre data set
# change some labels in the data frame to tidy it up
taxa <- 
  taxa %>% 
  filter(str_detect(group_var, "New_taxon")) %>%
  mutate(ID = if_else(group %in% comp_taxa, "1", "0")) %>%
  mutate(time_period = case_when(time_period == 1 ~ "T0-T1",
                                 time_period == 2 ~ "T0-T2",
                                 time_period == 3 ~ "T0-T3")) %>%
  dplyr::select(group, value, se, model, time_period, experiment, ID) %>%
  rename(Taxon = group, RE = value, RE_se = se) %>%
  filter(ID == "1")
names(taxa)
unique(taxa$Taxon) # 40 unique taxa here


################################################################################
# OPTION 1: Easy linear models
################################################################################

# quick and dirty linear regression using just the means of the demo metrics and REs

# first simplify the demo data to just the mean values
small_demo <- 
  demo %>% 
  group_by(Taxon, Demo_trait) %>%
  summarise(New_demo_value = mean(New_demo_value, na.rm = TRUE))

# check there is only one row per species / demo_trait combo 
small_demo %>% 
  group_by(Taxon, Demo_trait) %>% 
  count() %>%
  filter(n > 1) 

# join both data sets
both <- taxa %>% left_join(small_demo, by = "Taxon")
names(both)
unique(both$model)
unique(both$time_period)
unique(both$experiment)

# very rough crude linear regressions of mean values

both %>% group_by(model, time_period, experiment, Demo_trait) %>% count()
nested_both <- both %>% group_by(model, time_period, experiment, Demo_trait) %>% nest()

models <- 
  nested_both %>%
  mutate(data = map(data, ~ .x %>% 
                      filter(if_all(everything(), ~ is.finite(New_demo_value))) %>%  # keep only finite values
                      drop_na(New_demo_value))) %>%
  mutate(mod = map(data, ~ lm(RE ~ New_demo_value, data = .x)), 
         results = map(mod, tidy))

lm_results <- models %>%
  select(time_period, experiment, model, Demo_trait, results) %>%
  unnest(results) %>%
  filter(term == "New_demo_value")

# visualise these model results (p values)
lm_results %>% 
  filter(model == "Ordbeta") %>%
  ggplot(aes(p.value, Demo_trait)) + 
  theme_bw() +
  geom_point() + 
  facet_grid(time_period ~ experiment)

# also look at the Hurdle model results
lm_results %>% 
  filter(model == "Hurdle_cond") %>%
  ggplot(aes(p.value, Demo_trait)) + 
  geom_point() + 
  facet_grid(time_period ~ experiment)

lm_results %>%
  filter(model == "Ordbeta") %>%
  ggplot(aes(estimate, Demo_trait)) + 
  theme_bw()+
  geom_point() +
  geom_errorbar(aes(xmin = estimate - std.error, xmax = estimate + std.error), size = 0.3) +
  facet_grid(time_period ~ experiment)+
  geom_vline(xintercept = 0, colour = "red", linetype = "dashed")
# so mainly no relationships / there are some values not overlapping zero
# particularly Lambda


# the effects of Lambda are not help under the conditional Hurdle model as strongly
lm_results %>%
  filter(model == "Hurdle_cond") %>%
  ggplot(aes(estimate, Demo_trait)) + 
  theme_bw()+
  geom_point() +
  geom_errorbar(aes(xmin = estimate - std.error, xmax = estimate + std.error), size = 0.3) +
  facet_grid(time_period ~ experiment)+
  geom_vline(xintercept = 0, colour = "red", linetype = "dashed")

lm_results$p.value
lm_results %>% filter(p.value < 0.1)

# if we plot the actual data and look at the plots for the Lambda variable at T0-T3
# data is called both
names(both)
both %>% 
  filter(Demo_trait == "Lambda", model == "Ordbeta") %>%
  ggplot(aes(New_demo_value, RE))+
  theme_bw()+
  geom_point()+
  geom_smooth(method = "lm")+
  facet_grid(time_period ~ experiment)+
  xlab("Lambda")+
  ylab("Species cover change effect")
# this relationship is largely driven by an outlier
# save this plot as evidence
ggsave("figures/Lambda_mean_data_plot.jpeg", height = 5, width = 7)

unique(both$Demo_trait)
both %>% 
  filter(Demo_trait == "FirstStepAtt", model == "Ordbeta") %>%
  ggplot(aes(New_demo_value, RE))+
  geom_point()+
  geom_smooth(method = "lm")+
  facet_grid(time_period ~ experiment)+
  
# this is also very unconvincing



################################################################################
# OPTION 2: Weighted regression
################################################################################
# weights as inverse variance
names()
df <- df %>%
  mutate(weight = 1 / se^2)  # se = measurement standard error

fit <- lm(y ~ x, data = df, weights = weight)




################################################################################
# OPTION 3: Bayesian regression
################################################################################

# here we need a data frame with errors for both 




models_bayes <- nested_both %>%
  mutate(
    data = map(data, ~ .x %>%
                 filter(is.finite(New_demo_value), is.finite(RE)) %>%  # clean numeric columns
                 drop_na(New_demo_value, RE))),
mod = map(data, ~ brm(
  formula = RE ~ New_demo_value, 
  data = .x,
  family = gaussian(),
  chains = 2, iter = 2000, cores = 2,  # adjust for speed/memory
  silent = TRUE)), results = map(mod, ~ broom.mixed::tidy(.x, effects = "fixed")))