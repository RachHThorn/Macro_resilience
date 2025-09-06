# R Thornley
# 16/07/2025
# Project: P1_COMPADRE_DRAGNET
# Script: Get summary data for tables

library(tidyverse)

# calculate simple BACI equation for each of the experiments 
# for both the time points

# NOTE: This needs tidying as the code is dirty at the mo

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
Drag_taxa <- read_csv("results/common_species_drag_comp.csv")
Drag_taxa <- Drag_taxa %>% pull(x)
# Load the REs from the first round of modelling and tidy
effects <- read_csv("results/July_2025/RE_SE_Taxon_all_DRAGNet.csv") %>%
  dplyr::rename("Taxon" = "group") %>%
  filter(Taxon %in% Drag_taxa)

# join both data sets
both <- demo %>% left_join(effects)
names(both)
unique(both$Taxon)

Nos_taxa <-
  both %>% filter(model == "Ordbeta") %>% 
  group_by(time_period, experiment, Demo_trait) %>%
  summarise(n_taxa = n_distinct(Taxon), .groups = "drop")

Nos_taxa <- Nos_taxa %>%
  mutate(Demographic_Trait = case_when(Demo_trait == "FirstStepAtt" ~ "First Step Attenuation",
                                    Demo_trait == "Lambda" ~ "Population Growth Rate",
                                    Demo_trait == "Lmax" ~ "Maximum Length of Life",
                                    Demo_trait == "Lmean" ~ "Mean Life Expectancy",
                                    Demo_trait == "R0" ~ "Per Capita Reproduction",
                                    Demo_trait == "Reactivity" ~ "Reactivity",
                                    Demo_trait == "T_generation" ~ "Generation Time",
                                    Demo_trait == "age_repro" ~ "Age at first reproduction",
                                    Demo_trait == "percapita_repro" ~ "Per capita reproduction",
                                    TRUE ~ Demo_trait)) %>%
  mutate(Experiment = case_when(experiment == "DIST" ~ "Disturbance",
                                experiment == "NPK" ~ "NPK",
                                experiment == "INTER" ~ "Disturbance x NPK")) %>%
  mutate(Experiment = factor(Experiment, levels = c("Disturbance", "NPK", "Disturbance x NPK"))) %>%
  select(time_period, Experiment, Demographic_Trait, n_taxa) %>%
  rename(Time_period = time_period, Number_taxa = n_taxa)
levels(Nos_taxa$Experiment)

max(Nos_taxa$Number_taxa) # 30
min(Nos_taxa$Number_taxa) # 18


write_csv(Nos_taxa, "results/July_2025/number_of_taxa_per_analysis.csv")
  
################################################################################

# Load the REs from the first round of modelling and tidy
effects <- read_csv("results/July_2025/RE_SE_Taxon_site_all_DRAGNet.csv") 
head(effects)
effects <-
  effects %>%
  filter(model == "Ordbeta") %>%
  separate(group, into = c("Taxon_1", "Taxon_2", "Site_1", "Site_2"), sep = "_") %>%
  mutate(Taxon = paste0(Taxon_1, "_", Taxon_2)) %>%
  filter(Taxon %in% Drag_taxa)
head(effects)
n_distinct(effects$Site_1) # 41 sites of DRAGNet 

################################################################################

# get the number of sites that have done 5 years of data collection already 
# (ie heva started on the cessation part of the experiment)

blocks <- read_csv("results/dragnet_blocks_surveyed.csv")
# at year 5 only one site has data recorded at the moment (Potrok Aike)

###############################################################################

# get the number of matrices for each species 
demo <- read_csv("results/July_2025/all_compadre_metrics.csv")
head(demo)
n_distinct(demo$SpeciesAccepted)
n_distinct(demo$MatrixID) # 2962
drag <- demo %>% filter(DRAGNet == "TRUE") 
n_distinct(drag$SpeciesAccepted) # 36 species 
table <- 
  drag %>% group_by(SpeciesAccepted, demo_var) %>% 
  count() %>% 
  rename(Species = SpeciesAccepted, Demographic_metric = demo_var, 
         Number_of_matrices = n) %>%
  filter(Demographic_metric == "Lambda")
table

ggplot(table, aes(reorder(Species, Number_of_matrices), Number_of_matrices))+
  geom_col()+
  theme_bw()+
  xlab("Species Name")+
  ylab("Number of Matrices")+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))
ggsave("figures/July_2025/Number_matrices_per_species.jpeg", height = 5, width = 7)
  
###############################################################################

# we should probably show the results for the spread of the demo metrics 
# according to the data set we pass onto the actual analysis

dat <- read_csv("results/July_2025/effects_demo_all_times_joined.csv")
head(dat)
names(dat)
dat %>% count(Taxon, Demo_trait, model, time_period, experiment)

unique(dat$model)
unique(dat$experiment)
unique(dat$time_period)
unique(dat$Demo_trait)

dat %>% 
  filter(model == "Ordbeta") %>% 
  filter(experiment == "DIST") %>%
  filter(time_period == "T0-T1") %>%
  filter(Demo_trait == "Lambda") %>%
  ggplot(aes(reorder(Taxon, Demo_value), Demo_value)) + 
  theme_bw()+
  geom_col() + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))+
  xlab("taxon")+
  ylab("population growth rate")

unique(dat$Demo_trait)
dat %>% 
  filter(model == "Ordbeta") %>% 
  filter(experiment == "DIST") %>%
  filter(time_period == "T0-T1") %>%
  filter(Demo_trait == "Lmax") %>%
  ggplot(aes(reorder(Taxon, Demo_value), Demo_value)) + 
  theme_bw()+
  geom_col() + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))+
  xlab("taxon")+
  ylab("Maximum Life Expectancy")

dat %>% 
  filter(model == "Ordbeta") %>% 
  filter(experiment == "DIST") %>%
  filter(time_period == "T0-T1") %>%
  filter(Demo_trait == "percapita_repro") %>%
  ggplot(aes(reorder(Taxon, Demo_value), Demo_value)) + 
  theme_bw()+
  geom_col() + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))+
  xlab("taxon")+
  ylab("Per capita reproduction")

dat %>% 
  filter(model == "Ordbeta") %>% 
  filter(experiment == "DIST") %>%
  filter(time_period == "T0-T1") %>%
  filter(Demo_trait == "Lmean") %>%
  ggplot(aes(reorder(Taxon, Demo_value), Demo_value)) + 
  theme_bw()+
  geom_col() + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))+
  xlab("taxon")+
  ylab("Mean Life Expectancy")


################################################################################

meta <- read_csv("results/March_2025/meta_data_compadre.csv")
head(meta)
wanted <- drag$MatrixID
meta <- meta %>% filter(MatrixID %in% wanted)
n_distinct(meta$MatrixPopulation) # 76 sites

