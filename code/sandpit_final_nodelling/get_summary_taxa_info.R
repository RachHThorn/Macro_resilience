# R Thornley
# 15/07/2025
# Project: P1_COMPADRE_DRAGNET
# Script: get summary information for tables

library(tidyverse)

################################################################################

# read in the shared list of species found in both compadre and dragnet
Drag_taxa <- read_csv("results/July_2025/common_species_drag_comp.csv")
Drag_taxa <- Drag_taxa %>% pull(x)

effects_1 <- read_csv("results/July_2025/RE_SE_Taxon_all_DRAGNet.csv")
table_1a <-
  effects_1 %>% 
  filter(model == "Ordbeta") %>%
  dplyr::rename("Taxon" = "group") %>%
  filter(Taxon %in% Drag_taxa) %>%
  count(model, experiment, time_period) %>%
  rename(species_level_analysis = n)
table_1a
names(table_1a)

BACI_1 <- read_csv("results/July_2025/BACI_responses_simple_quadrat_all_DRAGNet.csv")
names(BACI_1)
head(BACI_1)
table_1b <-
  BACI_1 %>% 
  filter(Taxon %in% Drag_taxa) %>%
  group_by(trt, time) %>%
  mutate(taxon = as.character(Taxon)) %>%
  summarise(unique_taxa = n_distinct(Taxon))
names(table_1b)

table_1b <-
  table_1b %>% mutate(model = "BACI") %>% 
  rename(time_period = time, experiment = trt, species_level_analysis = unique_taxa)
names(table_1b)

# join the tables 
table <- rbind(table_1a, table_1b)
head(table)

# there seem to be the same number of taxa in the analysis no matter the approach -
# is this right or have I made an error somewhere?
# for now lets just leave it at that - but need looking at.

table %>% 
  select(time_period, experiment, model, species_level_analysis) %>%
  arrange(experiment, model, time_period)

# export the ordered beta model results only
write_csv(table_1a, "results/March_2025/number_species_per_analysis")


################################################################################

# number of effects for the community level analysis
# number of effects by species and site
effects_2 <- read_csv("results/July_2025/RE_SE_Taxon_site_all_DRAGNet.csv")
effects_2$group

head(effects)
table_2 <-
  effects_2 %>% 
  filter(model == "Ordbeta") %>%
  dplyr::rename("Taxon" = "group") %>%
  separate(Taxon, into = c("taxa_1", "taxa_2", "site")) %>%
  mutate(Taxon = paste0(taxa_1, "_", taxa_2)) %>%
  filter(Taxon %in% Drag_taxa) %>%
  mutate(taxa = paste0(taxa_1, "_", taxa_2)) %>%
  count(site, taxa)
table_2

unique(table_2$site)

################################################################################

# here I want to know how many sites are included in the analysis


# Load the REs from the GLMMs
effects %>%
  dplyr::rename("Taxon" = "group") %>%
  filter(Taxon %in% Drag_taxa)
# filter for just the ordered beta model and the disturbance experiment
plot_dat <- effects %>%
  filter(model == "Ordbeta") %>%
  filter(experiment == "DIST") 
unique(plot_dat$group) # 754 


