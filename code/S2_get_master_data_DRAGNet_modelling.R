# R Thornley 
# 28/10/2024
# get master DRAGNet data cleaned 
# Project: P1_COMPADRE_DRAGNET
# Script: S2_get_master_data_DRAGNET_modelling

# 1) Remove all non taxa entries from DRAGNet master data set
# 2) Filter for the list of sites that we have data for the time points we are interested in 
# 3) Add zeros so we have paired abundances for T0-T1 and T0-T2 and tidy up the data frame
# 4) Filter the DRAGNet data for the compadre overlap
# 5) Export the two data sets for teh two time frames for further modelling (T0-T1 and T0-T2)

library(tidyverse) 

################################################################################

# 1) Remove all non taxa entries from DRAGNet master data set

# create vector of none taxa entries in cover data
none_taxa <- c("Fungi", "Other_litter", "Other_standing_water", "Ground", "Bryophyte", 
               "Other_animal_diggings", "Other_woody_overstory", "Lichen",
               "Other_animal_droppings", "Other_rock",
               "Other_soil_biocrust", "Other_standing_dead")

# read in master raw data file and tidy taxon entries
drag <- read.csv("data/full-cover-2025-07-16.csv") %>%
  mutate(New_taxon = str_to_sentence(Taxon)) %>%
  mutate(New_taxon = str_replace_all(New_taxon, " ", "_")) %>%
  filter(!str_detect(New_taxon, ".sp")) %>% # get rid of entries not to taxon level
  filter(!str_detect(New_taxon, "_x_")) %>% # get rid of any hybrids
  filter(!str_detect(New_taxon, "Unknown")) %>% # get rid of unknown species
  filter(!New_taxon %in% none_taxa) %>% # get rid of non taxa entries 
  mutate(New_taxon = case_when(New_taxon == "Helianthemum_nummularium_var._Grandiflorum" ~ "Helianthemum_nummularium",
                               New_taxon == "Mimosa_quadrivalvis_var._Platycarpa" ~ "Mimosa_quadrivalvis",
                               New_taxon == "Sebaea_sedoides_var._Schoenlandii" ~ "Sebaea_sedoides", 
                               TRUE ~ New_taxon)) %>%
  arrange(New_taxon)

################################################################################

# T0-T1

# Filter for the list of sites that we have data for the time points we are interested in 
# use the table saved as 'dragnet blocks surveyed.csv' created by Script1
wanted <- read.csv("results/July_2025/dragnet_blocks_surveyed.csv", check.names = FALSE)
names(wanted) <- c("position", "site", "T1minus", "T0", "T1", "T2", "T3", "T4", "T5")

# select sites that have at least some data for T0 and T1 (before and after)
wanted <- wanted %>%
  select(site, T0, T1, T2, T3) %>%
  rowwise %>%
  filter(all(c_across(T0:T1) > 0)) %>%
  ungroup() %>%
  select(site, T0, T1)%>%
  pull(site)

# now filter the dragnet data with these sites
drag_1 <- drag %>% filter(site_name %in% wanted)
# filter with the 2 time points we want year_trt 0 and 1 for analysis 1
years_1 <- c("0", "1")
dat_1 <- drag_1 %>% filter(year_trt %in% years_1)

# replace all missing values with zero across the small dataset
dat_1 <- dat_1 %>% 
  select(site_name, New_taxon, block, trt, year_trt, max_cover)%>%
  group_by(site_name) %>% # group by site first
  complete(New_taxon, year_trt, trt, block) %>% # complete requires the unique combinations within each site
  replace(is.na(.),0)
dat_1

# for the first data set make sure the variables are ordered factors
# the cover values are percentages 
# the values of 100% are set to just less than this for modelling
dat_1 <- dat_1 %>%
  mutate(trt = factor(trt, levels = c("Control", "Disturbance", "NPK+Disturbance", "NPK", "NPK_Cessation"))) %>%
  filter(!trt == "NPK_Cessation") %>%
  mutate(year_trt = factor(year_trt, levels = c("0", "1"))) %>% 
  mutate(site_name = factor(site_name))%>%
  mutate(New_taxon = factor(New_taxon)) %>%
  mutate(taxon_site = paste0(New_taxon, "_", site_name)) %>% 
  mutate(taxon_site = factor(taxon_site)) %>%
  dplyr::select(site_name, New_taxon, taxon_site, year_trt, trt, max_cover, block) %>%
  mutate(new_max_cover = case_when(max_cover == 0 ~ 0,
                                   max_cover > 0 ~ max_cover/100,
                                   TRUE ~ max_cover)) %>%
  mutate(new_max_cover = if_else(new_max_cover >= 1, 0.99, new_max_cover))

# we are still left with some pairs of zeros so we need to filter these out
dat_1 <- dat_1 %>% group_by(site_name, New_taxon, taxon_site, trt, block) %>%
  filter(any(new_max_cover != 0)) %>%
  ungroup()

unique(dat_1$site_name) # 48 sites

# export this dataset - not filtered for the compadre species only
write.csv(dat_1, "results/July_2025/DRAGNet_all_clean_data_with_zeros_T0_T1.csv")

# 4) Filter the DRAGNet data for the compadre overlap
# list of shared species between dragnet and compadre as of OCT 2024
species <- 
  read.csv("results/July_2025/common_species_drag_comp.csv", header = TRUE) %>%
  pull(x)

# filter dragnet data for the species found in COMPADRE
dat_1 <- dat_1 %>% filter(New_taxon %in% species)
# look at the data we have to work with
dat_1 %>% ungroup %>% summarize(unique_values = n_distinct(site_name)) # 42 sites
dat_1 %>% ungroup %>% summarize(unique_values = n_distinct(New_taxon)) # 66 species

write.csv(dat_1, "results/July_2025/DRAGNet_overlap_clean_data_with_zeros_T0_T1.csv")

################################################################################

# T0-T2

# Filter for the list of sites that we have data for the time points we are interested in 
# use the table saved as 'dragnet blocks surveyed.csv' created by Script1
wanted <- read.csv("results/July_2025/dragnet_blocks_surveyed.csv", check.names = FALSE)
names(wanted) <- c("position", "site", "T1minus", "T0", "T1", "T2", "T3", "T4", "T5")
wanted
# select sites that have at least some data for T0 and T1 (before and after)
wanted <- wanted %>%
  select(site, T0, T1, T2, T3) %>%
  rowwise %>%
  filter(all(c_across(T0:T2) > 0)) %>%
  ungroup() %>%
  select(site, T0, T2)%>%
  pull(site)
wanted

# now filter the dragnet data with these sites
drag_2 <- drag %>% filter(site_name %in% wanted)
# filter with the 2 time points we want year_trt 0 and 2 for analysis 2
years_2 <- c("0", "2")
dat_2 <- drag_2 %>% filter(year_trt %in% years_2)

# replace all missing values with zero across the small dataset
dat_2 <- dat_2 %>% 
  select(site_name, New_taxon, block, trt, year_trt, max_cover)%>%
  group_by(site_name) %>% # group by site first
  complete(New_taxon, year_trt, trt, block) %>% # complete requires the unique combinations within each site
  replace(is.na(.),0)
dat_2

# for the second data set
dat_2 <- dat_2 %>%
  mutate(trt = factor(trt, levels = c("Control", "Disturbance", "NPK+Disturbance", "NPK", "NPK_Cessation"))) %>%
  filter(!trt == "NPK_Cessation") %>%
  mutate(year_trt = factor(year_trt, levels = c("0", "2"))) %>% 
  mutate(site_name = factor(site_name))%>%
  mutate(New_taxon = factor(New_taxon)) %>%
  mutate(taxon_site = paste0(New_taxon, "_", site_name)) %>% 
  mutate(taxon_site = factor(taxon_site)) %>%
  dplyr::select(site_name, New_taxon, taxon_site, year_trt, trt, max_cover, block) %>%
  mutate(new_max_cover = case_when(max_cover == 0 ~ 0,
                                   max_cover > 0 ~ max_cover/100,
                                   TRUE ~ max_cover)) %>%
  mutate(new_max_cover = if_else(new_max_cover >= 1, 0.99, new_max_cover))
range(dat_2$new_max_cover)

# remove any groups in the data that are all zeros
dat_2 <- dat_2 %>% group_by(site_name, New_taxon, taxon_site, trt, block) %>%
  filter(any(new_max_cover != 0)) %>%
  ungroup()

# export this data for modelling
write.csv(dat_2, "results/July_2025/DRAGNet_all_clean_data_with_zeros_T0_T2.csv")

# 4) Filter the DRAGNet data for the compadre overlap
# list of shared species between dragnet and compadre as of OCT 2024
species <- 
  read.csv("results/common_species_drag_comp.csv", header = TRUE) %>%
  pull(x)

# filter dragnet data for the species found in COMPADRE
dat_2 <- dat_2 %>% filter(New_taxon %in% species)
# look at the data we have to work with
dat_2 %>% ungroup %>% summarize(unique_values = n_distinct(site_name)) # 32 sites
dat_2 %>% ungroup %>% summarize(unique_values = n_distinct(New_taxon)) # 57 species

write.csv(dat_2, "results/July_2025/DRAGNet_overlap_clean_data_with_zeros_T0_T2.csv")

################################################################################

# T0-T3

# Filter for the list of sites that we have data for the time points we are interested in 
# use the table saved as 'dragnet blocks surveyed.csv' created by Script1
wanted <- read.csv("results/July_2025/dragnet_blocks_surveyed.csv", check.names = FALSE)
wanted$`3`
names(wanted) <- c("position", "site", "T1minus", "T0", "T1", "T2", "T3", "T4", "T5")
wanted
# select sites that have at least some data for T0 and T3 (before and after)
wanted <- wanted %>%
  select(site, T0, T1, T2, T3) %>%
  rowwise %>%
  filter(all(c_across(T0:T3) > 0)) %>%
  ungroup() %>%
  select(site, T0, T3)%>%
  pull(site)
wanted

# now filter the dragnet data with these sites
drag_2 <- drag %>% filter(site_name %in% wanted)
# filter with the 2 time points we want year_trt 0 and 2 for analysis 2
years_3 <- c("0", "3")
dat_3 <- drag_2 %>% filter(year_trt %in% years_3)

# replace all missing values with zero across the small dataset
dat_3 <- dat_3 %>% 
  select(site_name, New_taxon, block, trt, year_trt, max_cover)%>%
  group_by(site_name) %>% # group by site first
  complete(New_taxon, year_trt, trt, block) %>% # complete requires the unique combinations within each site
  replace(is.na(.),0)
dat_3

# for the third data set
dat_3 <- dat_3 %>%
  mutate(trt = factor(trt, levels = c("Control", "Disturbance", "NPK+Disturbance", "NPK", "NPK_Cessation"))) %>%
  filter(!trt == "NPK_Cessation") %>%
  mutate(year_trt = factor(year_trt, levels = c("0", "3"))) %>% 
  mutate(site_name = factor(site_name))%>%
  mutate(New_taxon = factor(New_taxon)) %>%
  mutate(taxon_site = paste0(New_taxon, "_", site_name)) %>% 
  mutate(taxon_site = factor(taxon_site)) %>%
  dplyr::select(site_name, New_taxon, taxon_site, year_trt, trt, max_cover, block) %>%
  mutate(new_max_cover = case_when(max_cover == 0 ~ 0,
                                   max_cover > 0 ~ max_cover/100,
                                   TRUE ~ max_cover)) %>%
  mutate(new_max_cover = if_else(new_max_cover >= 1, 0.99, new_max_cover))
range(dat_3$new_max_cover)

# remove any groups in the data that are all zeros
dat_3 <- dat_3 %>% group_by(site_name, New_taxon, taxon_site, trt, block) %>%
  filter(any(new_max_cover != 0)) %>%
  ungroup()

# export this data for modelling
write.csv(dat_3, "results/July_2025/DRAGNet_all_clean_data_with_zeros_T0_T3.csv")

# 4) Filter the DRAGNet data for the compadre overlap
# list of shared species between dragnet and compadre as of OCT 2024
species <- 
  read.csv("results/July_2025/common_species_drag_comp.csv", header = TRUE) %>%
  pull(x)
# 82 species that overlap

# filter dragnet data for the species found in COMPADRE
dat_3 <- dat_3 %>% filter(New_taxon %in% species)
# look at the data we have to work with
dat_3 %>% ungroup %>% summarize(unique_values = n_distinct(site_name)) # 19 sites
dat_3 %>% ungroup %>% summarize(unique_values = n_distinct(New_taxon)) # 51 species

write.csv(dat_3, "results/July_2025/DRAGNet_overlap_clean_data_with_zeros_T0_T3.csv")

# look at total number of sites we potentially have data for
T1_sites <- dat_1$site_name
T2_sites <- dat_2$site_name
T3_sites <- dat_3$site_name

both <- union(T1_sites, T2_sites) # 37 sites in T1 and T2
both <- union(both, T3_sites) # 42 sites in T1, T2, T3

write.csv(both, "results/July_2025/list_of_sites_in_analysis.csv")

################################################################################

# get list of sites used in study as site codes
drag <- read.csv("data/full-cover-2025-07-16.csv")
n_distinct(drag$site_code) # 56 sites
# filter for the sites where the COMPADRE species are found using the vector 'both'
drag <- drag %>% filter(site_name %in% both) %>% select(site_name, site_code)
drag <- drag %>% distinct(site_code, .keep_all = TRUE)
write_csv(drag, "results/July_2025/sites_used_codes.csv")

writeLines(drag$site_code, "results/July_2025/sites_used_codes.csv")
