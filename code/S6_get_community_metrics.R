# R Thornley
# 19/11/2024
# create community metrics to test for interactions with the demography metrics

# 1) Site species diversity at a site
# 2) Species native or non native status at a site
# 3) The average total cover per site per treatment and time point (proxy for competition)

rm(list = ls())

library(tidyverse) 
library(vegan)

################################################################################

# SITE LEVEL DIVERSITY # 

# we need to repeat the above for the other time point

# SITE LEVEL DIVERSITY # 

# if we do this at the quadrat scale - we need a list of quadrats that are 
# being used in the analysis
# read in the clean data 
data <- read.csv("results/DRAGNet_T0_T1_all.csv")
names(data)
# create a unique variable for each quadrat
data <- data %>% mutate(unique = paste(site_name, year_trt, trt, block, sep = "_"))
n_distinct(data$unique) # 1574 quadrats in the data set

# perform the diversity calculations
# to do this we need to pivot the data longer and create a matrix-like data frame
diversity <-  
  data %>%
  select(unique, max_cover, New_taxon) %>%
  pivot_wider(names_from = New_taxon, values_from = max_cover, values_fn = ~ mean(.x, na.rm = TRUE)) %>%
  replace(is.na(.), 0) %>%
  arrange(unique) %>%
  select(!unique) %>%
  mutate(shannon = vegan::diversity(., "shannon", MARGIN = 1)) %>%
  mutate(simpson = vegan::diversity(., "simpson", MARGIN = 1)) %>%
  mutate(invsimp = vegan::diversity(., "invsimpson", MARGIN = 1)) %>%
  select(shannon, simpson, invsimp)

# create species richness values for all the quadrats
richness <- data %>% 
  select(unique, max_cover, New_taxon) %>%
  group_by(unique) %>%
  mutate(richness = n_distinct(New_taxon)) %>% 
  select(unique, richness) %>%
  distinct(unique, richness)

# now join the quadrat level richness and diversiy metrics
all <- diversity %>% cbind(richness)

# now we have quadrat based diversity indices for our species context
# we may want to simplify these values later to mean species within site etc.

all <- all %>% 
  separate(unique, into = c("site", "time", "treatment", "block"), remove = FALSE, sep = "_") %>%
  mutate(site_time = paste(site, time, sep = "_")) %>% 
  group_by(site_time, treatment) %>%
  mutate(site_time_rich = max(richness)) %>%
  mutate(site_time_simp = mean(simpson)) %>%
  mutate(site_time_shan = mean(shannon)) %>%
  mutate(site_time_invsimp = mean(invsimp)) %>%
  ungroup() 
length(unique(all$site_time)) # 103 unique site time combinations

# tidy up some of the names for modelling later
all <- all %>% mutate(time_period = "T0-T1")

# summarise the data set to unique values only
# and keep only important columns
all <- all %>%
  distinct(site_time, site_time_simp, .keep_all = TRUE)

names(all)
T1_div <- all %>% select(site, time_period, treatment, block, unique, site_time_rich:site_time_invsimp)

################################################################################

# we need to repeat the above for the other time point

# SITE LEVEL DIVERSITY # 
# if we do this at the quadrat scale - we need a list of quadrats that are 
# being used in the analysis
# read in the clean data with the overlap of species
data <- read.csv("results/March_2025/DRAGNet_overlap_clean_data_with_zeros_T0_T2.csv")

wanted <- data %>% 
  mutate(trt = case_when(trt == "NPK_Cessation" ~ "NPK+Cesssation",
                         TRUE ~ trt)) %>%
  mutate(quadrats = paste(site_name, year_trt, trt, block, sep = "_")) %>%
  distinct(quadrats) %>%
  pull(quadrats) 
length(wanted)# 584 unique quadrats

# now get all our dragnet data
data <- read.csv("results/March_2025/DRAGNet_all_clean_data_with_zeros_T0_T2.csv")
data$year_trt
diversity <- data %>% 
  mutate(trt = case_when(trt == "NPK_Cessation" ~ "NPK+Cesssation",
                         TRUE ~ trt)) %>%
  mutate(unique = paste(site_name, year_trt, trt, block, sep = "_")) %>%
  filter(unique %in% wanted) %>%
  select(unique, max_cover, New_taxon) %>%
  pivot_wider(names_from = New_taxon, values_from = max_cover, values_fn = ~ mean(.x, na.rm = TRUE)) %>%
  replace(is.na(.), 0) %>%
  arrange(unique) %>%
  select(!unique) %>%
  mutate(shannon = vegan::diversity(., "shannon", MARGIN = 1)) %>%
  mutate(simpson = vegan::diversity(., "simpson", MARGIN = 1)) %>%
  mutate(invsimp = vegan::diversity(., "invsimpson", MARGIN = 1)) %>%
  select(shannon, simpson, invsimp)
# we have values for 558 quadrats
richness <- data %>% 
  mutate(trt = case_when(trt == "NPK_Cessation" ~ "NPK+Cesssation",
                         TRUE ~ trt)) %>%
  mutate(unique = paste(site_name, year_trt, trt, block, sep = "_")) %>%
  filter(unique %in% wanted) %>%
  select(unique, max_cover, New_taxon) %>%
  group_by(unique) %>%
  mutate(richness = n_distinct(New_taxon)) %>% 
  select(unique, richness) %>%
  distinct(unique, richness)

all <- diversity %>% cbind(richness)

# now we have quadrat based diversity indices for our species context
# we may want to simplify these values later to mean species within site etc.

all <- all %>% 
  separate(unique, into = c("site", "time", "treatment", "block"), remove = FALSE, sep = "_") %>%
  mutate(site_time = paste(site, time, sep = "_")) %>% 
  group_by(site_time, treatment) %>%
  mutate(site_time_rich = max(richness)) %>%
  mutate(site_time_simp = mean(simpson)) %>%
  mutate(site_time_shan = mean(shannon)) %>%
  mutate(site_time_invsimp = mean(invsimp)) %>%
  ungroup() 
length(unique(all$site_time)) # 48 unique site time combinations

# tidy up some of the names for modelling later
all <- all %>% mutate(time_period = "T0-T2")

# summarise the data set to unique values only
# and keep only important columns
all <- all %>%
  distinct(site_time, site_time_simp, .keep_all = TRUE)

T2_div <- all %>% select(site, treatment, time_period, site_time_rich:site_time_invsimp)

# we now have accompanying community data for each quadrat being used in our analysis
# with site based averages
both <- T1_div %>% rbind(T2_div)
write_csv(both, "results/March_2025/diversity_metrics_dragnet.csv")

################################################################################

# SPECIES PROVENANCE # 

# create species level metrics for dragnet focal species

library(tidyverse) 

# extract whether the species is native / non native
# we need to use the original raw data for this
# tidy the taxa names
# extract the site level variables for each taxa

# then we need to bind this to the tidy dataframes we made 

# GET NATIVE / NON-NATIVE RANGE

# Read in the raw DRAGNet data - we need this rather than the clean data as it has
# provenance column 
data <- read.csv("data/full-cover-drag-2024-07-27.csv")
names(data)

# create vector of none taxa entries in cover data
none_taxa <- c("Fungi", "Other_litter", "Other_standing_water", "Ground", "Bryophyte", 
               "Other_animal_diggings", "Other_woody_overstory", "Lichen",
               "Other_animal_digging", "Other_animal_droppings")

# create a tidy data frame of taxon name, site and whether native or non-native
native <- data %>% 
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
  dplyr::select(site_name, New_taxon, local_provenance)

# this now requires filtering as per the species in our compadre / dragnet data set
common <- read.csv("results/common_species_drag_comp.csv")
common <- common %>%
  pull(x)

# this is a list of the native and non-native status of the focal species across sites
all <- native %>% filter(New_taxon %in% common) %>% distinct(site_name, New_taxon, .keep_all = TRUE)

# check this makes sense

n_distinct(all$New_taxon) # 74 taxa
n_distinct(all$site_name) # 43 sites

# some of the provenances are either unknown or NULL
ggplot(all, aes(local_provenance))+
  geom_bar()

all %>% filter(local_provenance == "UNK")
# this can be corrected as we know that RB is native here

all %>% filter(local_provenance == "NULL")
# this can be corrected as we know that RB is native here

# there are 13 cases where NULL has been inserted
# we can look these up easily online

# Fragaria vesca - native to  Bayreuth DRAGNet
# Lactuca serriola - non-native in Blandy Experimental Farm
# Cerastium_fontanum - native in  CEREEP - Ecotron IDF
# Erigeron_canadensis - non native in CEREEP - Ecotron IDF 
# Alliaria_petiolata - native in CEREEP - Ecotron IDF 
# Achillea_millefolium - native in CEREEP - Ecotron IDF 
# Myosotis_ramosissima - native in CEREEP - Ecotron IDF
# Ranunculus_acris - native in Hazelrigg
# Verbascum_thapsus - non native in Piedmont Prairie
# Tragopogon_dubius - non native in University of Michigan Biological Station
# Tragopogon_dubius - non native in University of Wyoming - SAREC
# Bromus_tectorum - non native University of Wyoming - SAREC
# Latuca_serriola - non native University of Wyoming - SAREC

# now perform a case_when mutate to sort this data out
all$local_provenance
all$site_name
all$New_taxon
new <- all %>% 
  mutate(new_provenance = case_when(site_name == "Bayreuth DRAGNet" & New_taxon == "Fragaria_vesca" ~ "INT",
                                    site_name == "Blandy Experimental Farm" & New_taxon == "Lactuca_serriola" ~ "INT",
                                    site_name == "CEREEP - Ecotron IDF" & New_taxon == "Cerastium_fontanum" ~ "NAT",
                                    site_name == "CEREEP - Ecotron IDF" & New_taxon == "Erigeron_canadensis" ~ "INT",
                                    site_name == "CEREEP - Ecotron IDF" & New_taxon == "Alliaria_petiolata" ~ "NAT",
                                    site_name == "CEREEP - Ecotron IDF" & New_taxon == "Achillea_millefolium" ~ "NAT",
                                    site_name == "CEREEP - Ecotron IDF" & New_taxon == "Myosotis_ramosissima" ~ "NAT",
                                    site_name == "Hazelrigg" & New_taxon == "Ranunculus_acris" ~ "NAT",
                                    site_name == "Piedmont Prairie" & New_taxon == "Verbascum_thapsus" ~ "INT",
                                    site_name == "University of Michigan Biological Station" & New_taxon == "Tragopogon_dubius" ~ "INT",
                                    site_name == "University of Wyoming - SAREC" & New_taxon == "Tragopogon_dubius" ~ "INT",
                                    site_name == "University of Wyoming - SAREC" & New_taxon == "Bromus_tectorum" ~ "INT",
                                    site_name == "University of Wyoming - SAREC" & New_taxon == "Lactuca_serriola" ~ "INT",
                                    site_name == "Wytham Woods" & New_taxon == "Ranunculus_bulbosus" ~ "NAT",
                                    TRUE ~ local_provenance))

# some of the provenances are either unknown or NULL
ggplot(new, aes(new_provenance))+
  geom_bar()

# save data for next stage of modelling
write_csv(new, "results/March_2025/native_status.csv")

# visualise the proportion of data for native status for each of the species
new %>% group_by(New_taxon, new_provenance) %>% tally() %>%
  ggplot(aes(New_taxon, n, fill = new_provenance)) + 
  theme_bw()+
  geom_col() +
  theme(axis.text.x = element_blank())+
  xlab("Taxon")

#################################################################################

# TOTAL COVER #

rm(list=ls())

data <- read.csv("results/March_2025/DRAGNet_all_clean_data_with_zeros_T0_T1.csv")
head(data)
data <- data %>% 
  mutate(trt = case_when(trt == "NPK_Cessation" ~ "NPK+Cesssation",
                         TRUE ~ trt)) %>%
  select(site_name, year_trt:new_max_cover) %>%
  group_by(site_name, year_trt, trt, block) %>%
  summarise(total = sum(new_max_cover, na.rm = TRUE))%>%
  mutate(total = total *100)

# cover values per quadrat range from 27-641
range(data$total)
data <- data %>% filter(total < 400) # I think there is some sort of error here

# tidy up names
data <- data %>% mutate(experiment = case_when(trt == "Disturbance" ~ "DIST",
                                               trt == "NPK+Disturbance"  ~ "INTER", 
                                               trt == "NPK"  ~ "NPK", 
                                               trt == "Control" ~ "CONTROL"))

T1 <- data %>% filter(year_trt == "1") %>% rename(time_period = year_trt) %>% mutate(time_period = "T0-T1")
# now we have per quadrat cover values
# can we 
# write_csv(T1, "results/March_2025/total_cover_per_quadrat_T0_T1.csv")

# visualise these cover values
ggplot(T1, aes(x = reorder(site_name, -total), y = total))+
  theme_bw()+
  geom_boxplot()+
  facet_grid(experiment~time_period)+
  xlab("site")+
  theme(axis.text.x = element_text(angle = 90, hjust = 0.5))
# we still need to find a mean per treatment, site, time combination for modeling later


################################################################################

# year 2
data <- read.csv("results/March_2025/DRAGNet_all_clean_data_with_zeros_T0_T2.csv")
head(data)
data <- data %>% 
  mutate(trt = case_when(trt == "NPK_Cessation" ~ "NPK+Cesssation",
                         TRUE ~ trt)) %>%
  select(site_name, year_trt:new_max_cover) %>%
  group_by(site_name, year_trt, trt, block) %>%
  summarise(total = sum(new_max_cover,na.rm = TRUE))%>%
  mutate(total = total *100)

# cover values per quadrat range from 27-641
range(data$total)
data <- data %>% filter(total < 400) # I think there is some sort of error here

# tidy up names
data <- data %>% mutate(experiment = case_when(trt == "Disturbance" ~ "DIST",
                                               trt == "NPK+Disturbance"  ~ "INTER", 
                                               trt == "NPK"  ~ "NPK", 
                                               trt == "Control" ~ "CONTROL"))

T2 <- data %>% filter(year_trt == "2") %>% rename(time_period = year_trt) %>% mutate(time_period = "T0-T2")
# now we have per quadrat cover values
# can we 
# write_csv(T2, "results/March_2025/total_cover_per_quadrat_T0_T2.csv")

# visualise these cover values
ggplot(T2, aes(x = reorder(site_name, -total), y = total))+
  theme_bw()+
  geom_boxplot()+
  facet_grid(experiment~time_period)+
  xlab("site")+
  theme(axis.text.x = element_text(angle = 90, hjust = 0.5))

both <- T1 %>% rbind(T2)
both$trt <- NULL

write_csv(both, "results/March_2025/total_cover_per_quadrat.csv")

################################################################################
