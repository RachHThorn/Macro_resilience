# R Thornley
# 28/01/2024
# create focal species strategy
# create neighbourhood species strategy

# get list of DRAGNet species and desirable traits
# source the necessary leaf traits from TRY
# put this data into the Pierce et al. 2017 spreadsheet (StrateFy) to get CSR per species
# use these values to get a measure of community CSR distance using the function 'specialization' from Ricotta et al. 2023

library(tidyverse) 
# install.packages("rtry")
library(rtry) # package to load in TRY data easily
library(ggtern) # plotting CSR strategies in ternary plot

###############################################################################

# TRY DATA #

# Load in the data ordered through TRY
# use the rtry package to load and query data
# data sets are quite large (GBs) so are saved in my documents file not in project folder at the moment

dat_Nov <- rtry::rtry_import("data/37407.txt")
dat_Jan <- rtry_import("data/38952.txt")
dat <- rbind(dat_Nov, dat_Jan)
dat <- dat %>% arrange(TraitID)

n_distinct(dat$AccSpeciesID) # 928 species - DRAGNet / TRY overlap
unique(dat$AccSpeciesName)
unique(dat$TraitID) # this is a trait number
unique(dat$TraitName) # this is a trait descriptor

IDs <- dat %>% select(TraitID) %>% drop_na() %>% arrange(TraitID) %>% pull(TraitID)
unique(IDs) # we have selected 179 traits

# for the CSR strategies we need 3 leaf traits
# LDMC
# Leaf Area
# SLA

# look at a list of the trait names
unique(dat$TraitName)

# Leaf Dry Mass per Leaf fresh mass (LDMC)
LDMC <- dat %>% 
  dplyr::select(AccSpeciesName, AccSpeciesID, ObservationID, TraitID, TraitName, StdValue, UnitName) %>%
  filter(TraitName == "Leaf dry mass per leaf fresh mass (leaf dry matter content, LDMC)") %>%
  mutate(StdValue = as.numeric(StdValue))%>%
  group_by(AccSpeciesName) %>% 
  mutate(LDMC_mean = mean(StdValue)) %>%
  drop_na()
LDMC$UnitName # g/g - ratio or a percentage

LDMC_mean <- LDMC %>%
  group_by(AccSpeciesName) %>%
  summarise(LDMC_mean = mean(LDMC_mean))
# 546 obs

# this needs to be converted to a precentage not a ratio
LDMC_mean$LDMC_mean <- LDMC_mean$LDMC_mean*100

# Specific Leaf Area (SLA)
# this trait is captured by a few different categories
wanted_SLA <- c("Leaf area per leaf dry mass (specific leaf area, SLA or 1/LMA): petiole excluded",
                "Leaf area per leaf dry mass (specific leaf area, SLA or 1/LMA): petiole included",
                "Leaf area per leaf dry mass (specific leaf area, SLA or 1/LMA): undefined if petiole is in- or excluded") 
SLA <- dat %>% 
  dplyr::select(AccSpeciesName, AccSpeciesID, ObservationID, TraitID, TraitName, StdValue, UnitName) %>%
  filter(TraitName %in% wanted_SLA) %>%
  mutate(StdValue = as.numeric(StdValue))%>%
  group_by(AccSpeciesName) %>% 
  mutate(SLA_mean = mean(StdValue)) %>%
  drop_na()
SLA$UnitName # mm2 mg-1

SLA_mean <- SLA %>%
  group_by(AccSpeciesName) %>%
  summarise(SLA_mean = mean(SLA_mean))
# 622 obs

# Leaf Area codes (1, 3108, 3109, 3117, 3110, 3111, 3112, 3113, 3114)
wanted_LA <- c("Leaf area (not yet refined if leaf or leaflet)",
            "Leaf area (in case of compound leaves: leaf, petiole excluded)",
            "Leaf area (in case of compound leaves: leaflet, petiole and rachis excluded)",
            "Leaf area (in case of compound leaves: leaf, petiole included)",
            "Leaf area (in case of compound leaves: leaflet, petiole and rachis included)", 
            "Leaf area (in case of compound leaves: leaf, undefined if petiole in- or excluded)", 
            "Leaf area (in case of compound leaves: leaflet, undefined if petiole and rachis are in- or excluded)",
            "Leaf area (in case of compound leaves undefined: leaf or leaflet; petiole and rachis in- excluded)")
LA <- dat %>% 
  dplyr::select(AccSpeciesName, AccSpeciesID, ObservationID, TraitID, TraitName, StdValue, UnitName) %>%
  dplyr::filter(TraitName %in% wanted_LA) %>%
  mutate(StdValue = as.numeric(StdValue))%>%
  group_by(AccSpeciesName) %>% 
  mutate(LA_mean = mean(StdValue)) %>%
  drop_na()
LA$UnitName # mm2

LA_mean <- LA %>%
  group_by(AccSpeciesName) %>%
  summarise(LA_mean = mean(LA_mean))
# 514 obs

all <- LDMC_mean %>% left_join(SLA_mean) %>% 
  left_join(LA_mean) %>% drop_na() %>% # some of the species don't have data for all 3 categories
  dplyr::select(AccSpeciesName, LA_mean, LDMC_mean, SLA_mean)
# export this data for use with the excel spreadsheet StrateFy
write_csv(all, "results/March_2025/TRY_DRAGNet_StateFy_data_more.csv")
# 337 species with values for the 3 necessary traits

################################################################################

# CSR STRATEGIES FROM StrateFy #
# the file exported below was used in the StateFy Excel method 
# out side this code to get species strategies
# Look at strategies and overlap with compadre

# read in the CSR strategies as calculated by the StateFy spreadsheet
csr <- read_csv("results/March_2025/CSR_results_StateFy.csv")
names(csr)

# visualise these strategies
# the ternany values are all clustered around the S corner of the triangle
ggtern(data = csr, aes(x = percent_r, y = percent_c, z = percent_s)) + 
  geom_point() +
  ggtitle("CSR for 428 species in DRAGNet")
# colour the points by csr categorisation
ggtern(data = csr, aes(x = percent_r, y = percent_c, z = percent_s, colour = `Strategy class`)) + 
  geom_point() +
  theme_bw()+
  ggtitle("CSR for 428 species in DRAGNet")

# look at the representation of our focal species in the data set (ie the overlap between compadre and dragnet)
focal <- read_csv("results/March_2025/common_species_drag_comp.csv") 
dim(focal) # 74 species common between the two data sets
focal <- focal %>% pull(x) # get the vector of names
focal <- str_replace(focal, "_", " ")
focal # 74 species overlap
focal <- csr %>% filter(AccSpeciesName %in% focal)
dim(focal)
# for 55 of my focal species I have csr strategies for

# get a plot of those 55 strategies
ggtern(data = focal, aes(x = percent_r, y = percent_c, z = percent_s, colour = `Strategy class`)) + 
  geom_point() +
  theme_bw()+
  ggtitle("CSR for 55 focal species with compadre/dragnet overlap")

################################################################################

# DREGREE OF CSR SPECIALIASATION #

# to use the 'specialization' function from Ricotta et al. we need two matrices
# comm (a table of species abundances per quadrat) taxon are colnames / plots are rownames
# CSR table (a table of the csr values) with taxon as rownames; cols c, s, r values from Pierce et al. 

# make the CSR table object

# read in the CSR strategies as calculated by the StateFy spreadsheet
csr <- read_csv("results/March_2025/CSR_results_StateFy.csv")
names(csr)

# tidy the csr data frame into a matrix in the right format
CSRtable <- 
  csr %>%
  dplyr::select(AccSpeciesName, percent_c, percent_s, percent_r) %>%
  rename("Taxon" = AccSpeciesName, "C" = percent_c, "S" = percent_s, "R" = percent_r) %>%
  as.data.frame()
csr_taxa <- CSRtable$Taxon # get a vector of names for filtering the DRAGNet data with
rownames(CSRtable) <- CSRtable$Taxon
CSRtable$Taxon <- NULL
CSRtable <- as.matrix(CSRtable)

################################################################################

# make the comm object

# load dragnet data
drag_dat <- read.csv("data/full-cover-drag-2024-07-27.csv")
head(drag_dat)
names(drag_dat)
unique(drag_dat$site_name)
# we want the species data for each quadrat
# create a unique quadrat ID
drag_dat$ID <- paste(drag_dat$site_name, drag_dat$block, drag_dat$plot, drag_dat$subplot, drag_dat$trt, drag_dat$year_trt, sep = "_")
drag_dat$trt
# iterate over all sites
comm <-
  drag_dat %>% 
  filter(functional_group != "BRYOPHYTE") %>%
  mutate(Taxon = str_to_sentence(Taxon)) %>%
  filter(!trt == "NPK_Cessation") %>%
  dplyr::select(site_name, ID, Taxon, max_cover) %>%
  filter(Taxon %in% csr_taxa) %>%
  group_by(site_name) %>%
  nest()

# test function with one data frame
# comm <- comm %>% pull(data) %>% pluck(1)

###############################################################################

# load the specialization function
source("code/March_2025/specialization.R")

################################################################################

get_special_1 <- function(comm, CRStable) {
  comm <-
    comm %>%
    pivot_wider(names_from = "ID", values_from = "max_cover") %>%
    mutate(across(where(is.numeric), ~replace_na(., 0))) %>%
    as.data.frame() 
  rownames(comm) <- comm$Taxon
  comm$Taxon <- NULL
  # make sure the df is all numeric
  comm <- comm %>% mutate_all(as.numeric)
  # transpose and make a matrix
  comm <- t(comm)
  special <- specialization(comm, CSRtable, tol = 1e-8)
  dat <- special$sigmaCom # a df values for each community 
  dat$ID <- rownames(dat)
  return(dat)
}

# create the get_special data frame for each species
special_1 <- comm %>%
  mutate(comm_dat = map(data, ~ mutate(get_special_1(.x, CRStable)))) %>%
  select(site_name, comm_dat) %>%
  unnest(cols = c(comm_dat))

write_csv(special_1, "results/March_2025/strategy_plot_level_values.csv")

################################################################################
         
get_special_2 <- function(comm, CRStable) {
  comm <-
    comm %>%
    pivot_wider(names_from = "ID", values_from = "max_cover") %>%
    mutate(across(where(is.numeric), ~replace_na(., 0))) %>%
    as.data.frame() 
  rownames(comm) <- comm$Taxon
  comm$Taxon <- NULL
  # make sure the df is all numeric
  comm <- comm %>% mutate_all(as.numeric)
  # transpose and make a matrix
  comm <- t(comm)
  special <- specialization(comm, CSRtable, tol = 1e-8)
  dat <- special$sigmaSpe
  dat <- as.data.frame(dat)
  dat$Taxon <- rownames(dat)
  return(dat)
}

# iterate over all sites
comm <-
  drag_dat %>% 
  filter(functional_group != "BRYOPHYTE") %>%
  mutate(Taxon = str_to_sentence(Taxon)) %>%
  filter(!trt == "NPK_Cessation") %>%
  dplyr::select(site_name, ID, Taxon, max_cover) %>%
  filter(Taxon %in% csr_taxa) %>%
  group_by(site_name) %>% 
  nest()

# create the get_special data frame for each species
special_2 <- comm %>%
  mutate(comm_dat = map(data, ~ mutate(get_special_2(.x, CRStable)))) %>%
  select(site_name, comm_dat) %>%
  unnest(cols = c(site_name, comm_dat))

write_csv(special_2, "results/March_2025/strategy_species_level_values.csv")

