# R Thornley
# 22/07/2025
# visualise BACI concept

# Project: P1_COMPADRE_DRAGNET
# Script: Ax_BACI_example_visual

library(tidyverse)

################################################################################

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

dat_1 <- read_csv("results/July_2025/DRAGNet_overlap_clean_data_with_zeros_T0_T1.csv")
dat_2 <- read_csv("results/July_2025/DRAGNet_overlap_clean_data_with_zeros_T0_T2.csv")
dat_3 <- read_csv("results/July_2025/DRAGNet_overlap_clean_data_with_zeros_T0_T3.csv")

head(dat_1)
head(dat_2)
dat_2$year_trt
dat_2 <- dat_2 %>% filter(year_trt == 2)
dat_3 <- dat_3 %>% filter(year_trt == 3)
all <- dat_1 %>% rbind(dat_2) %>% rbind(dat_3)
unique(all$year_trt)
all %>% count(year_trt)

# filter the whole dataset down to something we can visualise
unique(drag$site_name)
unique(drag$Taxon)
treat_wanted <- c("Control", "Disturbance")
WW <- all %>% 
  filter(site_name == "Wytham Woods") %>%
  filter(trt %in% treat_wanted)
unique(WW$New_taxon) # 17 taxa
names(WW)
WW %>% count(New_taxon, trt, year_trt)

names(WW)

all %>%
  filter(New_taxon == "Agrimonia_eupatoria") %>%
  filter(site_name == "Wytham Woods") %>%
  filter(trt %in% c("Control", "Disturbance")) %>%
  ggplot(aes(year_trt, max_cover, colour = trt, group = trt)) +
  theme_bw()+
  facet_wrap(~block, ncol = 5)+
  geom_smooth(method = "loess")+
  xlab("Year")+
  ylab("relative cover")

ggsave("figures/DRAGNet_explain/BACI_AE_dist.jpeg", height = 3, width = 8)

all %>%
  filter(New_taxon == "Agrimonia_eupatoria") %>%
  filter(site_name == "Wytham Woods") %>%
  filter(trt %in% c("Control", "Disturbance")) %>%
  ggplot(aes(year_trt, max_cover, colour = trt, group = trt)) +
  theme_bw()+
  facet_wrap(~block, ncol = 5)+
  geom_point()+
  geom_line()+
  xlab("Year")+
  ylab("relative cover")
ggsave("figures/DRAGNet_explain/BACI_AE_dist.jpeg", height = 3, width = 8)

all %>%
  filter(New_taxon == "Agrimonia_eupatoria") %>%
  filter(site_name == "Wytham Woods") %>%
  filter(trt %in% c("Control", "NPK")) %>%
  ggplot(aes(year_trt, max_cover, colour = trt, group = trt)) +
  theme_bw()+
  facet_wrap(~block, ncol = 5)+
  geom_point()+
  geom_line()+
  xlab("Year")+
  ylab("relative cover")
ggsave("figures/DRAGNet_explain/BACI_AE_npk.jpeg", height = 3, width = 8)

unique(all$trt)

all %>%
  filter(New_taxon == "Agrimonia_eupatoria") %>%
  filter(site_name == "Wytham Woods") %>%
  filter(trt %in% c("Control", "NPK+Disturbance")) %>%
  ggplot(aes(year_trt, max_cover, colour = trt, group = trt)) +
  theme_bw()+
  facet_wrap(~block, ncol = 5)+
  geom_point()+
  geom_line()+
  xlab("Year")+
  ylab("relative cover")
ggsave("figures/DRAGNet_explain/BACI_AE_npk.jpeg", height = 3, width = 8)




names(all)
table <- all %>% count(site_name, New_taxon, trt)
table
write_csv(table, "results/July_2025/")
