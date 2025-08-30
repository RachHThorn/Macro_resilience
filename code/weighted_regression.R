# R Thornley
# 29/08/2025
# Regressions for the mean Life History Traits and mean RE outcomes

library(tidyverse)


# join the effect sizes 

# Load demography variables and find the mean and SE of each Taxon and trait

demo <- read_csv("results/all_COMPADRE_metrics.csv") %>% 
  filter(DRAGNet == TRUE) %>% 
  rename("Taxon" = "SpeciesAccepted", "Demo_trait" = "demo_var", "Demo_value" = "value") %>%
  mutate(Taxon = str_replace_all(Taxon, " ", "_")) %>%
  select(Taxon, Demo_trait, Demo_value) %>% 
  group_by(Taxon, Demo_trait) %>%
  summarise(Sample_size = n(),
            s = sd(Demo_value, na.rm = TRUE),
            se = s/sqrt(Sample_size),
            Demo_value = mean(Demo_value, na.rm = TRUE))

names(demo)
demo$Demo_trait

demo %>%
  filter(Demo_trait == "age_repro") %>%
  ggplot(aes(Demo_value, reorder(Taxon, Demo_value)))+
  theme_classic()+
  geom_point(size = 1)+
  geom_errorbar(aes(xmin = Demo_value - se, xmax = Demo_value + se), size = 0.3)+
  facet_wrap(~ Demo_trait, scales = "free_x")

demo %>%
  filter(Demo_trait == "T_generation") %>%
  ggplot(aes(Demo_value, reorder(Taxon, Demo_value)))+
  theme_classic()+
  geom_point(size = 1)+
  geom_errorbar(aes(xmin = Demo_value - se, xmax = Demo_value + se), size = 0.3)+
  facet_wrap(~ Demo_trait, scales = "free_x")


            
# drop_na(Demo_value)

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