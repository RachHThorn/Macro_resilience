# R Thornley
# 15/07/2025
# A2: random effects visualise 
# random effects plots to compare the GLMM results in the context of all the data
# NOTE: this has been simplified for the paper to show only T0-T1 for the Ordered Beta models

library(tidyverse)
library(ggpubr)

################################################################################

# read in effect sizes from the whole of DRAGNet to contextualise the effects from
all <- read_csv("results/July_2025/RE_SE_Taxon_all_DRAGNet.csv")
unique(all$model)
n_distinct(all$group) # 1086 species 

# read in the taxon RE from the compadre and dragnet overlap
taxon <- read_csv("results/July_2025/RE_SE_Taxon_overlap_DRAGNet.csv")
head(taxon)
# get a df of the taxa only
overlap_species <- 
  taxon %>% select(group) %>% 
  arrange(group) %>% unique()
overlap_species <- overlap_species %>% mutate(id = rownames(overlap_species))

# code a column for the overlap species
all <- all %>% left_join(overlap_species) %>%
  mutate(id = case_when(is.na(id) ~ "0",
                        TRUE ~ "1"))

################################################################################

# make sure the labels for the experiments are ordered correctly
unique(all$experiment)
custom_labels = c(DIST = "Disturbance", NPK = "NPK", INTER = "Interaction")


plot <- 
  all %>% filter(model == "Ordbeta") %>%
  filter(time_period == "T0-T1") %>%
  mutate(experiment = factor(experiment, levels = c("DIST", "NPK", "INTER"))) %>%
  ggplot(aes(value, reorder(group, value)))+
  theme_classic()+
  geom_point(size = 1)+
  geom_errorbar(aes(xmin = value - se, xmax = value + se, colour = id), size = 0.3)+
  ylab("Taxon") +
  xlab("Standardised effects with SE")+
  facet_wrap(~experiment, scales = "free_y", labeller = labeller(experiment = custom_labels))+
  geom_vline(aes(xintercept = 0), colour ="red", linetype = "dashed")+
  theme(axis.text.y = element_blank())+
  scale_colour_manual(values = c("black", "red"), labels = c("all DRAGNet", "COMPADRE overlap"))+
  theme(legend.title = element_blank(), axis.ticks.y = element_blank(), legend.position = "top")
plot

#############################################################################

ggsave("figures/Dragnet_draft_July/A2_Caterpillar_context_plots.jpeg", plot, height = 4, width = 7)
