# R Thornley
# 30/08/2025
# Regressions for the mean Life History Traits and mean RE outcomes

library(tidyverse)
library(ggpubr)

# Load demography variables and find the mean and SE of each Taxon and trait
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
  mutate(Demo_value_log = log(Demo_value))
names(demo)
demo$Taxon_label
demo$Demo_value_log
# list the demography traits we are dealng with
print(unique(demo$Demo_trait))
# "T_generation"    
# "R0"              
# "percapita_repro" 
# "Lmean"           
# "age_repro"       
# "Lambda"         
# "Reactivity"      
# "FirstStepAtt"    
# "Lmax"  

################################################################################

# make plots for each demo trait
# label x axis with the sample size and the Taxon

# "T_generation"   
p1 <-  
  demo %>%
  filter(Demo_trait == "T_generation") %>%
  drop_na(Demo_value) %>%
  ggplot(aes(reorder(Taxon_label, Demo_value_mean), Demo_value))+
  geom_boxplot()+
  theme_classic()+
  theme(axis.text.x = element_text(angle = 90, hjust = 1),
        axis.title.x = element_blank()) +
  ylab("Generation Time")
p1
# ggsave("figures/T_generation_boxplot.jpeg", p1, height = 5, width = 6)


# "R0" 
p2 <-  
  demo %>%
  filter(Demo_trait == "R0") %>%
  drop_na(Demo_value) %>%
  ggplot(aes(reorder(Taxon_label, Demo_value_mean), Demo_value))+
  geom_boxplot()+
  theme_classic()+
  theme(axis.text.x = element_text(angle = 90, hjust = 1),
        axis.title.x = element_blank()) +
  ylab("Reproductive output")
p2
# ggsave("figures/Reproductive_rate_boxplot.jpeg", p2, height = 5, width = 6)

# "percapita_repro" 
p3 <-  
  demo %>%
  filter(Demo_trait == "percapita_repro") %>%
  drop_na(Demo_value) %>%
  ggplot(aes(reorder(Taxon_label, Demo_value_mean), Demo_value))+
  geom_boxplot()+
  theme_classic()+
  theme(axis.text.x = element_text(angle = 90, hjust = 1),
        axis.title.x = element_blank()) +
  ylab("Per capita reproduction")
p3
# ggsave("figures/Percapita_reproduction_boxplot.jpeg", p3, height = 5, width = 6)

# Lmean
p4 <-  
  demo %>%
  filter(Demo_trait == "Lmean") %>%
  drop_na(Demo_value) %>%
  ggplot(aes(reorder(Taxon_label, Demo_value_mean), Demo_value))+
  geom_boxplot()+
  theme_classic()+
  theme(axis.text.x = element_text(angle = 90, hjust = 1),
        axis.title.x = element_blank()) +
  ylab("Mean Life Expectancy")
p4
# ggsave("figures/Mean_life_expectancy_boxplot.jpeg", p4, height = 5, width = 6)

# "age_repro"
p5 <-  
  demo %>%
  filter(Demo_trait == "R0") %>%
  drop_na(Demo_value) %>%
  ggplot(aes(reorder(Taxon_label, Demo_value_mean), Demo_value))+
  geom_boxplot()+
  theme_classic()+
  theme(axis.text.x = element_text(angle = 90, hjust = 1),
        axis.title.x = element_blank()) +
  ylab("Reproductive output")
p5
# ggsave("figures/Reproductive_output_boxplot.jpeg", p5, height = 5, width = 6)

# "Lambda"   
p6 <- 
  demo %>%
  filter(Demo_trait == "Lambda") %>%
  drop_na(Demo_value) %>%
  ggplot(aes(reorder(Taxon_label, Demo_value_mean), Demo_value))+
  geom_boxplot()+
  theme_classic()+
  theme(axis.text.x = element_text(angle = 90, hjust = 1),
        axis.title.x = element_blank()) +
  ylab("Lambda")
p6
# ggsave("figures/Lambda_boxplot.jpeg", p6, height = 5, width = 6)

# "Reactivity"    
p7 <- 
  demo %>%
  filter(Demo_trait == "Reactivity") %>%
  drop_na(Demo_value) %>%
  ggplot(aes(reorder(Taxon_label, Demo_value_mean), Demo_value))+
  geom_boxplot()+
  theme_classic()+
  theme(axis.text.x = element_text(angle = 90, hjust = 1),
        axis.title.x = element_blank()) +
  ylab("Reactivity")
p7
# ggsave("figures/Reactivity_boxplot.jpeg", p7, height = 5, width = 6)

#  "FirstStepAtt"   
p8 <- 
  demo %>%
  filter(Demo_trait == "FirstStepAtt") %>%
  drop_na(Demo_value) %>%
  ggplot(aes(reorder(Taxon_label, Demo_value_mean), Demo_value))+
  geom_boxplot()+
  theme_classic()+
  theme(axis.text.x = element_text(angle = 90, hjust = 1),
        axis.title.x = element_blank()) +
  ylab("FirstStepAtt")
p8
# ggsave("figures/FirstStepAtt_boxplot.jpeg", p8, height = 5, width = 6)

# "Lmax" 
p9 <- 
  demo %>%
  filter(Demo_trait == "Lmax") %>%
  drop_na(Demo_value) %>%
  ggplot(aes(reorder(Taxon_label, Demo_value_mean), Demo_value))+
  geom_boxplot()+
  theme_classic()+
  theme(axis.text.x = element_text(angle = 90, hjust = 1),
        axis.title.x = element_blank()) +
  ylab("Maximum Life Expectancy")
p9
# ggsave("figures/Maximum_life_expectancy_boxplot.jpeg", p9, height = 5, width = 6)

################################################################################

combined <- ggarrange(p1, p2, p3, p4, p5, p6, p7, p8, p9,
                      ncol = 3, nrow = 3,
                      common.legend = FALSE)

combined <- annotate_figure(
  combined,
  bottom = text_grob("Taxon", color = "black", face = "bold", size = 14),
  left = text_grob("Demographic Metric", color = "black", face = "bold", size = 14, rot = 90)
)
combined
ggsave("figures/all_demo_metrics_boxplots.jpeg", combined, height = 18, width = 12)


################################################################################

# Try some of the values logged as they have some very large values

# "R0" 
p2_again <-  
  demo %>%
  filter(Demo_trait == "R0") %>%
  drop_na(Demo_value) %>%
  ggplot(aes(reorder(Taxon_label, Demo_value_mean), Demo_value_log))+
  geom_boxplot()+
  theme_classic()+
  theme(axis.text.x = element_text(angle = 90, hjust = 1),
        axis.title.x = element_blank()) +
  ylab("Reproductive output (log)")
p2_again

# "percapita repro"
p3_again <-  
  demo %>%
  filter(Demo_trait == "percapita_repro") %>%
  drop_na(Demo_value) %>%
  ggplot(aes(reorder(Taxon_label, Demo_value_mean), Demo_value_log))+
  geom_boxplot()+
  theme_classic()+
  theme(axis.text.x = element_text(angle = 90, hjust = 1),
        axis.title.x = element_blank()) +
  ylab("Per capita reproduction")
p3_again

# "Lmean"
p4_again <-  
  demo %>%
  filter(Demo_trait == "Lmean") %>%
  drop_na(Demo_value) %>%
  ggplot(aes(reorder(Taxon_label, Demo_value_mean), Demo_value_log))+
  geom_boxplot()+
  theme_classic()+
  theme(axis.text.x = element_text(angle = 90, hjust = 1),
        axis.title.x = element_blank()) +
  ylab("Mean Life Expectancy (log)")
p4

# "Reactivity"    
p7_again <- 
  demo %>%
  filter(Demo_trait == "Reactivity") %>%
  drop_na(Demo_value) %>%
  ggplot(aes(reorder(Taxon_label, Demo_value_mean), Demo_value_log))+
  geom_boxplot()+
  theme_classic()+
  theme(axis.text.x = element_text(angle = 90, hjust = 1),
        axis.title.x = element_blank()) +
  ylab("Reactivity (log)")
p7_again

# "Lmax"
p9_again <- 
  demo %>%
  filter(Demo_trait == "Lmax") %>%
  drop_na(Demo_value) %>%
  ggplot(aes(reorder(Taxon_label, Demo_value_mean), Demo_value_log))+
  geom_boxplot()+
  theme_classic()+
  theme(axis.text.x = element_text(angle = 90, hjust = 1),
        axis.title.x = element_blank()) +
  ylab("Maximum Life Expectancy (log)")
p9_again

################################################################################


demo %>%
  filter(Demo_trait == "age_repro") %>%
  ggplot(aes(Demo_value, reorder(Taxon, Demo_value)))+
  theme_classic()+
  geom_point(size = 1)+
  geom_errorbar(aes(xmin = Demo_value - se, xmax = Demo_value + se), size = 0.3)+
  facet_wrap(~ Demo_trait, scales = "free_x")


demo %>%
  ggplot(aes(Demo_value, reorder(Taxon, Demo_value)))+
  theme_classic()+
  geom_point(size = 1)+
  geom_errorbar(aes(xmin = Demo_value - se, xmax = Demo_value + se), size = 0.3)+
  facet_wrap(~ Demo_trait, scales = "free_x")+
  xlab("Years")


