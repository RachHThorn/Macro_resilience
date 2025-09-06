# R Thornley
# 22/07/2025
# Interrogate demographic metrics 

library(tidyverse)

# look at the spread of values of the demographic metrics
demo <- read_csv("results/July_2025/all_compadre_metrics.csv")
head(demo)
nos_matrices <-
  demo %>% filter(DRAGNet == TRUE) %>% 
  count(SpeciesAccepted, demo_var) %>% 
  filter(demo_var == "Lmax")

drag <- demo %>% filter(DRAGNet == TRUE) 
n_distinct(drag$MatrixID) # 334 matrices

ggplot(drag, aes(value))+
  geom_density()+
  facet_wrap(~)



