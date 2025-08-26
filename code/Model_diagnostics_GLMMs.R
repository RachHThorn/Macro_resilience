# R Thornley
# 31/07/2025
# Project: P1_COMPADRE_DRAGNET
# Script: Model diagnostics for GLMMs
# calculate R2 marginal (fixed only) and conditional (fixed and random)
# calculate variance components (fixed, random and residual)

# compared the total explained variance between the beta zi hurdle and the ordbeta model
# appendix figure to explain why this model was used 

# use the glmmTMB package
library(tidyverse)
library(easystats)
library(glmmTMB)
library(lme4)
library(mixedup)

# this function gets the variance associated with the different parts of the model
# fixed, random, residual
create_model_summary <- function(model){
  d1 <- r2_nakagawa(model)
  d1a <- (d1[[1]])
  d1b <- (d1[[2]])
  d2 <-icc(model)
  d3 <- get_variance_residual(model)
  d4 <- get_variance_fixed(model)
  d5 <- get_variance_random(model)
  dat <- c(d1a, d1b, d2, d3, d4, d5)
  dat <- as.data.frame(dat)
  dat$var_total <- dat$var.residual + dat$var.fixed + dat$var.random
  dat$prop_fixed <- dat$var.fixed/dat$var_total *100
  dat$prop_random <- dat$var.random/dat$var_total *100
  dat$prop_residual <- dat$var.residual/dat$var_total*100
  return(dat)
}

################################################################################

# Time period 1

dat_1 <- read_csv("results/July_2025/DRAGNet_overlap_clean_data_with_zeros_T0_T1.csv") %>% 
  mutate(trt = factor(trt, levels = c("Control", "Disturbance", "NPK+Disturbance", "NPK", "NPK_Cessation"))) %>%
  mutate(year_trt = factor(year_trt, levels = c("0", "1"))) %>% 
  mutate(site_name = factor(site_name))%>%
  mutate(New_taxon = factor(New_taxon)) %>%
  mutate(taxon_site = paste0(New_taxon, "_", site_name)) %>% 
  mutate(taxon_site = factor(taxon_site)) %>%
  dplyr::select(site_name, New_taxon, taxon_site, year_trt, trt, max_cover, block) %>%
  mutate(new_max_cover = case_when(max_cover == 0 ~ 0,
                                   max_cover > 0 ~ max_cover/100,
                                   TRUE ~ max_cover)) %>%
  mutate(new_max_cover = if_else(new_max_cover == 1, 0.99, new_max_cover))


# DIST
wanted <- c("Disturbance","Control")
DIST <- dat_1 %>% filter(trt %in% wanted)

# build the model and 
Hurdle <- DIST %>%
  glmmTMB(new_max_cover ~ trt * year_trt + (1|New_taxon),
          family = beta_family(link = "logit"), ziformula = ~ 1 + trt * year_trt, data = .) %>%
  create_model_summary() %>%
  mutate(experiment = "DIST") %>%
  mutate(random = "taxa") %>%
  mutate(model = "Hurdle")

Ord <- DIST %>%
  glmmTMB(new_max_cover ~ trt * year_trt + (1|New_taxon), family = ordbeta(link = "logit"), data = .) %>% 
  create_model_summary() %>%
  mutate(experiment = "DIST") %>%
  mutate(random = "taxa") %>%
  mutate(model = "OrdBeta")

DIST_T1 <- rbind(Ord, Hurdle)

# NPK
wanted <- c("NPK","Control")
NPK <- dat_1 %>% filter(trt %in% wanted)


Hurdle <- NPK %>%
  glmmTMB(new_max_cover ~ trt * year_trt + (1|New_taxon),
          family = beta_family(link = "logit"), ziformula = ~ 1 + trt * year_trt, data = .) %>%
  create_model_summary() %>%
  mutate(experiment = "NPK") %>%
  mutate(random = "taxa") %>%
  mutate(model = "Hurdle")

Ord <- dat_1 %>%
  glmmTMB(new_max_cover ~ trt * year_trt + (1|New_taxon), family = ordbeta(link = "logit"), data = .) %>% 
  create_model_summary() %>%
  mutate(experiment = "NPK") %>%
  mutate(random = "taxa") %>%
  mutate(model = "OrdBeta")


NPK_T1 <- rbind(Ord, Hurdle)

# INTER
wanted <- c("NPK+Disturbance", "Control")
INTER <- dat_1 %>% filter(trt %in% wanted)

Hurdle <- INTER %>%
  glmmTMB(new_max_cover ~ trt * year_trt + (1|New_taxon),
          family = beta_family(link = "logit"), ziformula = ~ 1 + trt * year_trt, data = .) %>%
  create_model_summary() %>%
  mutate(experiment = "INTER") %>%
  mutate(random = "taxa") %>%
  mutate(model = "Hurdle")

Ord <- INTER %>%
  glmmTMB(new_max_cover ~ trt * year_trt + (1|New_taxon), family = ordbeta(link = "logit"), data = .) %>% 
  create_model_summary() %>%
  mutate(experiment = "INTER") %>%
  mutate(random = "taxa") %>%
  mutate(model = "OrdBeta")


INTER_T1 <- rbind(Ord, Hurdle)

T1 <- rbind(DIST_T1, NPK_T1, INTER_T1)
T1$time_period <- "T0-T1"

################################################################################

# Time period 2

# read in data and clean up
dat_2 <- read_csv("results/July_2025/DRAGNet_overlap_clean_data_with_zeros_T0_T2.csv") %>% 
  mutate(trt = factor(trt, levels = c("Control", "Disturbance", "NPK+Disturbance", "NPK", "NPK_Cessation"))) %>%
  mutate(year_trt = factor(year_trt, levels = c("0", "2"))) %>% 
  mutate(site_name = factor(site_name))%>%
  mutate(New_taxon = factor(New_taxon)) %>%
  mutate(taxon_site = paste0(New_taxon, "_", site_name)) %>% 
  mutate(taxon_site = factor(taxon_site)) %>%
  dplyr::select(site_name, New_taxon, taxon_site, year_trt, trt, max_cover, block) %>%
  mutate(new_max_cover = case_when(max_cover == 0 ~ 0,
                                   max_cover > 0 ~ max_cover/100,
                                   TRUE ~ max_cover)) %>%
  mutate(new_max_cover = if_else(new_max_cover == 1, 0.99, new_max_cover))

# DIST
wanted <- c("Disturbance","Control")
DIST <- dat_2 %>% filter(trt %in% wanted)

Hurdle <- DIST %>%
  glmmTMB(new_max_cover ~ trt * year_trt + (1|New_taxon),
          family = beta_family(link = "logit"), ziformula = ~ 1 + trt * year_trt, data = .) %>%
  create_model_summary() %>%
  mutate(experiment = "DIST") %>%
  mutate(random = "taxa") %>%
  mutate(model = "Hurdle")

Ord <- DIST %>%
  glmmTMB(new_max_cover ~ trt * year_trt + (1|New_taxon), family = ordbeta(link = "logit"), data = .) %>% 
  create_model_summary() %>%
  mutate(experiment = "DIST") %>%
  mutate(random = "taxa") %>%
  mutate(model = "OrdBeta")


DIST_T2 <- rbind(Ord, Hurdle)

# NPK
wanted <- c("NPK","Control")
NPK <- dat_2 %>% filter(trt %in% wanted)


Hurdle <- NPK %>%
  glmmTMB(new_max_cover ~ trt * year_trt + (1|New_taxon),
          family = beta_family(link = "logit"), ziformula = ~ 1 + trt * year_trt, data = .) %>%
  create_model_summary() %>%
  mutate(experiment = "NPK") %>%
  mutate(random = "taxa") %>%
  mutate(model = "Hurdle")

Ord <- NPK %>%
  glmmTMB(new_max_cover ~ trt * year_trt + (1|New_taxon), family = ordbeta(link = "logit"), data = .) %>% 
  create_model_summary() %>%
  mutate(experiment = "NPK") %>%
  mutate(random = "taxa") %>%
  mutate(model = "OrdBeta")

NPK_T2 <- rbind(Ord, Hurdle)

# INTER
wanted <- c("NPK+Disturbance", "Control")
INTER <- dat_2 %>% filter(trt %in% wanted)

Hurdle <- INTER %>%
  glmmTMB(new_max_cover ~ trt * year_trt + (1|New_taxon),
          family = beta_family(link = "logit"), ziformula = ~ 1 + trt * year_trt, data = .) %>%
  create_model_summary() %>%
  mutate(experiment = "INTER") %>%
  mutate(random = "taxa") %>%
  mutate(model = "Hurdle")

Ord <- INTER %>%
  glmmTMB(new_max_cover ~ trt * year_trt + (1|New_taxon), family = ordbeta(link = "logit"), data = .) %>% 
  create_model_summary() %>%
  mutate(experiment = "INTER") %>%
  mutate(random = "taxa") %>%
  mutate(model = "OrdBeta")


INTER_T2 <- rbind(Ord, Hurdle)

T2 <- rbind(DIST_T2, NPK_T2, INTER_T2)
T2$time_period <- "T0-T2"

###############################################################################

# Time period 3

# read in data and clean up
dat_3 <- read_csv("results/July_2025/DRAGNet_overlap_clean_data_with_zeros_T0_T3.csv") %>% 
  mutate(trt = factor(trt, levels = c("Control", "Disturbance", "NPK+Disturbance", "NPK", "NPK_Cessation"))) %>%
  mutate(year_trt = factor(year_trt, levels = c("0", "3"))) %>% 
  mutate(site_name = factor(site_name))%>%
  mutate(New_taxon = factor(New_taxon)) %>%
  mutate(taxon_site = paste0(New_taxon, "_", site_name)) %>% 
  mutate(taxon_site = factor(taxon_site)) %>%
  dplyr::select(site_name, New_taxon, taxon_site, year_trt, trt, max_cover, block) %>%
  mutate(new_max_cover = case_when(max_cover == 0 ~ 0,
                                   max_cover > 0 ~ max_cover/100,
                                   TRUE ~ max_cover)) %>%
  mutate(new_max_cover = if_else(new_max_cover == 1, 0.99, new_max_cover))

# DIST
wanted <- c("Disturbance","Control")
DIST <- dat_3 %>% filter(trt %in% wanted)

Hurdle <- DIST %>%
  glmmTMB(new_max_cover ~ trt * year_trt + (1|New_taxon),
          family = beta_family(link = "logit"), ziformula = ~ 1 + trt * year_trt, data = .) %>%
  create_model_summary() %>%
  mutate(experiment = "DIST") %>%
  mutate(random = "taxa") %>%
  mutate(model = "Hurdle")

Ord <- DIST %>%
  glmmTMB(new_max_cover ~ trt * year_trt + (1|New_taxon), family = ordbeta(link = "logit"), data = .) %>% 
  create_model_summary() %>%
  mutate(experiment = "DIST") %>%
  mutate(random = "taxa") %>%
  mutate(model = "OrdBeta")


DIST_T3 <- rbind(Ord, Hurdle)

# NPK
wanted <- c("NPK","Control")
NPK <- dat_3 %>% filter(trt %in% wanted)


Hurdle <- NPK %>%
  glmmTMB(new_max_cover ~ trt * year_trt + (1|New_taxon),
          family = beta_family(link = "logit"), ziformula = ~ 1 + trt * year_trt, data = .) %>%
  create_model_summary() %>%
  mutate(experiment = "NPK") %>%
  mutate(random = "taxa") %>%
  mutate(model = "Hurdle")

Ord <- NPK %>%
  glmmTMB(new_max_cover ~ trt * year_trt + (1|New_taxon), family = ordbeta(link = "logit"), data = .) %>% 
  create_model_summary() %>%
  mutate(experiment = "NPK") %>%
  mutate(random = "taxa") %>%
  mutate(model = "OrdBeta")

NPK_T3 <- rbind(Ord, Hurdle)
dat_3$trt
# INTER
wanted <- c("NPK+Disturbance", "Control")
INTER <- dat_3 %>% filter(trt %in% wanted)

Hurdle <- INTER %>%
  glmmTMB(new_max_cover ~ trt * year_trt + (1|New_taxon),
          family = beta_family(link = "logit"), ziformula = ~ 1 + trt * year_trt, data = .) %>%
  create_model_summary() %>%
  mutate(experiment = "INTER") %>%
  mutate(random = "taxa") %>%
  mutate(model = "Hurdle")

Ord <- INTER %>%
  glmmTMB(new_max_cover ~ trt * year_trt + (1|New_taxon), family = ordbeta(link = "logit"), data = .) %>% 
  create_model_summary() %>%
  mutate(experiment = "INTER") %>%
  mutate(random = "taxa") %>%
  mutate(model = "OrdBeta")


INTER_T3 <- rbind(Ord, Hurdle)

T3 <- rbind(DIST_T3, NPK_T3, INTER_T3)
T3$time_period <- "T0-T3"

################################################################################

all <- rbind(T1, T2, T3)
write_csv(all, "results/July_2025/GLMMs_ICC_R2_varcomp_summary.csv")

################################################################################

# visualise this as a barchart

all$experiment
all$model
plot_1 <-
  all %>%
  dplyr::select(prop_fixed:time_period) %>%
  pivot_longer(cols = prop_fixed:prop_residual, names_to = "part", values_to = "value") %>%
  ggplot(aes(x = experiment, y = value, fill= part, label = value)) +
  theme_classic()+
  theme(legend.title = element_blank())+
  geom_bar(position = "stack", stat = "identity") +
  facet_grid(model~time_period)+
  scale_fill_manual(values = c("black", "darkolivegreen2", "grey80"), 
                    labels = c("Fixed", "Random", "Residual"))+
  ylab("Variance Explained")+
  xlab("Experiment")+
  theme(legend.position = "bottom")+
  theme(text = element_text(size = 14))+
  theme(legend.text=element_text(size=rel(1)))+
  ggtitle("Hurdle Beta and Ordered beta models -\nVariance components when random intercept = taxa")
plot_1

ggsave("figures/July_2025/variance_components_hurdle_ordBeta_models.jpeg", plot_1, height = 6, width = 6)

################################################################################
all


plot_2 <- 
  all %>% ggplot(aes(model, ICC_adjusted)) +
  theme_bw()+
  geom_col() + 
  facet_grid(experiment~time_period)
plot_2  
  