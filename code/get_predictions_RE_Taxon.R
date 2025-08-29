# what happens if we try and predict for each level of the RE


# get predictions for each level of a RE model

# load packages
packages <- c("tidyverse", "glmmTMB", "rstan", "assertthat", "purrr", "remotes", 
              "mixedup", "easystats", "purrr")

for (p in packages) {
  if (!requireNamespace(p, quietly = TRUE)) {
    install.packages(p)
  }
  library(p, character.only = TRUE)
}

################################################################################

# test model code for a single experiment and time period

data <- read_csv("results/DRAGNet_T0_T1_all.csv")
experiments <- c("Control", "Disturbance")
data <- data %>% filter(trt %in% experiments)
names(data)
data <- data %>% mutate(site_block = paste0(site_name, "_", block))
data <- data %>% mutate(site_block_taxon = paste0(site_name, "_", block, "_", New_taxon))
n_distinct(data$site_block) # there are 216 unique combinations of site and block
n_distinct(data$site_block_taxon) # there are 4002 unique combinations of site, block, taxon

data %>% group_by(site_block) %>% count()
data %>% group_by(site_block_taxon) %>% count() # these data sets look very small to converge to a RE
# but they specify the paired nature of this data

################################################################################

# run the betreg model
betareg_mod <- glmmTMB(
  new_max_cover ~ trt * year_trt + (1|New_taxon) + + (1|site_block_taxon),
  family = ordbeta(link = "logit"),
  data = data)
  
# not try and predict from it for each level of the RE (New_Taxon)
# these are called Unit-level predictions
# use the ggeffects model to do this
me <- predict_response(betareg_mod, terms = c("New_taxon"), type = "random")
plot(me, show_ci = TRUE) # rough plot


################################################################################

# plot these proportions as a caterpillar plot

################################################################################

# read in the taxon RE from the compadre and dragnet overlap
comp_taxa <- read_csv("results/common_species_drag_comp.csv") %>% pull(x)

# select the random effects for the var of interest only
# create a new ID variable
# change some labels in the data frame to tidy it up
names(me)
me <- 
  me %>% 
  as.data.frame() %>%
  dplyr::select(x, predicted, std.error, conf.low, conf.high) %>%
  mutate(ID = if_else(x %in% comp_taxa, "1", "0")) %>%
  mutate(ID = factor(ID, levels = c("0", "1"))) %>%
  rename(taxon = x) %>%
  mutate(label = if_else(ID == "1", taxon, NA_character_))

# plot all the effects with str error
plot <- 
  me %>% 
  ggplot(aes(predicted, reorder(taxon, predicted)))+
  theme_classic()+
  geom_point(size = 1)+
  geom_errorbar(aes(xmin = predicted - std.error, xmax = predicted + std.error, colour = ID), size = 0.3)+
  geom_text_repel(
    aes(label = label), size = 3, 
    nudge_x = -0.1,
    direction = "y",
    hjust = 1,
    na.rm = TRUE) +  # labels only where ID==1
  ylab("Taxon") +
  xlab("Standardised effects with SE")+
  geom_vline(aes(xintercept = 0), colour ="red", linetype = "dashed")+
  theme(axis.text.y = element_blank())+
  scale_colour_manual(values = c("black", "red"), labels = c("all DRAGNet", "COMPADRE overlap"))+
  theme(legend.title = element_blank(), axis.ticks.y = element_blank(), legend.position = "top")
plot

names(me)
# plot all the effects with CIs
plot <- 
  me %>% 
  ggplot(aes(predicted, reorder(taxon, predicted)))+
  theme_classic()+
  geom_point(size = 1)+
  geom_errorbar(aes(xmin = conf.low, xmax = conf.high, colour = ID), size = 0.3)+
  geom_text_repel(
    aes(label = label), size = 3, 
    nudge_x = -0.1,
    direction = "y",
    hjust = 1,
    na.rm = TRUE) +  # labels only where ID==1
  ylab("Taxon") +
  xlab("Standardised effects with SE")+
  geom_vline(aes(xintercept = 0), colour ="red", linetype = "dashed")+
  theme(axis.text.y = element_blank())+
  scale_colour_manual(values = c("black", "red"), labels = c("all DRAGNet", "COMPADRE overlap"))+
  theme(legend.title = element_blank(), axis.ticks.y = element_blank(), legend.position = "top")
plot

# just show the compadre and dragnet overlap
plot <- 
  me %>% 
  filter(ID == "1") %>%
  ggplot(aes(predicted, reorder(taxon, predicted)))+
  theme_classic()+
  geom_point(size = 1)+
  geom_errorbar(aes(xmin = predicted - std.error, xmax = predicted + std.error, colour = ID), size = 0.3)+
  ylab("Taxon") +
  xlab("Standardised effects with SE")+
  geom_vline(aes(xintercept = 0), colour ="red", linetype = "dashed")+
  theme(axis.text.y = element_blank())+
  scale_colour_manual(values = c("black", "red"), labels = c("all DRAGNet", "COMPADRE overlap"))+
  theme(legend.title = element_blank(), axis.ticks.y = element_blank(), legend.position = "top")
plot 

# how about the confidence intervals
plot <- 
  me %>% 
  filter(ID == "1") %>%
  ggplot(aes(predicted, reorder(taxon, predicted)))+
  theme_classic()+
  geom_point(size = 1)+
  geom_errorbar(aes(xmin = conf.low, xmax = conf.high, colour = ID), size = 0.3)+
  ylab("Taxon") +
  xlab("Standardised effects with CI")+
  geom_vline(aes(xintercept = 0), colour ="red", linetype = "dashed")+
  theme(axis.text.y = element_blank())+
  scale_colour_manual(values = c("black", "red"), labels = c("all DRAGNet", "COMPADRE overlap"))+
  theme(legend.title = element_blank(), axis.ticks.y = element_blank(), legend.position = "top")
plot

