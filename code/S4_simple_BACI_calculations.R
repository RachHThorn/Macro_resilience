# R Thornley
# 16/07/2025
# Project: P1_COMPADRE_DRAGNET
# Script: S4_simple_BACI_calculations

# calculate simple BACI equation for each of the experiments 
# for both the time points

# NOTE: This needs tidying as the code is dirty at the mo

################################################################################

# T0-T1 

# read in all DRAGNet clean data
T1 <- read_csv("results/July_2025/DRAGNet_all_clean_data_with_zeros_T0_T1.csv")
names(T1)
unique(T1$trt)

# the simple BACIs are row wise calculations so we need to convert the data
wanted <- c("Control", "Disturbance")
DIST <- T1 %>% 
  mutate(site_taxon = paste0(site_name, sep = "_", New_taxon, sep = "_", year_trt)) %>%
  filter(trt %in% wanted) %>%
  group_by(site_name, New_taxon, year_trt, trt, block) %>% 
  summarise(mean_cover = mean(new_max_cover)) %>%
  mutate(trt_time = paste0(trt, sep = "_", year_trt)) %>%
  ungroup() %>%
  dplyr::select(site_name, New_taxon, block, trt_time, mean_cover) %>%
  group_by(site_name, New_taxon, block) %>%
  pivot_wider(values_from = mean_cover, names_from = trt_time)

# rename the columns to easier names
names(DIST) <- c("site_name", "Taxon", "block", "Cont_0", "Dist_0",  "Cont_1", "Dist_1")
# some rows now have NA values in them as we have done pivot wider
# replace with 0
colSums(is.na(DIST)) # there are a lot of NAs 
DIST <- DIST %>%
  mutate(across(everything(), ~ replace_na(., 0)))

# also get rid of rows where all the inputs are zeros them as they add nothing to the analysis
DIST <- DIST %>%
  rowwise() %>%
  filter(any(across(where(is.numeric), ~ . != 0))) %>%
  ungroup()
colSums(DIST == 0, na.rm = TRUE) # counts the number of zeros in each of the columns

# calculate the BACI at the quadrat scale
Quadrat_T1_DIST <- 
  DIST %>% 
  group_by(Taxon, site_name, block) %>% 
  mutate(t1_dist = Dist_1 - Cont_1) %>%
  mutate(t0_dist = Dist_0 - Cont_0) %>%
  mutate(BACI_dist = t1_dist-t0_dist) %>%
  select(Taxon, site_name, block, BACI_dist) %>%
  mutate(group_var = "quadrat") %>%
  mutate(time = "T0-T1") %>%
  mutate(trt = "DIST") %>%
  rename(BACI = BACI_dist)

n_distinct(Quadrat_T1_DIST$Taxon) # 716

# the simple BACIs are row wise calculations so we need to convert the data
wanted <- c("Control", "NPK")
NPK <- T1 %>% 
  mutate(site_taxon = paste0(site_name, sep = "_", New_taxon, sep = "_", year_trt)) %>%
  filter(trt %in% wanted) %>%
  group_by(site_name, New_taxon, year_trt, trt, block) %>% 
  summarise(mean_cover = mean(new_max_cover)) %>%
  mutate(trt_time = paste0(trt, sep = "_", year_trt)) %>%
  ungroup() %>%
  dplyr::select(site_name, New_taxon, block, trt_time, mean_cover) %>%
  group_by(site_name, New_taxon, block) %>%
  pivot_wider(values_from = mean_cover, names_from = trt_time)

# rename the columns to easier names
names(NPK) <- c("site_name", "Taxon", "block", "Cont_0", "NPK_0",  "Cont_1", "NPK_1")
# some rows now have NA values in them as we have done pivot wider
# replace with 0
colSums(is.na(NPK)) # there are a lot of NAs 
NPK <- NPK %>%
  mutate(across(everything(), ~ replace_na(., 0)))

# also get rid of rows where all the inputs are zeros them as they add nothing to the analysis
NPK <- NPK %>%
  rowwise() %>%
  filter(any(across(where(is.numeric), ~ . != 0))) %>%
  ungroup()
colSums(NPK == 0, na.rm = TRUE) # counts the number of zeros in each of the columns

# calculate the BACI at the quadrat scale
Quadrat_T1_NPK <- 
  NPK %>% 
  group_by(Taxon, site_name, block) %>% 
  mutate(t1_npk = NPK_1 - Cont_1) %>%
  mutate(t0_npk = NPK_0 - Cont_0) %>%
  mutate(BACI_npk = t1_npk - t0_npk) %>%
  select(Taxon, site_name, block, BACI_npk) %>%
  mutate(group_var = "quadrat") %>%
  mutate(time = "T0-T1") %>%
  mutate(trt = "NPK") %>%
  rename(BACI = BACI_npk) 

n_distinct(Quadrat_T1_NPK$Taxon) # 679

# the simple BACIs are row wise calculations so we need to convert the data
wanted <- c("Control", "NPK+Disturbance")
INTER <- T1 %>% 
  mutate(site_taxon = paste0(site_name, sep = "_", New_taxon, sep = "_", year_trt)) %>%
  filter(trt %in% wanted) %>%
  group_by(site_name, New_taxon, year_trt, trt, block) %>% 
  summarise(mean_cover = mean(new_max_cover)) %>%
  mutate(trt_time = paste0(trt, sep = "_", year_trt)) %>%
  ungroup() %>%
  dplyr::select(site_name, New_taxon, block, trt_time, mean_cover) %>%
  group_by(site_name, New_taxon, block) %>%
  pivot_wider(values_from = mean_cover, names_from = trt_time)

# rename the columns to easier names
names(INTER) <- c("site_name", "Taxon", "block", "Cont_0", "INTER_0",  "Cont_1", "INTER_1")
# some rows now have NA values in them as we have done pivot wider
# replace with 0
colSums(is.na(INTER)) # there are a lot of NAs 
INTER <- INTER %>%
  mutate(across(everything(), ~ replace_na(., 0)))

# also get rid of rows where all the inputs are zeros them as they add nothing to the analysis
INTER <- INTER %>%
  rowwise() %>%
  filter(any(across(where(is.numeric), ~ . != 0))) %>%
  ungroup()
colSums(INTER == 0, na.rm = TRUE) # counts the number of zeros in each of the columns

# calculate the BACI at the quadrat scale
Quadrat_T1_INTER <- 
  INTER %>% 
  group_by(Taxon, site_name, block) %>% 
  mutate(t1_inter = INTER_1 - Cont_1) %>%
  mutate(t0_inter = INTER_0 - Cont_0) %>%
  mutate(BACI_inter = t1_inter - t0_inter) %>%
  select(Taxon, site_name, block, BACI_inter) %>%
  mutate(group_var = "quadrat") %>%
  mutate(time = "T0-T1") %>%
  mutate(trt = "INTER") %>%
  rename(BACI = BACI_inter) 

n_distinct(Quadrat_T1_INTER$Taxon) # 711

T1_result <- rbind(Quadrat_T1_DIST, Quadrat_T1_NPK, Quadrat_T1_INTER)
write_csv(T1_result, "results/July_2025/BACI_T1_results.csv")

################################################################################

# T0-T2

# read in all DRAGNet clean data
T2 <- read_csv("results/July_2025/DRAGNet_all_clean_data_with_zeros_T0_T2.csv")
names(T2)
unique(T2$trt)

# the simple BACIs are row wise calculations so we need to convert the data
wanted <- c("Control", "Disturbance")
DIST <- T2 %>% 
  mutate(site_taxon = paste0(site_name, sep = "_", New_taxon, sep = "_", year_trt)) %>%
  filter(trt %in% wanted) %>%
  group_by(site_name, New_taxon, year_trt, trt, block) %>% 
  summarise(mean_cover = mean(new_max_cover)) %>%
  mutate(trt_time = paste0(trt, sep = "_", year_trt)) %>%
  ungroup() %>%
  dplyr::select(site_name, New_taxon, block, trt_time, mean_cover) %>%
  group_by(site_name, New_taxon, block) %>%
  pivot_wider(values_from = mean_cover, names_from = trt_time)

# rename the columns to easier names
names(DIST) <- c("site_name", "Taxon", "block", "Cont_0", "Dist_0",  "Cont_1", "Dist_1")
# some rows now have NA values in them as we have done pivot wider
# replace with 0
colSums(is.na(DIST)) # there are a lot of NAs 
DIST <- DIST %>%
  mutate(across(everything(), ~ replace_na(., 0)))

# also get rid of rows where all the inputs are zeros them as they add nothing to the analysis
DIST <- DIST %>%
  rowwise() %>%
  filter(any(across(where(is.numeric), ~ . != 0))) %>%
  ungroup()
colSums(DIST == 0, na.rm = TRUE) # counts the number of zeros in each of the columns

# calculate the BACI at the quadrat scale
Quadrat_T2_DIST <- 
  DIST %>% 
  group_by(Taxon, site_name, block) %>% 
  mutate(t1_dist = Dist_1 - Cont_1) %>%
  mutate(t0_dist = Dist_0 - Cont_0) %>%
  mutate(BACI_dist = t1_dist-t0_dist) %>%
  select(Taxon, site_name, block, BACI_dist) %>%
  mutate(group_var = "quadrat") %>%
  mutate(time = "T0-T2") %>%
  mutate(trt = "DIST") %>%
  rename(BACI = BACI_dist)

n_distinct(Quadrat_T2_DIST$Taxon) # 627

# the simple BACIs are row wise calculations so we need to convert the data
wanted <- c("Control", "NPK")
NPK <- T2 %>% 
  mutate(site_taxon = paste0(site_name, sep = "_", New_taxon, sep = "_", year_trt)) %>%
  filter(trt %in% wanted) %>%
  group_by(site_name, New_taxon, year_trt, trt, block) %>% 
  summarise(mean_cover = mean(new_max_cover)) %>%
  mutate(trt_time = paste0(trt, sep = "_", year_trt)) %>%
  ungroup() %>%
  dplyr::select(site_name, New_taxon, block, trt_time, mean_cover) %>%
  group_by(site_name, New_taxon, block) %>%
  pivot_wider(values_from = mean_cover, names_from = trt_time)

# rename the columns to easier names
names(NPK) <- c("site_name", "Taxon", "block", "Cont_0", "NPK_0",  "Cont_1", "NPK_1")
# some rows now have NA values in them as we have done pivot wider
# replace with 0
colSums(is.na(NPK)) # there are a lot of NAs 
NPK <- NPK %>%
  mutate(across(everything(), ~ replace_na(., 0)))

# also get rid of rows where all the inputs are zeros them as they add nothing to the analysis
NPK <- NPK %>%
  rowwise() %>%
  filter(any(across(where(is.numeric), ~ . != 0))) %>%
  ungroup()
colSums(NPK == 0, na.rm = TRUE) # counts the number of zeros in each of the columns

# calculate the BACI at the quadrat scale
Quadrat_T2_NPK <- 
  NPK %>% 
  group_by(Taxon, site_name, block) %>% 
  mutate(t1_npk = NPK_1 - Cont_1) %>%
  mutate(t0_npk = NPK_0 - Cont_0) %>%
  mutate(BACI_npk = t1_npk - t0_npk) %>%
  select(Taxon, site_name, block, BACI_npk) %>%
  mutate(group_var = "quadrat") %>%
  mutate(time = "T0-T2") %>%
  mutate(trt = "NPK") %>%
  rename(BACI = BACI_npk) 

n_distinct(Quadrat_T2_NPK$Taxon) # 584

# the simple BACIs are row wise calculations so we need to convert the data
wanted <- c("Control", "NPK+Disturbance")
INTER <- T2 %>% 
  mutate(site_taxon = paste0(site_name, sep = "_", New_taxon, sep = "_", year_trt)) %>%
  filter(trt %in% wanted) %>%
  group_by(site_name, New_taxon, year_trt, trt, block) %>% 
  summarise(mean_cover = mean(new_max_cover)) %>%
  mutate(trt_time = paste0(trt, sep = "_", year_trt)) %>%
  ungroup() %>%
  dplyr::select(site_name, New_taxon, block, trt_time, mean_cover) %>%
  group_by(site_name, New_taxon, block) %>%
  pivot_wider(values_from = mean_cover, names_from = trt_time)

# rename the columns to easier names
names(INTER) <- c("site_name", "Taxon", "block", "Cont_0", "INTER_0",  "Cont_1", "INTER_1")
# some rows now have NA values in them as we have done pivot wider
# replace with 0
colSums(is.na(INTER)) # there are a lot of NAs 
INTER <- INTER %>%
  mutate(across(everything(), ~ replace_na(., 0)))

# also get rid of rows where all the inputs are zeros them as they add nothing to the analysis
INTER <- INTER %>%
  rowwise() %>%
  filter(any(across(where(is.numeric), ~ . != 0))) %>%
  ungroup()
colSums(INTER == 0, na.rm = TRUE) # counts the number of zeros in each of the columns

# calculate the BACI at the quadrat scale
Quadrat_T2_INTER <- 
  INTER %>% 
  group_by(Taxon, site_name, block) %>% 
  mutate(t1_inter = INTER_1 - Cont_1) %>%
  mutate(t0_inter = INTER_0 - Cont_0) %>%
  mutate(BACI_inter = t1_inter - t0_inter) %>%
  select(Taxon, site_name, block, BACI_inter) %>%
  mutate(group_var = "quadrat") %>%
  mutate(time = "T0-T2") %>%
  mutate(trt = "INTER") %>%
  rename(BACI = BACI_inter) 

n_distinct(Quadrat_T2_INTER$Taxon) # 622 species

T2_result <- rbind(Quadrat_T2_DIST, Quadrat_T2_NPK, Quadrat_T2_INTER)
write_csv(T2_result, "results/July_2025/BACI_T2_results.csv")

################################################################################

# T0-T3

# read in all DRAGNet clean data
T3 <- read_csv("results/July_2025/DRAGNet_all_clean_data_with_zeros_T0_T3.csv")

# the simple BACIs are row wise calculations so we need to convert the data
wanted <- c("Control", "Disturbance")
DIST <- T3 %>% 
  mutate(site_taxon = paste0(site_name, sep = "_", New_taxon, sep = "_", year_trt)) %>%
  filter(trt %in% wanted) %>%
  group_by(site_name, New_taxon, year_trt, trt, block) %>% 
  summarise(mean_cover = mean(new_max_cover)) %>%
  mutate(trt_time = paste0(trt, sep = "_", year_trt)) %>%
  ungroup() %>%
  dplyr::select(site_name, New_taxon, block, trt_time, mean_cover) %>%
  group_by(site_name, New_taxon, block) %>%
  pivot_wider(values_from = mean_cover, names_from = trt_time)

# rename the columns to easier names
names(DIST) <- c("site_name", "Taxon", "block", "Cont_0", "Dist_0",  "Cont_1", "Dist_1")
# some rows now have NA values in them as we have done pivot wider
# replace with 0
colSums(is.na(DIST)) # there are a lot of NAs 
DIST <- DIST %>%
  mutate(across(everything(), ~ replace_na(., 0)))

# also get rid of rows where all the inputs are zeros them as they add nothing to the analysis
DIST <- DIST %>%
  rowwise() %>%
  filter(any(across(where(is.numeric), ~ . != 0))) %>%
  ungroup()
colSums(DIST == 0, na.rm = TRUE) # counts the number of zeros in each of the columns

# calculate the BACI at the quadrat scale
Quadrat_T3_DIST <- 
  DIST %>% 
  group_by(Taxon, site_name, block) %>% 
  mutate(t1_dist = Dist_1 - Cont_1) %>%
  mutate(t0_dist = Dist_0 - Cont_0) %>%
  mutate(BACI_dist = t1_dist-t0_dist) %>%
  select(Taxon, site_name, block, BACI_dist) %>%
  mutate(group_var = "quadrat") %>%
  mutate(time = "T0-T3") %>%
  mutate(trt = "DIST") %>%
  rename(BACI = BACI_dist)

n_distinct(Quadrat_T3_DIST$Taxon) # 627

# the simple BACIs are row wise calculations so we need to convert the data
wanted <- c("Control", "NPK")
NPK <- T3 %>% 
  mutate(site_taxon = paste0(site_name, sep = "_", New_taxon, sep = "_", year_trt)) %>%
  filter(trt %in% wanted) %>%
  group_by(site_name, New_taxon, year_trt, trt, block) %>% 
  summarise(mean_cover = mean(new_max_cover)) %>%
  mutate(trt_time = paste0(trt, sep = "_", year_trt)) %>%
  ungroup() %>%
  dplyr::select(site_name, New_taxon, block, trt_time, mean_cover) %>%
  group_by(site_name, New_taxon, block) %>%
  pivot_wider(values_from = mean_cover, names_from = trt_time)

# rename the columns to easier names
names(NPK) <- c("site_name", "Taxon", "block", "Cont_0", "NPK_0",  "Cont_1", "NPK_1")
# some rows now have NA values in them as we have done pivot wider
# replace with 0
colSums(is.na(NPK)) # there are a lot of NAs 
NPK <- NPK %>%
  mutate(across(everything(), ~ replace_na(., 0)))

# also get rid of rows where all the inputs are zeros them as they add nothing to the analysis
NPK <- NPK %>%
  rowwise() %>%
  filter(any(across(where(is.numeric), ~ . != 0))) %>%
  ungroup()
colSums(NPK == 0, na.rm = TRUE) # counts the number of zeros in each of the columns

# calculate the BACI at the quadrat scale
Quadrat_T3_NPK <- 
  NPK %>% 
  group_by(Taxon, site_name, block) %>% 
  mutate(t1_npk = NPK_1 - Cont_1) %>%
  mutate(t0_npk = NPK_0 - Cont_0) %>%
  mutate(BACI_npk = t1_npk - t0_npk) %>%
  select(Taxon, site_name, block, BACI_npk) %>%
  mutate(group_var = "quadrat") %>%
  mutate(time = "T0-T3") %>%
  mutate(trt = "NPK") %>%
  rename(BACI = BACI_npk) 

n_distinct(Quadrat_T3_NPK$Taxon) # 402

# the simple BACIs are row wise calculations so we need to convert the data
wanted <- c("Control", "NPK+Disturbance")
INTER <- T3 %>% 
  mutate(site_taxon = paste0(site_name, sep = "_", New_taxon, sep = "_", year_trt)) %>%
  filter(trt %in% wanted) %>%
  group_by(site_name, New_taxon, year_trt, trt, block) %>% 
  summarise(mean_cover = mean(new_max_cover)) %>%
  mutate(trt_time = paste0(trt, sep = "_", year_trt)) %>%
  ungroup() %>%
  dplyr::select(site_name, New_taxon, block, trt_time, mean_cover) %>%
  group_by(site_name, New_taxon, block) %>%
  pivot_wider(values_from = mean_cover, names_from = trt_time)

# rename the columns to easier names
names(INTER) <- c("site_name", "Taxon", "block", "Cont_0", "INTER_0",  "Cont_1", "INTER_1")
# some rows now have NA values in them as we have done pivot wider
# replace with 0
colSums(is.na(INTER)) # there are a lot of NAs 
INTER <- INTER %>%
  mutate(across(everything(), ~ replace_na(., 0)))

# also get rid of rows where all the inputs are zeros them as they add nothing to the analysis
INTER <- INTER %>%
  rowwise() %>%
  filter(any(across(where(is.numeric), ~ . != 0))) %>%
  ungroup()
colSums(INTER == 0, na.rm = TRUE) # counts the number of zeros in each of the columns

# calculate the BACI at the quadrat scale
Quadrat_T3_INTER <- 
  INTER %>% 
  group_by(Taxon, site_name, block) %>% 
  mutate(t1_inter = INTER_1 - Cont_1) %>%
  mutate(t0_inter = INTER_0 - Cont_0) %>%
  mutate(BACI_inter = t1_inter - t0_inter) %>%
  select(Taxon, site_name, block, BACI_inter) %>%
  mutate(group_var = "quadrat") %>%
  mutate(time = "T0-T3") %>%
  mutate(trt = "INTER") %>%
  rename(BACI = BACI_inter) 

n_distinct(Quadrat_T3_INTER$Taxon) # 402 species

T3_result <- rbind(Quadrat_T3_DIST, Quadrat_T3_NPK, Quadrat_T3_INTER)
write_csv(T3_result, "results/July_2025/BACI_T3_results.csv")

################################################################################


both <- rbind(T1_result, T2_result, T3_result)
head(both)
write_csv(both, "results/July_2025/BACI_responses_simple_quadrat_all_DRAGNet.csv")
 