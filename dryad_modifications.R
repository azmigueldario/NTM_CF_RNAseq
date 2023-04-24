####################################### Dryad submission clean-up ###################################### 

# Load the necessary libraries, mostly tidyverse
library(tidyverse)

 
# We want to remove any identifiers of participants for submission to the Dryad repository. 
# In this script, I clean-up all primary data and let only the necessary variables to produce all the plots and 
# reproduce all analyses done in the manuscript. 
# 
# 
# We will use the backup dataset backup_ntm_data.csv as the building block

ntm_data <- 
  read_csv(file = "./backup/backup_ntm_data.csv") %>% 
  mutate(across(where(is.character), 
                as_factor))

############################################# Culture Data ############################################# 

# ---------------------  Data wrangling primary culture data        

temp_data <-  list()
temp_data$wide_culture<- read_csv("./backup/NTM_cultures.csv", 
                                  col_types = cols(CFB_study_id = col_factor(),
                                                   consent_signed = col_skip(),
                                                   Other_ID = col_skip(),
                                                   ntm_growth = col_skip(),
                                                   NTM_ID = col_skip()
                                                   )
                                  ) 

# remove those without RNA seq data
temp_data$wide_culture <- 
  temp_data$wide_culture %>% 
  filter(!is.na(rna_sample_date))

# pivot to long 
temp_data$culture_long_data <-
  temp_data$wide_culture %>% 
  pivot_longer(cols= -c(CFB_study_id, ntm_disease),
               names_to=c(".value", "set"),
               # matching everything before and numbers after it
               names_pattern="(.*)_([0-9])",
               values_drop_na = T) %>% 
  select(-set) # eliminate the set variable, corresponds to set # of cultures 

# recode labels
temp_data$culture_long_data <- 
  temp_data$culture_long_data %>% 
  mutate(culture=recode_factor(culture,
                               "0"="Negative",
                               "1"="MAC",
                               "2"="MABs",
                               "3"="M_gordonae",
                               "4"="M_chelonae",
                               "5"="M_fortuitum",
                               "6"="M_cosmeticum",
                               "7"="M_peregrinum",
                               "8"="M_simiae",
                               "9"="M_aurepense",
                               "10"="M_scrofulaceum", 
                               "11"= "Unspecified"))


# ---------------------  Create culture var with all NTM species 

# create var with all ntm species
temp_data$culture_dryad <- 
  temp_data$culture_long_data %>%
  filter(culture!="Negative") %>%
  group_by(CFB_study_id) %>%
  count(culture) %>% 
  top_n(1, n) %>% 
  mutate(
    mycobacteria_species_ungrouped = culture, 
    mycobacteria_species_grouped = case_when(culture == "MABs" ~ "MABs",
                                             culture == "MAC" ~ "MAC",
                                             TRUE ~ "other"))

# merge that variable into primary dataset
temp_data$culture_dryad <- 
  temp_data$culture_dryad %>% 
  select(CFB_study_id, mycobacteria_species_ungrouped)

# base_data$mycobacteria_species
culture_dryad <- 
  left_join(ntm_data, temp_data$culture_dryad, by="CFB_study_id")

# recode original ntm_species var to remove M. gordonae

culture_dryad <- 
  culture_dryad %>%
  mutate(
    mycobacteria_species = case_when(mycobacteria_species == "M_gordonae" ~ "other",
                                     TRUE ~ mycobacteria_species))


########################################## Lung function data ########################################## 

# ------------------------  Prepare lung_function and clinical data

# lung function percent predicted requires specific age, so load dataset with identifiers 
temp_data$pft_clinical <- 
  read_csv(file = "./backup/NTM_clinical.csv") %>%
  select(CFB_study_id, mycobacteria_species, first_course_date,
         CF_related_disorder, ethnicity, age_years,
         transient_inf, ntm_disease, ntm_disease_date, sample_date) %>%
  rename(rna_sample_date="sample_date")

# load dataset with lung function data
temp_data$lung_function_dat <- 
  read_csv(file="./backup/NTM_lung_function.csv") %>%
  mutate(across(where(is.character), as_factor)) %>% 
  mutate(female=as_factor(female))

# match everything in lung_function_dat object using left join
pft_merged <- 
  left_join(temp_data$lung_function_dat,
            temp_data$pft_clinical,
            by = c("CFB_study_id", "ethnicity")) %>% 
  # eliminate if participant does not have RNAseq data
  filter(!is.na(rna_sample_date)) 

# verify that there are 42 participants/levels 
nlevels(as_factor(pft_merged$CFB_study_id))

# ------------------------- Prepare data for use in the rspiro package 

# mutate ntm_groups and ethnicity as required by rspiro
pft_merged <-  
  pft_merged %>% 
  mutate(height_mts = height_cms / 100,
         ethnicity = case_when(ethnicity == "caucasian" ~ 1,
                               ethnicity == "asian" ~ 4,
                               TRUE ~ 5)
         ) %>%
  mutate(across(where(is.character), 
                as_factor))

# check where lung function data ends
colnames(pft_merged)[190:200]

# format dates
pft_merged <- 
  pft_merged %>% 
  mutate(across(contains("date"), ymd)) 

# pivot to long by CFB_ID for analysis
pft_long <- 
  pft_merged %>%
  pivot_longer(cols= starts_with(c("FEV","FVC","date")),
               # new variables to create from variable names
               names_to=c(".value","set"),
               # patterns to match to value and set
               names_pattern="(FEV_L|FVC_L|date)_(.*)",
               values_to="value_L",
               values_drop_na = T) %>%
  # for later functions we factorize CFTR_modulator
  mutate(modulator_ind = case_when(!is.na(start_modulator) ~ "yes",
                                   TRUE ~ "no")) 

# ---------------------- Apply GLI equations in rspiro package 


# install.packages("devtools")
# devtools::install_github("thlytras/rspiro")


library(rspiro)

# create variable showing interval between pft measure and birth
pft_long <- 
  mutate(pft_long,
         age_years_pft = time_length(difftime(date, year_birth), "years"),
         age_years_pft = floor(age_years_pft)) 

# calculate percent predicted values by age/ethnicity/height/sex
pft_long<-  
  pft_long %>%
  mutate(FEVpp = pctpred_GLI(age = age_years_pft,
                             ethnicity = ethnicity,
                             height = height_mts,
                             gender = female,
                             FEV1 = FEV_L),
         FVCpp = pctpred_GLI(age = age_years_pft,
                             height = height_mts,
                             gender = female,
                             FEV1 = FVC_L )) %>% 
  select(-c(height_mts, ethnicity, female, year_birth, 
            FEV_L, FVC_L, set, rna_sample_date))

# -------------------- Establish baseline lung function 

# Using measurements before the `first_pos_date`, we calculate a mean (\~2 years)
# of baseline lung function.

temp_data$baseline_pft <-
  pft_long %>%
  filter(date <= first_pos_date) %>%
  group_by(CFB_study_id) %>%
  summarise(base_ppFEV = mean(FEVpp, na.rm=T),
            base_ppFVC = mean(FVCpp, na.rm=T),
            base_ppFEV = round(base_ppFEV, 2),
            base_ppFVC = round(base_ppFVC, 2))

# merge with primary dataset from cultures
pft_dryad <- 
  left_join(culture_dryad, temp_data$baseline_pft) %>% 
  # remove height as identifier
  select(-height_cms)

# QC - no apparent outliers in ppPFT
skimr::skim(pft_dryad)

# QC - missing data at baseline for 1 patient
pft_dryad %>% filter(is.na(base_ppFEV))

# ---------------------- Additional removal of possible identifiers

pft_dryad <- 
  pft_dryad %>% 
  mutate(batch_extraction=case_when(extraction_date == "2020-10-01" ~ 1,
                                    extraction_date == "2020-10-09" ~ 2,
                                    extraction_date == "2020-10-14" ~ 3,
                                    extraction_date == "2021-01-20" ~ 4,
                                    extraction_date == "2021-04-28" ~ 5))

pft_dryad <- 
  pft_dryad %>% 
  select(-c(first_course_date, extraction_date, ntm_disease_date))


########################################### Build CBC dataset ######################################### 

cbc_temp <- 
  read_csv(file="./backup/NTM_CBC_counts.csv") %>% 
  rename("CFB_study_id"="Lab ID")

cbc_data <-  
  read_csv(file = "./backup/clinical_rnaseq.csv") %>%
  inner_join(x = cbc_temp, y = ., by = "CFB_study_id") %>%
  mutate(discordant = factor(discordant),
         batch_extraction = factor(batch_extraction)) %>%
  mutate(across(where(is.character), as_factor))

colnames(cbc_data)

cbc_data <- 
  cbc_data %>% 
  select(CFB_study_id, ntm_disease, discordant,
         matches("abs|perc")) 


####################################### Clean variables and Export #################################### 

# -------------------------- Bin BMI_baseline

pft_dryad <- 
  pft_dryad %>%
  mutate(bmi_baseline = case_when(bmi_baseline < 18.5 ~ "< 18.5",
                                  bmi_baseline < 25 ~ "18.5 to 24.9",
                                  bmi_baseline < 30 ~ "25 to 29.9",
                                  bmi_baseline >= 30 ~ "\u2265 30",
                                  TRUE ~ NA_character_ )) 

# -------------------------  Export datasets

write_csv(pft_dryad, 
          file = "ntm_analysis_dryad.csv")

write_csv(cbc_data,
          file = "ntm_cbc_dryad.csv")

