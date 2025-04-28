library(tidyverse)
library(arrow)
library(glue)

source('scripts/helpers.R')

### Directory where EHR data is stored
ehr_dir <- '/n/haneuse_ehr_l3/V1.0'
data_dir <- '/n/haneuse_ehr_l3/V1.0/clean_datasets'

### Load in Pre-Saved Data
trial_range <- 1:84
df_trials <- 
  read_parquet(glue('{data_dir}/bariatric_tte/all_trials_combined.parquet')) %>% 
  mutate('elig_a1c_missing' = !is.na(baseline_a1c)) %>% 
  filter(trial_id %in% trial_range) 



### Stage 0: 
df0 <- 
  df_trials %>% 
  group_by(surgery) %>% 
  summarise('n_subjects' = n_distinct(subject_id),
            'n_trials' = n())

### Elig 1
df1 <- 
  df_trials %>% 
  filter(elig_age, elig_a1c_missing, elig_diabetes, elig_surgery, elig_bmi_missing, elig_bmi_range, elig_cvd, elig_cancer, elig_pregnancy) %>% 
  group_by(surgery) %>% 
  summarise('n_subjects' = n_distinct(subject_id),
            'n_trials' = n())

### Elig 2
df2 <- 
  df_trials %>% 
  filter(elig_age, elig_a1c_missing, elig_diabetes, elig_surgery, elig_bmi_missing, elig_bmi_range, elig_cvd, elig_cancer, elig_pregnancy) %>% 
  filter(elig_smoking, elig_bmi_change) %>% 
  group_by(surgery) %>% 
  summarise('n_subjects' = n_distinct(subject_id),
            'n_trials' = n())

### Elig 3
df3 <- 
  df_trials %>% 
  filter(elig_age, elig_a1c_missing, elig_diabetes, elig_surgery, elig_bmi_missing, elig_bmi_range, elig_cvd, elig_cancer, elig_pregnancy) %>% 
  filter(elig_smoking, elig_bmi_change) %>% 
  filter(elig_exclusions) %>% 
  group_by(surgery) %>% 
  summarise('n_subjects' = n_distinct(subject_id),
            'n_trials' = n())
