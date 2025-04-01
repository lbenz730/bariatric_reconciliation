library(tidyverse)
library(arrow)
library(glue)
library(haven)
library(lubridate)

### Directory where EHR data is stored
ehr_dir <- '/n/haneuse_ehr_l3/V1.0'
data_dir <- '/n/haneuse_ehr_l3/V1.0/clean_datasets'

### Read in Hopsital DX File for Several Exclusions
hospital_cases <- 
  read_parquet(glue('{ehr_dir}/parquet_files/cases/hospital_dx.parquet')) %>% 
  rename('subject_id' = 'durable_studyid')

hospital_controls <- 
  read_parquet(glue('{ehr_dir}/parquet_files/controls/hospital_dx_controls.parquet')) %>% 
  rename('subject_id' = 'control_studyid')


exclude_hospital <- 
  hospital_controls %>% 
  bind_rows(hospital_cases) %>% 
  mutate('dx_numeric' = as.numeric(dx)) %>% 
  filter( (dx_numeric >= 490 & dx_numeric < 493) | ### COPD
            grepl('^533', dx) | ### Peptic Ulcer
            (dx_numeric >= 555 & dx_numeric < 557) | ### Inflammatory Bowel Disease
            grepl('^570', dx) | ### Liver Failure
            (dx %in% c('070.2', '070.3')) | ###  Hepatitis B
            (grepl('^806', dx) | grepl('^952', dx))) %>%  ### Spinal Cord Injury 
  mutate('exclusion' = case_when(
    (dx_numeric >= 490 & dx_numeric < 493) ~ 'Chronic Obstructive Pulmonary Disease',
    grepl('^533', dx) ~ 'Peptic Ulcer',
    (dx_numeric >= 555 & dx_numeric < 557) ~ 'Inflammatory Bowel Disease',
    grepl('^570', dx) ~ 'Liver Failure',
    dx %in% c('070.2', '070.3') ~ 'Hepatitis B',
    grepl('^806', dx) | grepl('^952', dx) ~ 'Spinal Cord Injury')) %>% 
  select(subject_id, adate, dx, exclusion)

### MH Exlcusions
mh_controls <- 
  read_parquet(glue('{ehr_dir}/parquet_files/controls/mental_health_dx_controls.parquet')) %>% 
  rename('subject_id' = 'control_studyid')
mh_cases <- 
  read_parquet(glue('{ehr_dir}/parquet_files/cases/mental_health_dx_cases.parquet')) %>% 
  rename('subject_id' = 'durable_studyid')

exclude_mh <- 
  mh_controls %>% 
  bind_rows(mh_cases) %>% 
  filter((grepl('^296', dx) | grepl('^311', dx)) | ### Depression/BiPolor
           (grepl('^303', dx) | grepl('^304', dx)) | ### Substance Use Disorder
           (grepl('^295', dx) | grepl('^297', dx) | grepl('^298', dx) | grepl('^299', dx))) %>%  ### Psycosis
  mutate('exclusion' = case_when(
    grepl('^296', dx) | grepl('^311', dx) ~ 'Depression/Bipolar',
    grepl('^303', dx) | grepl('^304', dx) ~ 'Substance Use Disorder',
    grepl('^295', dx) | grepl('^297', dx) | grepl('^298', dx) | grepl('^299', dx) ~ 'Psycosis')) %>% 
  select(subject_id, adate, dx, exclusion)

### Cancer Exclusions
case_cancer <- 
  read_sas(glue('{ehr_dir}/_data_for_eric_cases/raw_cancer_dx_cases.sas7bdat')) %>% 
  rename('subject_id' = 'durable_studyid')
control_cancer <- 
  read_sas(glue('{ehr_dir}/_data_for_eric_controls_001/raw_cancer_dx_controls.sas7bdat')) %>% 
  rename('subject_id' = 'control_studyid')

exclude_cancer <- 
  control_cancer %>% 
  bind_rows(case_cancer) %>% 
  filter(as.numeric(dx) >= 140 & as.numeric(dx) < 210) %>% 
  mutate('exclusion' = case_when(
    as.numeric(dx) >= 173 & as.numeric(dx) < 174 ~ 'Non-Melanoma Skin Cancer',
    T ~ 'Cancer')) %>% 
  select(subject_id, adate, dx, exclusion)

df_exclude <- 
  bind_rows(exclude_cancer, exclude_hospital, exclude_mh) %>% 
  arrange(adate) %>% 
  group_by(subject_id) %>% 
  slice(1) %>% 
  ungroup()

write_parquet(df_exclude, glue('{data_dir}/bariatric_tte/exclusion_criteria.parquet'))


### Raw Exclusions
raw_exclusions <- bind_rows(exclude_cancer, exclude_hospital, exclude_mh)
write_parquet(raw_exclusions, glue('{data_dir}/bariatric_tte/raw_exclusion_criteria.parquet'))
