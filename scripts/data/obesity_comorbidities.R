library(tidyverse)
library(arrow)
library(glue)
library(haven)
library(lubridate)

### Directory where EHR data is stored
ehr_dir <- '/n/haneuse_ehr_l3/V1.0'
data_dir <- '/n/haneuse_ehr_l3/V1.0/clean_datasets'


### Obesity Comborbidities: Dsylipedmia, hypertension, osteoarthritis, obstructive sleep apnea
### Obstructive Sleep Apnea
cardiac_dx_cases <- 
  read_sas(glue('{ehr_dir}/_data_for_eric_cases/raw_cardiac_dx.sas7bdat')) %>% 
  rename('subject_id' = 'durable_studyid')
cardiac_dx_controls <- 
  read_sas(glue('{ehr_dir}/_data_for_eric_controls_001/raw_cardiac_dx_controls.sas7bdat')) %>% 
  rename('subject_id' = 'control_studyid')


df_osa <- 
  bind_rows(cardiac_dx_cases, cardiac_dx_controls) %>% 
  mutate('obstructive_sleep_apnea' = as.numeric(dx == '327.23')) %>% 
  filter(obstructive_sleep_apnea == 1) %>% 
  select(subject_id, adate, obstructive_sleep_apnea)

### Dyslipidemia + Hypertension
dht_cases <- 
  read_sas(glue('{ehr_dir}/_data_for_eric_cases/hypertension_lipid_dx_case.sas7bdat')) %>% 
  rename('subject_id' = 'durable_studyid')

dht_controls <- 
  read_sas(glue('{ehr_dir}/_data_for_eric_controls_001/hypertension_lipid_dx_controls.sas7bdat')) %>% 
  rename('subject_id' = 'control_studyid')

df_dht <- 
  dht_cases %>% 
  rename('dx' = 'DX', 'adate' = 'ADATE') %>%
  bind_rows(dht_controls) %>% 
  mutate('dx_numeric' = as.numeric(dx)) %>% 
  mutate('hypertension' = as.numeric(dx_numeric >= 401 & dx_numeric < 406)) %>% 
  mutate('dyslipidemia' = as.numeric((dx_numeric >= 272 & dx_numeric <= 272.2) | (dx_numeric >= 272.4 & dx_numeric <= 272.5))) %>%
  filter(hypertension == 1 | dyslipidemia == 1) %>% 
  select(subject_id, 'adate' = adate, dyslipidemia, hypertension)

### Osteoarthritis (using hospital_dx for now)
hospital_cases <- 
  read_parquet(glue('{ehr_dir}/parquet_files/cases/hospital_dx.parquet')) %>% 
  rename('subject_id' = 'durable_studyid')

hospital_controls <- 
  read_parquet(glue('{ehr_dir}/parquet_files/controls/hospital_dx_controls.parquet')) %>% 
  rename('subject_id' = 'control_studyid')

df_arthritis <- 
  bind_rows(hospital_controls, hospital_cases) %>% 
  mutate('osteoarthritis' = as.numeric(grepl('^715', dx))) %>% 
  filter(osteoarthritis == 1) %>% 
  select(subject_id, adate, osteoarthritis)

### Final Dataset
df_obesity <- 
  df_osa %>% 
  bind_rows(df_dht) %>% 
  bind_rows(df_arthritis) %>% 
  group_by(subject_id, adate) %>% 
  summarise('obstructive_sleep_apnea' = max(0, obstructive_sleep_apnea, na.rm = T),
            'hypertension' = max(0, hypertension, na.rm = T),
            'dyslipidemia' = max(0, dyslipidemia, na.rm = T),
            'osteoarthritis' = max(0, osteoarthritis, na.rm = T)) %>% 
  ungroup()

if(!dir.exists(glue('{data_dir}/bariatric_tte'))) {
  dir.create(glue('{data_dir}/bariatric_tte'))
}

write_parquet(df_obesity, glue('{data_dir}/bariatric_tte/obesity_comorbidities.parquet')) 
