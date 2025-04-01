library(tidyverse)
library(arrow)
library(glue)
library(haven)
library(lubridate)

### Directory where EHR data is stored
ehr_dir <- '/n/haneuse_ehr_l3/V1.0'
data_dir <- '/n/haneuse_ehr_l3/V1.0/clean_datasets'


### Outcome: Coronary Artery Disease, Stroke, Coronary Artery Stenting, Coronary Artery Bypass Grafting
cardiac_dx_cases <- 
  read_sas(glue('{ehr_dir}/_data_for_eric_cases/raw_cardiac_dx.sas7bdat')) %>% 
  rename('subject_id' = 'durable_studyid')
cardiac_px_cases <- read_sas(glue('{ehr_dir}/_data_for_eric_cases/raw_cardiac_px_cases.sas7bdat')) %>% 
  rename('subject_id' = 'durable_studyid')

cardiac_dx_controls <- 
  read_sas(glue('{ehr_dir}/_data_for_eric_controls_001/raw_cardiac_dx_controls.sas7bdat')) %>% 
  rename('subject_id' = 'control_studyid')
cardiac_px_controls <- 
  read_sas(glue('{ehr_dir}/_data_for_eric_controls_001/raw_cardiac_px_controls.sas7bdat')) %>% 
  rename('subject_id' = 'control_studyid')

df_cvd <- 
  bind_rows(
    bind_rows(cardiac_dx_cases, cardiac_dx_controls) %>% 
      mutate('congestive_heart_failure' = as.numeric(grepl('^428', dx)),
             'coronary_artery_disease' = as.numeric(as.numeric(dx) >= 410 & as.numeric(dx) < 415),
             'cerebrovascular_accident' = as.numeric(grepl('^431', dx) | grepl('^433', dx) | grepl('^434', dx) | grepl('^435', dx))) %>% 
      filter(cerebrovascular_accident == 1 | coronary_artery_disease == 1 | congestive_heart_failure == 1) %>% 
      select(subject_id, adate) %>% 
      mutate('cvd' = 1),
    
    bind_rows(cardiac_px_cases, cardiac_px_controls) %>% 
      mutate('coronary_artery_bypass_graft' = as.numeric(px %in% c('36.03', '36.04') | grepl('^36.1', px) | grepl('^36.2', px) | px %in% as.character(c(33510:33514, 33516:33523, 33533:33536, 33572))),
             'cornary_stent_placement' = as.numeric(px %in% c('36.06', '36.07', '36.09', '00.66') | px %in% as.character(c(92980:92982, 92984:92996, 92929, 92933, 92934, 92941, 92943, 92944, 92973:92975, 92977)))) %>% 
      filter(coronary_artery_bypass_graft == 1 | cornary_stent_placement == 1) %>% 
      select(subject_id, adate) %>% 
      mutate('cvd' = 1)
  ) %>% 
  distinct()


if(!dir.exists(glue('{data_dir}/bariatric_tte'))) {
  dir.create(glue('{data_dir}/bariatric_tte'))
}

write_parquet(df_cvd, glue('{data_dir}/bariatric_tte/cvd_outcomes.parquet')) 
              