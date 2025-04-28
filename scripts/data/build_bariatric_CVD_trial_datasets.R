library(tidyverse)
library(arrow)
library(glue)
library(haven)
library(lubridate)

### Directory where EHR data is stored
ehr_dir <- '/n/haneuse_ehr_l3/V1.0'
data_dir <- '/n/haneuse_ehr_l3/V1.0/clean_datasets'

### Eligibility Criteria for Trial 1 (Surgery Onwards)
### (1) Enrollment during years 2007-2019
### (2) Age 18 to 65
### (3) TD2M 
### (4) No Prior Surgery 
### (5) No Current Smoking Status [Not in Trial 2]
### (6) No recent weight gain (increase in BMI < 2 kg/m^2 compared to 1 year pre-baseline) [Not in Trial 2]
### (7) No history of primary outcome 
### (8) Lab Values w/in certain ranges (including A1c change)
### (9) Vital signs w/in certain ranges (such as BMI < 60, and BP range)
### (10) BMI w/ in certain categories 
### (11) One of the exlcusion criteria

### Load in all data
df_subjects <- read_parquet(glue('{data_dir}/microvascular_tte/subjects.parquet'))
df_enrollment <- read_parquet(glue('{data_dir}/microvascular_tte/enrollment.parquet'))
diabetes_rx <-  
  bind_rows(
    read_parquet(glue('{ehr_dir}/parquet_files/controls/rx_diabetes_controls.parquet')) %>% 
      rename('subject_id' = control_studyid),
    read_parquet(glue('{ehr_dir}/parquet_files/cases/all_diabetes_rx.parquet')) %>% 
      rename('subject_id' = durable_studyid)
  )
diabetes_dx <-  read_parquet(glue('{data_dir}/microvascular_tte/diabetes_dx.parquet'))
diabetes_labs <- read_parquet(glue('{data_dir}/microvascular_tte/diabetes_labs.parquet'))
smoking <- read_parquet(glue('{data_dir}/microvascular_tte/smoking.parquet'))
weights <- read_parquet(glue('{data_dir}/all_weights_further_cleaned.parquet')) ### Cleaned Weights w/ outliers removed
cvd_outcomes <- read_parquet(glue('{data_dir}/bariatric_tte/cvd_outcomes.parquet'))
obesity_comorbidities <- read_parquet(glue('{data_dir}/bariatric_tte/obesity_comorbidities.parquet'))
exclusions <- read_parquet(glue('{data_dir}/bariatric_tte/exclusion_criteria.parquet'))
exclusions_cancer <- read_parquet(glue('{data_dir}/bariatric_tte/exclusion_cancer.parquet'))
pregnancy <- read_parquet(glue('{data_dir}/microvascular_tte/pregnancy.parquet'))
death <- read_parquet(glue('{data_dir}/microvascular_tte/death.parquet'))
df_surgery <- read_parquet(glue('{data_dir}/bs_types_reviewed.parquet')) ### bs types

### Pre-Process Outcomes
df_cvd <- 
  cvd_outcomes %>% 
  group_by(subject_id) %>% 
  summarise('cvd_date' = min(adate))

build_cvd_trial_1 <- function(trial_id, study_start, study_end) {
  options(dplyr.summarise.inform = F)
  
  ### Trial "Enrollment" Period
  trial_start <- study_start %m+% months(trial_id - 1)
  trial_end <- study_start %m+% months(trial_id) - 1
  
  ### BMI/A1c Look Backs of 1 Year
  a1c_lookback <- years(1)
  a1c_start <- trial_start %m-% a1c_lookback
  
  bmi_lookback <- years(1)
  bmi_start <- trial_start %m-% bmi_lookback
  
  ### Obesity Comordbit Lookback
  obesity_lookback <- years(1)
  obs_start <- trial_start %m-% obesity_lookback
  
  ### Pre-Processing of Diabetes Info 
  ### Most Recent Lab within 1 year and most stale value w/in a year
  df_a1c <-
    diabetes_labs %>% 
    filter(lab_date <= trial_end, lab_date >= a1c_start) %>% 
    # Remove extreme A1c Values
    filter((test_type == 'HGBA1C' & result <= 11) | (test_type == 'GLU_F' & (result + 46.7)/28.7 <= 11)) %>% 
    arrange(lab_date) %>% 
    group_by(subject_id, test_type) %>% 
    summarise('baseline_value' = last(result),
              '1yr_value' = first(result)) %>% 
    ungroup() %>% 
    pivot_wider(names_from = test_type,
                values_from = c('baseline_value', '1yr_value')) %>% 
    rename('baseline_gluF' = baseline_value_GLU_F,
           'baseline_a1c' = baseline_value_HGBA1C,
           'lag_1yr_gluF' = `1yr_value_GLU_F`,
           'lag_1yr_a1c' = `1yr_value_HGBA1C`) %>% 
    mutate('baseline_a1c' = ifelse(is.na(baseline_a1c), (baseline_gluF + 46.7)/28.7, baseline_a1c),
           'baseline_gluF' = ifelse(is.na(baseline_gluF), baseline_a1c * 28.7 - 46.7, baseline_gluF),
           'lag_1yr_a1c' = ifelse(is.na(lag_1yr_a1c), (lag_1yr_gluF + 46.7)/28.7, lag_1yr_a1c),
           'lag_1yr_gluF' = ifelse(is.na(lag_1yr_gluF), lag_1yr_a1c * 28.7 - 46.7, lag_1yr_gluF)) %>% 
    mutate('a1c_change_1yr' = baseline_a1c - lag_1yr_a1c)
  
  df_rx <- 
    diabetes_rx %>% 
    filter(rxdate <= trial_end, rxdate >= a1c_start) %>% 
    mutate('rx_end_date' = rxdate + rxsup) %>% 
    group_by(subject_id) %>% 
    summarise('baseline_insulin' = any(insulin_flg == 1 & rx_end_date >= trial_start, na.rm = T),
              'baseline_diabetes_rx' = any(rx_end_date >= trial_start, na.rm = T),
              'past_1yr_insulin' = as.numeric(any(insulin_flg == 1, na.rm = T)),
              'past_1yr_diabetes_rx' = 1) %>% 
    mutate('diabetes_rx' = as.numeric(baseline_diabetes_rx)) %>% 
    ungroup()
  
  df_dx <-
    diabetes_dx %>% 
    filter(adate <= trial_end, adate >= a1c_start) %>% 
    filter(grepl('250', dx)) %>% 
    group_by(subject_id) %>% 
    summarise('diabetes_dx' = 1) %>% 
    ungroup()
  
  ### Pre-Processing Smoking 
  df_smoking <- 
    smoking %>% 
    filter(contact_date <= trial_end) %>% 
    arrange(desc(contact_date)) %>% 
    group_by(subject_id)  %>% 
    slice(1) %>% 
    ungroup()
  
  ### Pre-Process Weights
  df_bmi <- 
    weights %>% 
    filter(bmi >= 15, bmi <= 70) %>% ### Remove Non-Sensical BMI
    filter(measure_date <= trial_end, measure_date >= bmi_start) %>% 
    filter(measure_date <= index_date | is.na(index_date)) %>% 
    arrange(measure_date) %>% 
    group_by(subject_id) %>% 
    summarise('baseline_bmi' = last(bmi),
              'lag_1yr_bmi' = first(bmi),
              'max_bmi' = max(bmi),
              'min_bmi' = min(bmi))%>% 
    mutate('bmi_change_1yr' = baseline_bmi - lag_1yr_bmi)
  
  ### Pre-Process Obesity Comorbidities
  df_obesity <- 
    obesity_comorbidities %>% 
    filter(adate <= trial_end, adate >= obs_start) %>% 
    group_by(subject_id) %>% 
    summarise('obstructive_sleep_apnea' = max(obstructive_sleep_apnea),
              'hypertension' = max(hypertension),
              'dyslipidemia' = max(dyslipidemia),
              'osteoarthritis' = max(osteoarthritis))
  
  ### Pre-Process Exclusions
  df_exclude <- 
    exclusions %>% 
    filter(adate < trial_start) %>% 
    rename('exclude_date' = adate)
  
  df_exclude_cancer <- 
    exclusions_cancer %>% 
    filter(adate < trial_start) %>% 
    rename('cancer_date' = adate)
  
  ### Pre-Process Pregnancy
  df_pregnancy <- 
    pregnancy %>% 
    filter(adate < trial_start, adate >= trial_start %m-% years(1)) %>% 
    group_by(subject_id) %>% 
    summarise('pregnancy_date' = min(adate))
  
  ### Build Trial Dataset
  df_trial <- 
    df_subjects %>% 
    mutate('trial_id' = trial_id) %>% 
    mutate('race_black' = as.numeric(race == 'BA')) %>% 
    
    ### Obesity Comorbidities
    left_join(df_obesity, by = 'subject_id') %>% 
    mutate('obstructive_sleep_apnea' = ifelse(is.na(obstructive_sleep_apnea), 0, obstructive_sleep_apnea),
           'hypertension' = ifelse(is.na(hypertension), 0, hypertension),
           'dyslipidemia' = ifelse(is.na(dyslipidemia), 0, dyslipidemia),
           'osteoarthritis' = ifelse(is.na(osteoarthritis), 0, osteoarthritis)) %>% 
    
    ### (1) Ensure at least 1 year of continuous enrollment before trial start 
    inner_join(df_enrollment, by = 'subject_id') %>% 
    filter(trial_start >= enr_1yr, trial_start <= enr_end) %>% 
    
    ### (2) Age 18 to 65
    mutate('baseline_age' = as.numeric(trial_start - birth_date)/365.25) %>% 
    mutate('elig_age' = baseline_age >= 18 & baseline_age < 65) %>% 
    
    ### (3) Diabetes Status
    left_join(df_dx, by = 'subject_id') %>% 
    left_join(df_rx, by = 'subject_id') %>% 
    left_join(df_a1c, by = 'subject_id') %>% 
    mutate('diabetes' = case_when(!is.na(diabetes_dx) & diabetes_dx == 1 ~ 1,
                                  !is.na(diabetes_rx) & diabetes_rx == 1 ~ 1,
                                  !is.na(baseline_a1c) & baseline_a1c >= 6.5 ~ 1,
                                  !is.na(baseline_gluF) & baseline_gluF >= 126 ~ 1,
                                  T ~ 0)) %>% 
    mutate('elig_diabetes' = diabetes == 1) %>% 
    mutate('past_1yr_insulin' = replace(past_1yr_insulin, is.na(past_1yr_insulin), 0),
           'past_1yr_diabetes_rx' = replace(past_1yr_diabetes_rx, is.na(past_1yr_diabetes_rx), 0)) %>% 
    
    ### (4) No Prior Surgery (and code treatment) 
    left_join(df_surgery, by = c('subject_id', 'index_date')) %>% 
    mutate('surgery' = as.numeric(!is.na(index_date) & index_date >= trial_start & index_date <= trial_end)) %>% 
    mutate('elig_surgery' = is.na(index_date) | (index_date >= trial_start & bs_type != 'EXCLUDE')) %>% 
    
    ### (5) Smoking Status. Trial 1 Includes Never/Former, Trial 2 has no restrictions
    left_join(df_smoking, by = 'subject_id') %>% 
    mutate('smoking_status' = ifelse(is.na(smoking_status), 'no_self_report', smoking_status)) %>% 
    mutate('elig_smoking' = smoking_status != 'current') %>% 
    
    ### (6)/(10) BMI
    left_join(df_bmi, by = 'subject_id') %>% 
    mutate('elig_bmi_missing' = !is.na(baseline_bmi)) %>% 
    mutate('elig_bmi_range' = max_bmi >= 35 & baseline_bmi <= 60 & baseline_bmi >= 30) %>%  ### They have diabetes as a comorbidity so that's already covered
    mutate('elig_bmi_change' = bmi_change_1yr < 2) %>%
    
    ### (7) No History of Outcome/Code Time to Outcome
    left_join(df_cvd, by = 'subject_id') %>% 
    mutate('elig_cvd' = cvd_date >= trial_end | is.na(cvd_date)) %>% 
    
    ### (11) Exclusions 
    left_join(df_exclude, by = 'subject_id') %>% 
    left_join(df_exclude_cancer, by = 'subject_id') %>% 
    mutate('elig_exclusions' = is.na(exclude_date)) %>% 
    mutate('elig_cancer' = is.na(cancer_date)) %>% 
    
    ### (12) Pregnant in Last 12 Months
    left_join(df_pregnancy, by = 'subject_id')  %>% 
    mutate('elig_pregnancy' = is.na(pregnancy_date)) %>% 
    
    ### Censoring due to Death
    left_join(death, by = 'subject_id') %>% 
    filter(is.na(censor_death) | censor_death >= trial_start) %>%
    
    ### Compute Event Time
    mutate('study_end' = study_end,
           'trial_start' = study_start %m+% months(trial_id - 1)) %>% 
    mutate('censor_date' = 
             case_when(surgery == 0 ~ pmin(index_date, censor_death, enr_end, study_end, na.rm = T),
                       surgery == 1 ~ pmin(enr_end, censor_death, study_end, na.rm = T)),
           'cvd_event' = case_when(is.na(cvd_date) ~ 0,
                                   cvd_date > censor_date ~ 0,
                                   cvd_date <= censor_date ~ 1),
           'cvd_time' = case_when(cvd_event == 0 ~ as.numeric(censor_date - trial_start)/30.25,
                                  cvd_event == 1 ~ as.numeric(cvd_date - trial_start)/30.25)) %>% 
    
    
    ### Select Covariates 
    select(subject_id, trial_id, surgery, trial_start, index_date, contains('elig'),
           baseline_age, gender, race_black, smoking_status,
           baseline_bmi, bmi_change_1yr, baseline_a1c, a1c_change_1yr,
           past_1yr_insulin, past_1yr_diabetes_rx, osteoarthritis, 
           hypertension, dyslipidemia, obstructive_sleep_apnea,
           enr_start, enr_end, study_end, censor_death, cvd_date, censor_date, cvd_time, cvd_event)
  
  return(df_trial)
}

### Build Trial Dataset
args <- commandArgs(trailingOnly = T)
trial_id <- as.numeric(args[1])

df_trial <- 
  build_cvd_trial_1(trial_id, 
                    study_start = as.Date('2005-01-01'), 
                    study_end = as.Date('2015-09-30'))

write_parquet(df_trial, glue('{data_dir}/bariatric_tte/trial_{trial_id}_tte1.parquet'))
