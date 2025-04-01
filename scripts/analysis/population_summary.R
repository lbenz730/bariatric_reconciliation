library(tidyverse)
library(arrow)
library(glue)
library(haven)
library(lubridate)
library(table1)
library(gt)
library(survival)
library(survminer)
library(patchwork)
library(splines)
source('scripts/helpers.R')

### Directory where EHR data is stored
ehr_dir <- '/n/haneuse_ehr_l3/V1.0'
data_dir <- '/n/haneuse_ehr_l3/V1.0/clean_datasets'

### Load in Pre-Saved Data
n_trials <- 96
# df_trials <- map_dfr(1:n_trials, ~read_parquet(glue('{data_dir}/bariatric_tte/trial_{.x}_tte1.parquet')))
# write_parquet(df_trials, glue('{data_dir}/bariatric_tte/all_trials_combined.parquet'))
df_trials <- read_parquet(glue('{data_dir}/bariatric_tte/all_trials_combined.parquet'))

trial_range <- 1:84

### Naive Eligibility Criteria
df1_naive <- 
  df_trials %>% 
  mutate('elig_a1c_missing' = !is.na(baseline_a1c)) %>% 
  filter(trial_id %in% trial_range) %>% 
  filter(elig_age, elig_a1c_missing, elig_diabetes, elig_surgery, elig_bmi_missing, elig_bmi_range, elig_cvd) %>% 
  mutate('elig' = 'No Pre-Op Restrictions')

### Pre-Operative Eligibility Criteria
df1_preop <- 
  df_trials %>% 
  mutate('elig_a1c_missing' = !is.na(baseline_a1c)) %>% 
  filter(trial_id %in% trial_range) %>% 
  filter(elig_age, elig_a1c_missing, elig_diabetes, elig_surgery, elig_bmi_missing, elig_bmi_range, elig_cvd) %>%
  filter(elig_smoking, elig_bmi_change) %>% 
  mutate('elig' = 'Pre-Op Restriction Proxy')

### Pre-Operative Eligibility Criteria + Additional Exclusion Criteria
df1_extra <- 
  df_trials %>% 
  mutate('elig_a1c_missing' = !is.na(baseline_a1c)) %>% 
  filter(trial_id %in% trial_range) %>% 
  filter(elig_age, elig_a1c_missing, elig_diabetes, elig_surgery, elig_bmi_missing, elig_bmi_range, elig_cvd) %>%
  filter(elig_smoking, elig_bmi_change) %>% 
  filter(elig_exclusions, elig_pregnancy) %>% 
  mutate('elig' = 'Pre-Op Restriction Proxy + Additional Exclusion Criteria')

df1 <-
  df1_naive %>% 
  bind_rows(df1_preop) %>% 
  bind_rows(df1_extra) %>% 
  filter(!is.na(baseline_a1c)) %>% 
  mutate('surgery' = factor(ifelse(surgery == 1, 'Bariatric Surgery', 'No Bariatric Surgery')),
         'race_black' = factor(ifelse(race_black == 1, 'African-American', 'Non African-American')),
         'past_1yr_insulin' = fct_rev(factor(ifelse(past_1yr_insulin == 1, 'Yes', 'No'))),
         'past_1yr_diabetes_rx' = fct_rev(factor(ifelse(past_1yr_diabetes_rx == 1, 'Yes', 'No'))),
         'hypertension' = fct_rev(factor(ifelse(hypertension == 1, 'Yes', 'No'))),
         'dyslipidemia' = fct_rev(factor(ifelse(dyslipidemia == 1, 'Yes', 'No'))),
         'osteoarthritis' = fct_rev(factor(ifelse(osteoarthritis == 1, 'Yes', 'No'))),
         'obstructive_sleep_apnea' = fct_rev(factor(ifelse(obstructive_sleep_apnea == 1, 'Yes', 'No'))),
         'cvd_event_7yr' = case_when(cvd_event == 1 & cvd_time <= 84 ~ 'CVD',
                                     cvd_event == 1 & cvd_time > 84 ~ 'No CVD',
                                     cvd_event == 0 & cvd_time >= 84 ~ 'No CVD',
                                     cvd_event == 0 & cvd_time <= 84 ~ 'Censored'))

label(df1$baseline_age) <- 'Baseline Age'
label(df1$gender) <- 'Sex'
label(df1$race_black) <- 'Race'
label(df1$baseline_bmi) <- 'Baseline BMI'
label(df1$baseline_a1c) <- 'Baseline HbA1c %'
label(df1$bmi_change_1yr) <- '1-Year Change in BMI'
label(df1$past_1yr_insulin) <- 'Insulin Usage w/in Past Year'
label(df1$past_1yr_diabetes_rx) <- 'Any Diabetes Medication Rx w/in Past Year'
label(df1$hypertension) <- 'Hypertension'
label(df1$obstructive_sleep_apnea) <- 'Obstructive Sleep Apnea'
label(df1$osteoarthritis) <- 'Osteoarthritis'
label(df1$dyslipidemia) <- 'Dyslipidemia'
label(df1$smoking_status) <- 'Self-Reported Smoking Status'
label(df1$cvd_event_7yr) <- 'CVD Status at 7 Years'

tbl_1 <- 
  table1(~baseline_age + gender + race_black + smoking_status + baseline_bmi + bmi_change_1yr +  baseline_a1c + a1c_change_1yr +
           past_1yr_insulin + past_1yr_diabetes_rx + hypertension + osteoarthritis + dyslipidemia + 
           obstructive_sleep_apnea + cvd_event_7yr | elig + surgery,
         data = df1, 
         render.continuous = render.continuous,
         render.strat = render.strat.subj_trials,
         render.categorical = render.categorical,
         overall = F)
tbl_1
write_lines(tbl_1, 'figures/tables/population_summary_table.html')


df_summary <- 
  df_trials %>% 
  filter(trial_id %in% trial_range) %>% 
  mutate('elig_a1c_missing' = !is.na(baseline_a1c)) %>% 
  # mutate('elig_bmi_range' = ifelse(is.na(elig_bmi_range), F, elig_bmi_range),
  # 'elig_bmi_change' = ifelse(is.na(elig_bmi_change), F, elig_bmi_change)) %>% 
  select(surgery, contains('elig')) %>% 
  pivot_longer(cols = contains('elig'),
               names_to = 'elig_criteria',
               values_to = 'elig_status') %>% 
  group_by(surgery, elig_criteria) %>% 
  summarise('n' = n(),
            'n_fail' = sum(1-elig_status, na.rm = T),
            'pct_retained' = mean(elig_status, na.rm = T)) 

df_summary <- 
  df_summary %>% 
  mutate('elig_criteria' = case_when(elig_criteria == 'elig_age' ~ '18 <= Baseline Age < 65',
                                     elig_criteria == 'elig_bmi_change' ~ '1 Year BMI Change < 2 Kg/m^2 (Among those w/ available BMI)',
                                     elig_criteria == 'elig_cvd' ~ 'No History of CVD',
                                     elig_criteria == 'elig_diabetes' ~ 'Type II Diabetes',
                                     elig_criteria == 'elig_smoking' ~ 'Not Currently Smoker',
                                     elig_criteria == 'elig_surgery' ~ 'No Previous Bariatric Surgery',
                                     elig_criteria == 'elig_bmi_missing' ~ 'Available BMI w/in 1 Year',
                                     elig_criteria == 'elig_a1c_missing' ~ 'Available A1c w/in 1 Year',
                                     elig_criteria == 'elig_pregnancy' ~ 'No Pregnancy w/in 1 Year',
                                     elig_criteria == 'elig_exclusions' ~ 'No Exclusions',
                                     elig_criteria == 'elig_bmi_range' ~ 'Basline BMI <= 60 kg/m^2, Max BMI in Last Year >= 35 kg/m^2 (Among those w/ available BMI)')) %>% 
  mutate('surgery' = paste0(ifelse(surgery == 1, 'Surgery', 'No Surgery'), ' (# of Possible Subject Trials = ', prettyNum(n, big.mark = ','), ')')) %>% 
  select(-n)

gt_elig <- 
  gt(df_summary) %>% 
  cols_align('center') %>% 
  fmt_number(n_fail, decimals = 0) %>% 
  fmt_percent(pct_retained, decimals = 1) %>% 
  cols_label('elig_criteria' = 'Eligibility Criterion',
             'pct_retained' = '% of Subject-Trials Meeting Criterion',
             'n_fail' = '# of Subject-Trials Failing to Meet Criterion') %>% 
  tab_header(title = md('**Bariatric Surgery Target Trial Eligibility Retention Table**'),
             subtitle = md('**After Restricting to Continued Enrollment for 1+ Year**')) %>% 
  tab_source_note('NOTE: Subject-Trials may fail multiple criteria at the same time') %>% 
  data_color(columns = pct_retained,
             fn = scales::col_numeric(palette = 'RdYlGn', reverse = F, domain = range(df_summary$pct_retained)),
             autocolor_text = T) %>%
  tab_options(column_labels.font.size = 20,
              heading.title.font.size = 36,
              heading.subtitle.font.size = 30,
              heading.title.font.weight = 'bold',
              heading.subtitle.font.weight = 'bold',
              column_labels.font.weight = 'bold',
              row_group.font.weight = 'bold')

cluster_gtsave(gt_elig, 'figures/tables/eligibility_summary_table.html') 


### Eligibility Breakdown by Additional Exclusions among subjects meeting pre-operative restrictions
pregnancy <- read_parquet(glue('{data_dir}/microvascular_tte/pregnancy.parquet'))
df_exclusions_raw <- 
  read_parquet(glue('{data_dir}/bariatric_tte/raw_exclusion_criteria.parquet')) %>% 
  bind_rows(
    pregnancy %>% 
      mutate('exclusion' = 'Pregnancy')
  )


df_surg <- 
  df1_preop %>% 
  filter(surgery == 1) %>% 
  select(subject_id, index_date, trial_id, trial_start)

df_exclude <- 
  df_surg %>% 
  left_join(df_exclusions_raw, by = 'subject_id', relationship = "many-to-many") %>% 
  filter((adate < trial_start & exclusion != 'Pregnancy') | 
           (adate < trial_start & adate >= trial_start %m-% years(1) & exclusion == 'Pregnancy')) %>% 
  distinct(subject_id, trial_id, exclusion)

exclude_summary <- 
  df_exclude %>% 
  group_by(exclusion) %>% 
  count() %>% 
  ungroup() %>% 
  mutate('pct' = n/nrow(df_surg)) %>% 
  arrange(-n)

gt_exclude <-
  gt(exclude_summary) %>% 
  cols_align('center') %>% 
  fmt_number('n', decimals = 0) %>% 
  fmt_percent('pct', decimals = 1) %>% 
  cols_label('exclusion' = 'Exclusion',
             'n' = '# of Patients',
             'pct' = '%') %>% 
  tab_header(title = md('**Table of Exclusions for Patients Undergoing Bariatric Surgery**'),
             subtitle = md('**Among 4,969 Patients Meeting All Other Eligibility Criteria<br>(Including Pre-Operative Restrictions)**'))
cluster_gtsave(gt_exclude, 'figures/tables/exclusion_summary_table.png', vheight = 100, vwidth = 700)


### Same thing but among all bariatric surgery patients
df_subjects <- read_parquet(glue('{data_dir}/microvascular_tte/subjects.parquet'))   
df_surg_all <- 
  df_subjects %>% 
  filter(!is.na(index_date)) 

df_exclude_all <- 
  df_surg_all %>% 
  left_join(df_exclusions_raw, by = 'subject_id') %>% 
  filter((adate < index_date & exclusion != 'Pregnancy') | 
           (adate < index_date & adate >= index_date %m-% years(1) & exclusion == 'Pregnancy')) %>% 
  distinct(subject_id, exclusion)

exclude_summary_all <- 
  df_exclude_all %>% 
  group_by(exclusion) %>% 
  count() %>% 
  ungroup() %>% 
  mutate('pct' = n/nrow(df_surg_all)) %>% 
  arrange(-n)

gt_exclude_all <-
  gt(exclude_summary_all) %>% 
  cols_align('center') %>% 
  fmt_number('n', decimals = 0) %>% 
  fmt_percent('pct', decimals = 1) %>% 
  cols_label('exclusion' = 'Exclusion',
             'n' = '# of Patients',
             'pct' = '%') %>% 
  tab_header(title = md('**Table of Exclusions for Patients Undergoing Bariatric Surgery**'),
             subtitle = md('**Among 45,956 Patients in DURABLE Database**'))
cluster_gtsave(gt_exclude_all, 'figures/tables/exclusion_summary_table_DURABLE.png', vheight = 100, vwidth = 700)
