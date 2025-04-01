library(tidyverse)
library(arrow)
library(glue)
library(haven)
library(lubridate)
source('scripts/helpers.R')


### Directory where EHR data is stored
ehr_dir <- '/n/haneuse_ehr_l3/V1.0'
data_dir <- '/n/haneuse_ehr_l3/V1.0/clean_datasets'

### Load in Weight Data
df_subjects <- read_parquet(glue('{data_dir}/microvascular_tte/subjects.parquet'))
df_enrollment <- read_parquet(glue('{data_dir}/microvascular_tte/enrollment.parquet'))
case_px <- read_sas(glue('{ehr_dir}/_data_for_eric_cases/durable_cases.sas7bdat'))
weights <- read_parquet(glue('{data_dir}/all_weights_further_cleaned.parquet')) ### Cleaned Weights w/ outliers removed

### Drop measures on the same date
weights <-
  weights %>% 
  distinct(subject_id, measure_date, .keep_all = T)


### Filter to People who underwent BS and were enrolled for 1 year before 
df_surg1 <- 
  df_subjects %>% 
  filter(!is.na(index_date)) %>% 
  mutate('index_date_minus1yr' = index_date %m-% years(1)) %>% 
  inner_join(df_enrollment, by = 'subject_id') %>% 
  filter(index_date >= enr_start %m+% years(1), enr_end >= index_date) %>% 
  left_join(select(weights, subject_id, bmi, measure_date),
            join_by(subject_id, closest(index_date >= measure_date))) %>% 
  left_join(select(weights, subject_id, bmi, measure_date),
            join_by(subject_id, closest(index_date_minus1yr <= measure_date)),
            suffix = c('', '_1yr')) %>% 
  filter(!is.na(bmi), measure_date >= index_date_minus1yr) %>% 
  mutate('bmi_diff' = bmi - bmi_1yr,
         'measure_diff' = as.numeric(measure_date - measure_date_1yr)) %>% 
  mutate('bmi_diff_cat' = case_when(bmi_diff > 2 ~ 'Significant BMI Gain (Exclude)',
                                    bmi_diff > 0 ~ 'Minimal BMI Gain (Include)',
                                    T ~ 'BMI Loss (Include)'))

### Recreate this chart for every measure, for a1c
ggplot(df_surg1, aes(x = measure_diff, y = bmi_diff)) + 
  geom_point(alpha = 0.5, aes(color = bmi_diff_cat)) + 
  labs(x = 'Days Before Measurement Used as Baseline BMI',
       y = 'BMI Change From Date Until Index Date',
       # title = 'BMI Change From Measure Closest to 1 Year Before Index Date',
       # subtitle = 'Among Patients Undergoing Bariatric Surgery',
       color = '')
ggsave('figures/results/bmi_change.png', height = 9/1.2, width = 16/1.2)


### Filter to People who underwent BS and were enrolled for 1 year before 
df_surg2 <- 
  df_subjects %>% 
  filter(!is.na(index_date)) %>% 
  mutate('index_date_minus1yr' = index_date %m-% years(1)) %>% 
  inner_join(df_enrollment, by = 'subject_id') %>% 
  filter(index_date >= enr_start %m+% years(1), enr_end >= index_date) %>% 
  left_join(select(weights, subject_id, bmi, measure_date),
            join_by(subject_id, closest(index_date >= measure_date))) %>% ### Closest Measure Before Index Date w/in 
  mutate('measure_date_minus1yr' = measure_date %m-% years(1)) %>% 
  left_join(select(weights, subject_id, bmi, measure_date),
            join_by(subject_id, closest(measure_date_minus1yr <= measure_date)),
            suffix = c('', '_1yr')) %>% 
  filter(!is.na(bmi), measure_date >= index_date_minus1yr) %>% 
  mutate('bmi_diff' = bmi - bmi_1yr,
         'measure_diff' = as.numeric(measure_date - measure_date_1yr)) %>% 
  mutate('bmi_diff_cat' = case_when(bmi_diff > 2 ~ 'Significant BMI Gain (Exclude)',
                                    bmi_diff > 0 ~ 'Minimal BMI Gain (Include)',
                                    T ~ 'BMI Loss (Include)'))

### Recreate this chart for every measure, for a1c
ggplot(df_surg2, aes(x = measure_diff, y = bmi_diff)) + 
  geom_point(alpha = 0.5, aes(color = bmi_diff_cat)) + 
  labs(x = 'Days Before Measurement Used as Baseline BMI',
       y = 'BMI Change From Date Until Index Date',
       title = 'BMI Change From Measure Closest to 1 Year Before Date of BMI Used at Baseline',
       subtitle = 'Among Patients Undergoing Bariatric Surgery',
       color = '')
ggsave('figures/eda/bmi_change_alt.png', height = 9/1.2, width = 16/1.2)