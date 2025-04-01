library(tidyverse)
library(arrow)
library(glue)
library(haven)
library(lubridate)

### Directory where EHR data is stored
ehr_dir <- '/n/haneuse_ehr_l3/V1.0'
data_dir <- '/n/haneuse_ehr_l3/V1.0/clean_datasets'
weights <- read_parquet(glue('{data_dir}/all_weights.parquet')) ### Cleaned Weights w/ outliers removed

### Further Weight Outliers (not model based) 
weights <-
  weights %>% 
  group_by(subject_id) %>%
  mutate('bmi_lag1' = lag(bmi),
         'bmi_lead1' = lead(bmi),
         'max_bmi' = max(bmi),
         'min_bmi' = min(bmi),
         'second_max_bmi' = bmi[kit::topn(bmi, n = 2, decreasing = T, hasna = T)][2],
         'second_min_bmi' = bmi[kit::topn(bmi, n = 2, decreasing = F, hasna = T)][2]) %>%
  ungroup() %>%
  mutate('outlier_skip' = abs(bmi - bmi_lag1) >= 10 & abs(bmi_lead1 - bmi_lag1) <= 5) %>%
  mutate('outlier_endpoint' =
           (!is.na(second_max_bmi) & (bmi == max_bmi & (max_bmi - second_max_bmi) >= 10)) |
           (!is.na(second_min_bmi) & (bmi == min_bmi & (min_bmi - second_max_bmi) <= -10))) %>%
  filter(is.na(outlier_skip) | !outlier_skip) %>%
  filter(!outlier_endpoint) %>%
  select(-bmi_lag1, -bmi_lead1, -max_bmi, -min_bmi, -contains('second'))
write_parquet(weights, glue('{data_dir}/all_weights_further_cleaned.parquet'))
