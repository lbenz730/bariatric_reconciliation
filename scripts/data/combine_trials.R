library(tidyverse)
library(arrow)
library(glue)

### Directory where EHR data is stored
ehr_dir <- '/n/haneuse_ehr_l3/V1.0'
data_dir <- '/n/haneuse_ehr_l3/V1.0/clean_datasets'


n_trials <- 96
df_trials <- map_dfr(1:n_trials, ~read_parquet(glue('{data_dir}/bariatric_tte/trial_{.x}_tte1.parquet')))
write_parquet(df_trials, glue('{data_dir}/bariatric_tte/all_trials_combined.parquet'))

n_trials <- 84
durable_trials <- map_dfr(1:n_trials, ~read_parquet(glue('{data_dir}/bariatric_tte/trial_{.x}_DURABLE.parquet')))
write_parquet(durable_trials, glue('{data_dir}/bariatric_tte/DURABLE_trials_combined.parquet'))
