library(tidyverse)
library(arrow)
library(glue)
library(lubridate)
library(splines)
source('scripts/helpers.R')


### Set Seed for boot strapping
args <- commandArgs(trailingOnly = T)
boot_id <- as.numeric(args[1])
n_boot <- 1000
set.seed(73097)
seeds <- sample(1:1000000, n_boot, replace = F)
set.seed(seeds[boot_id])

### Directory where EHR data is stored
ehr_dir <- '/n/haneuse_ehr_l3/V1.0'
data_dir <- '/n/haneuse_ehr_l3/V1.0/clean_datasets'

### Load in Pre-Saved Data
max_follow_up <- 84
n_trials <- 96
trial_range <- 1:84
# df_trials <- map_dfr(1:n_trials, ~read_parquet(glue('{data_dir}/bariatric_tte/trial_{.x}_tte1.parquet')))
df_trials <- read_parquet(glue('{data_dir}/bariatric_tte/all_trials_combined.parquet'))

criteria <- c('elig_age', 'elig_diabetes', 'elig_surgery', 'elig_bmi_missing', 'elig_bmi_range', 'elig_cvd')
pre_op <- c('elig_smoking', 'elig_bmi_change')
exclusions <- c('elig_exclusions', 'elig_pregnancy')

formulae <- 
  list('Linear Main Effects' = 
         surgery ~ 
         ### Categorical
         race_black + gender + smoking_status + 
         hypertension + dyslipidemia + past_1yr_insulin +
         
         ### Continuous
         baseline_age + baseline_bmi + baseline_a1c + 
         bmi_change_1yr + a1c_change_1yr + trial_id,
       
       
       'Madenci Proxy (NCS on Continuous Covariates, No Smoking or 1-Year A1c Change)' = 
         surgery ~ 
         ### Categorical
         race_black + gender + past_1yr_insulin + 
         hypertension + dyslipidemia + 
         
         ### Continuous
         ns(trial_id, knots = quantile(trial_id, 0.5), Boundary.knots = quantile(trial_id, c(0.1, 0.9))) + 
         ns(baseline_age, knots = quantile(baseline_age, 0.5), Boundary.knots = quantile(baseline_age, c(0.1, 0.9))) +
         ns(baseline_bmi, knots = quantile(baseline_bmi, 0.5), Boundary.knots = quantile(baseline_bmi, c(0.1, 0.9))) + 
         ns(baseline_a1c, knots = quantile(baseline_a1c, 0.5), Boundary.knots = quantile(baseline_a1c, c(0.1, 0.9))) +
         ns(bmi_change_1yr, knots = quantile(bmi_change_1yr, 0.5), Boundary.knots = quantile(bmi_change_1yr, c(0.1, 0.9))),
       
       'Best BIC Model (Splines on BMI/Age, No A1c 1yr Change)' = 
         surgery ~ 
         ### Categorical
         race_black + gender + smoking_status + 
         hypertension + dyslipidemia + past_1yr_insulin +
         
         ### Continuous
         baseline_a1c + trial_id + bmi_change_1yr +
         
         ns(baseline_bmi, knots = quantile(baseline_bmi, 0.5), Boundary.knots = quantile(baseline_bmi, c(0.1, 0.9))) +
         ns(baseline_age, knots = 33)
  )


### Loop Over Various Elig Criteria to run analysis
df_ci_all <- NULL
for(elig in 1:3) {
  cat('Eligibility Criteria', elig, '\n') 
  if(elig == 1) {
    elig_criteria <- criteria 
  } else if(elig == 2) {
    elig_criteria <- c(criteria, pre_op)
  } else if(elig == 3) {
    elig_criteria <- c(criteria, pre_op, exclusions)
  }
  
  df_analysis <- 
    df_trials %>% 
    filter(trial_id %in% trial_range) %>% 
    filter_at(.vars = c(elig_criteria), ~.x) %>% 
    filter(!is.na(baseline_a1c)) ### Remove Missing Baseline a1c for now
  
  for(i in 1:length(formulae)) {
    ### Bootstrap at (Eligible) Subject Level
    subject_ids <- unique(df_analysis$subject_id)
    ix_subj <- sample(1:length(subject_ids), size = length(subject_ids), replace = T)
    df_subj <- tibble('subject_id' = subject_ids[ix_subj])
    
    df_sample <- 
      df_analysis %>% 
      inner_join(df_subj, by = 'subject_id', relationship = 'many-to-many') %>% 
      group_by(subject_id, trial_id) %>% 
      mutate('replicate_id' = 1:n()) %>% 
      ungroup()
    
    
    ### IPW Model for Baseline Confounding
    cat('IPW Model Eligibility Criteria', elig, names(formulae)[i], '\n') 
    ipw_model <-
      speedglm::speedglm(formulae[[i]],
                         family = binomial(link = 'logit'),
                         data = df_sample)
    
    df_sample$p_surgery <- predict(ipw_model, newdata = df_sample, type = 'response')
    
    df_sample <- 
      df_sample %>% 
      mutate('w_treatment' = case_when(surgery == 0 ~ 1/(1-p_surgery),
                                       surgery == 1 ~ 1/p_surgery),
             'sw_treatment' = case_when(surgery == 0 ~ mean(1-surgery)/(1-p_surgery),
                                        surgery == 1 ~ mean(surgery)/(p_surgery))) %>% 
      ### Trim Weights 
      group_by(surgery) %>% 
      mutate('w_treatment' = winsorize(w_treatment, q = c(0, 0.99)),
             'sw_treatment' = winsorize(sw_treatment, q = c(0, 0.99))) %>% 
      ungroup()
    
    
    ### Expanded Dataset for Outcome Model
    ### Compute Time to Event Outcpmes
    df_times <-   
      df_sample %>% 
      mutate('cvd_time' = pmax(1, ceiling(cvd_time))) %>% 
      select(trial_id, subject_id, replicate_id, cvd_time, cvd_event) %>% 
      group_by(trial_id, subject_id, replicate_id) %>% 
      reframe('time' = 1:cvd_time,
              'outcome' = as.numeric(cvd_event == 1 & time == cvd_time)) %>% 
      ungroup() 
    
    cat('Outcome Model Eligibility Criteria', elig, names(formulae)[i], '\n') 
    df_expanded <- 
      df_sample %>% 
      select(subject_id, trial_id, replicate_id, surgery, w_treatment, sw_treatment) %>% 
      inner_join(df_times, by = c('subject_id', 'trial_id', 'replicate_id')) 
    
    df_model <- 
      df_expanded %>% 
      group_by(time, surgery) %>% 
      summarise('n' = n(),
                'n_weighted' = sum(sw_treatment),
                'n_outcomes' = sum(outcome),
                'n_outcomes_weighted' = sum(outcome * sw_treatment)) 
    
    
    outcome_model <- 
      speedglm::speedglm(cbind(n_outcomes_weighted, n_weighted - n_outcomes_weighted) ~ surgery * ns(time, knots = 12 * c(2, 3.5, 5), Boundary.knots = 12 * c(0.5, 6)), 
                         family =  binomial(link = 'logit'),
                         data = df_model)
    
    
    
    ### Get CI   
    cat('Compute Cumulative Incidence', elig, '\n') 
    df_ci <- 
      crossing('surgery' = c(0,1),
               'time' = 0:max_follow_up) %>% 
      mutate('p_event' = predict(outcome_model, newdata = ., type = 'response')) %>% 
      mutate('p_event' = ifelse(time == 0, 0, p_event)) %>% 
      arrange(time) %>% 
      group_by(surgery) %>% 
      mutate('cum_surv' = cumprod(1-p_event),
             'cum_incidence' = 1 - cum_surv) %>% 
      mutate('elig_criteria' = case_when(elig == 1 ~ 'No Pre-Op Restrictions',
                                         elig == 2 ~ 'Pre-Op Restrictions',
                                         elig == 3 ~ 'Pre-Op Restrictions + Additional Exclusions'),
             'ipw_formula' = names(formulae)[i],
             'outcome_model' = '(Madenci et al., 2024) Knot Placement',
             'outcome_formula' = 'CVD ~ surgery * spline(time, knots = (2, 3.5, 5), boundary.knots = (0.5, 6))') %>% 
      ungroup() %>% 
      mutate('boot_id' = boot_id)
    
    df_ci_all <- 
      df_ci_all %>%
      bind_rows(df_ci)
    
  }
}

if(!dir.exists(glue('data/bootstraps/model_comp'))) {
  dir.create(glue('data/bootstraps/model_comp'))
}
write_parquet(df_ci_all, glue('data/bootstraps/model_comp/cvd_ci_{boot_id}.parquet'))