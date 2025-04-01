library(tidyverse)
library(arrow)
library(glue)
library(lubridate)
library(splines)
library(gt)
source('scripts/helpers.R')

### Directory where EHR data is stored
ehr_dir <- '/n/haneuse_ehr_l3/V1.0'
data_dir <- '/n/haneuse_ehr_l3/V1.0/clean_datasets'

### Load in Pre-Saved Data
max_follow_up <- 84
n_trials <- 96
trial_range <- 
  list('Trials 2005 - 2012' = 1:96,
       'Trials 2007 - 2012' = 25:96,
       'Trials 2005 - 2011' = 1:84)

df_trials <- map_dfr(1:n_trials, ~read_parquet(glue('{data_dir}/bariatric_tte/trial_{.x}_tte1.parquet')))

### Various eligibility criteria
criteria <- c('elig_age', 'elig_diabetes', 'elig_surgery', 'elig_bmi_missing', 'elig_bmi_range', 'elig_cvd')
pre_op <- c('elig_smoking', 'elig_bmi_change')
exclusions <- c('elig_exclusions', 'elig_pregnancy')

elig_criteria <- c(criteria)

df_analysis <- 
  df_trials %>% 
  filter(trial_id <= 84) %>% 
  filter_at(.vars = c(elig_criteria), ~.x) %>% 
  filter(!is.na(baseline_a1c)) ### Remove Missing Baseline a1c for now

df_continuous <- 
  df_analysis %>% 
  select(surgery, trial_id, baseline_age, baseline_a1c, baseline_bmi, bmi_change_1yr, a1c_change_1yr) %>% 
  pivot_longer(-surgery, 
               names_to = 'variable',
               values_to = 'value') %>% 
  mutate('bucket' = case_when(variable == 'trial_id' ~ value,
                              variable == 'baseline_age' ~ round(value),
                              variable == 'baseline_bmi' ~ round(value),
                              variable == 'baseline_a1c' ~ plyr::round_any(value, 0.5),
                              variable == 'bmi_change_1yr' ~ round(value),
                              variable == 'a1c_change_1yr' ~ round(value))) %>% 
  group_by(variable, bucket) %>% 
  summarise('pct_surgery' = mean(surgery),
            'n_obs' = n()) %>% 
  ungroup() %>% 
  mutate('variable' = case_when(variable == 'trial_id' ~ 'Calendar Time (Trial #)',
                                variable == 'baseline_age' ~ 'Baseline Age',
                                variable == 'baseline_bmi' ~ 'Baseline BMI',
                                variable == 'baseline_a1c' ~ 'Baseline A1c',
                                variable == 'bmi_change_1yr' ~ '1 Year BMI Change',
                                variable == 'a1c_change_1yr' ~ '1 Year A1c Change'))

### BIC analysis
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
       
       'NCS on Calendar Time (Single Knot at Median)' = 
         surgery ~ 
         ### Categorical
         race_black + gender + smoking_status + 
         hypertension + dyslipidemia + past_1yr_insulin +
         
         ### Continuous
         baseline_age + baseline_bmi + baseline_a1c + 
         bmi_change_1yr + a1c_change_1yr +
         
         ns(trial_id, knots = quantile(trial_id, 0.5), Boundary.knots = quantile(trial_id, c(0.1, 0.9))),
       
       'NCS on Calendar Time (Single User Chosen Knot)' = 
         surgery ~ 
         ### Categorical
         race_black + gender + smoking_status + 
         hypertension + dyslipidemia + past_1yr_insulin +
         
         ### Continuous
         baseline_age + baseline_bmi + baseline_a1c + 
         bmi_change_1yr + a1c_change_1yr +
         
         ns(trial_id, knots = c(30)),
       
       'NCS on Baseline Age (Single Knot at Median)' = 
         surgery ~ 
         ### Categorical
         race_black + gender + smoking_status + 
         hypertension + dyslipidemia + past_1yr_insulin +
         
         ### Continuous
         baseline_bmi + baseline_a1c + trial_id + 
         bmi_change_1yr + a1c_change_1yr +
         
         ns(baseline_age, knots = quantile(baseline_age, 0.5), Boundary.knots = quantile(baseline_age, c(0.1, 0.9))),
       
       'NCS on Baseline Age (Single User Chosen Knot)' = 
         surgery ~ 
         ### Categorical
         race_black + gender + smoking_status + 
         hypertension + dyslipidemia + past_1yr_insulin +
         
         ### Continuous
         baseline_bmi + baseline_a1c + trial_id + 
         bmi_change_1yr + a1c_change_1yr +
         
         ns(baseline_age, knots = 33),
       
       'NCS on Baseline BMI (Single Knot at Median)' = 
         surgery ~ 
         ### Categorical
         race_black + gender + smoking_status + 
         hypertension + dyslipidemia + past_1yr_insulin +
         
         ### Continuous
         baseline_age + trial_id + baseline_a1c + 
         bmi_change_1yr + a1c_change_1yr +
         
         ns(baseline_bmi, knots = quantile(baseline_bmi, 0.5), Boundary.knots = quantile(baseline_bmi, c(0.1, 0.9))),
       
       'NCS on Baseline BMI (Single User Chosen Knot)' = 
         surgery ~ 
         ### Categorical
         race_black + gender + smoking_status + 
         hypertension + dyslipidemia + past_1yr_insulin +
         
         ### Continuous
         baseline_age + trial_id + baseline_a1c + 
         bmi_change_1yr + a1c_change_1yr +
         
         ns(baseline_bmi, knots = 35),
       
       'NCS on Baseline A1c (Single Knot at Median)' = 
         surgery ~ 
         ### Categorical
         race_black + gender + smoking_status + 
         hypertension + dyslipidemia + past_1yr_insulin +
         
         ### Continuous
         baseline_age + baseline_bmi + trial_id + 
         bmi_change_1yr + a1c_change_1yr +
         
         ns(baseline_a1c, knots = quantile(baseline_a1c, 0.5), Boundary.knots = quantile(baseline_a1c, c(0.1, 0.9))),
       
       
       'NCS on BMI Change 1 Year (Single Knot at Median)' = 
         surgery ~ 
         ### Categorical
         race_black + gender + smoking_status + 
         hypertension + dyslipidemia + past_1yr_insulin +
         
         ### Continuous
         baseline_age + baseline_bmi + trial_id + 
         baseline_a1c + a1c_change_1yr +
         
         ns(bmi_change_1yr, knots = quantile(bmi_change_1yr, 0.5), Boundary.knots = quantile(bmi_change_1yr, c(0.1, 0.9))),
       
       
       'NCS on A1c Change 1 Year (Single Knot at Median)' = 
         surgery ~ 
         ### Categorical
         race_black + gender + smoking_status + 
         hypertension + dyslipidemia + past_1yr_insulin +
         
         ### Continuous
         baseline_age + baseline_bmi + trial_id + 
         baseline_a1c + bmi_change_1yr +
         
         ns(a1c_change_1yr, knots = quantile(a1c_change_1yr, 0.5), Boundary.knots = quantile(a1c_change_1yr, c(0.1, 0.9))),
       
       'NCS on A1c Change 1 Year (Single User Chosen Knot)' = 
         surgery ~ 
         ### Categorical
         race_black + gender + smoking_status + 
         hypertension + dyslipidemia + past_1yr_insulin +
         
         ### Continuous
         baseline_age + baseline_bmi + trial_id + 
         baseline_a1c + bmi_change_1yr +
         
         ns(a1c_change_1yr, knots = c(-3)),
       
       'NCS on Baseline BMI and Age (Best Single Knots)' = 
         surgery ~ 
         ### Categorical
         race_black + gender + smoking_status + 
         hypertension + dyslipidemia + past_1yr_insulin +
         
         ### Continuous
         baseline_a1c + trial_id + bmi_change_1yr +
         
         ns(baseline_bmi, knots = quantile(baseline_bmi, 0.5), Boundary.knots = quantile(baseline_bmi, c(0.1, 0.9))) +
         ns(baseline_age, knots = 33),
       
       'NCS on Baseline BMI and Calendar Time, No A1c Change (Best Single Knots)' = 
         surgery ~ 
         ### Categorical
         race_black + gender + smoking_status + 
         hypertension + dyslipidemia + past_1yr_insulin +
         
         ### Continuous
         baseline_a1c + baseline_age + bmi_change_1yr +
         
         ns(baseline_bmi, knots = quantile(baseline_bmi, 0.5), Boundary.knots = quantile(baseline_bmi, c(0.1, 0.9))) +
         ns(trial_id, knots = quantile(trial_id, 0.5), Boundary.knots = quantile(trial_id, c(0.1, 0.9))),
       
       'NCS on Baseline Age and Calendar Time, No A1c Change  (Best Single Knots)' = 
         surgery ~ 
         ### Categorical
         race_black + gender + smoking_status + 
         hypertension + dyslipidemia + past_1yr_insulin +
         
         ### Continuous
         baseline_a1c + baseline_bmi + bmi_change_1yr +
         
         ns(trial_id, knots = quantile(trial_id, 0.5), Boundary.knots = quantile(trial_id, c(0.1, 0.9))) +
         ns(baseline_age, knots = 33)
       
       # 'NCS on Baseline BMI, Age and Calendar Time, No A1c Change (Best Single Knots)' = 
       #   surgery ~ 
       #   ### Categorical
       #   race_black + gender + smoking_status + 
       #   hypertension + dyslipidemia + past_1yr_insulin +
       #   
       #   ### Continuous
       #   baseline_a1c + bmi_change_1yr +
       #   
       #   ns(trial_id, knots = quantile(trial_id, 0.5), Boundary.knots = quantile(trial_id, c(0.1, 0.9))) +
       #   ns(baseline_bmi, knots = quantile(baseline_bmi, 0.5), Boundary.knots = quantile(baseline_bmi, c(0.1, 0.9))) +
       #   ns(baseline_age, knots = 33)
       
  )


evaluate_model <- function(ipw_formula) {
  
  ipw_model <-
    speedglm::speedglm(ipw_formula,
                       family = binomial(link = 'logit'),
                       data = df_analysis)
  
  
  df_analysis_ <- 
    df_analysis %>% 
    mutate('p_surgery' = predict(ipw_model, newdata = ., type = 'response')) %>% 
    mutate('w_treatment' = case_when(surgery == 0 ~ 1/(1-p_surgery),
                                     surgery == 1 ~ 1/p_surgery),
           'sw_treatment' = case_when(surgery == 0 ~ mean(1-surgery)/(1-p_surgery),
                                      surgery == 1 ~ mean(surgery)/(p_surgery))) %>% 
    ### Trim Weights 
    group_by(surgery) %>% 
    mutate('w_treatment' = winsorize(w_treatment, q = c(0, 0.99)),
           'sw_treatment' = winsorize(sw_treatment, q = c(0, 0.99))) %>% 
    ungroup()
  
  df_stats <-
    tibble('bic' = BIC(ipw_model),
           'log_lik' = ipw_model$logLik,
           'n_covariates' = length(ipw_model$coefficients),
           'max_SW' = max(df_analysis_$sw_treatment))
  
  return(df_stats)
  
}


### Fit all models and show model summary
df_eval <- 
  map_dfr(formulae, evaluate_model) %>% 
  mutate('model' = names(formulae)) %>% 
  select(model, everything()) 

write_csv(df_eval, 'data/results/spline_BIC.csv') 

gt_bic <- 
  gt(df_eval %>% arrange(bic)) %>% 
  cols_align('center') %>% 
  fmt_number(c('bic', 'log_lik'), decimals = 1, use_seps =  F) %>% 
  fmt_number(c('max_SW'), decimals = 2) %>% 
  cols_label('model' = 'Model', 
             'bic' = 'BIC',
             'log_lik' = 'Log-Likelihood',
             'n_covariates' = 'Design Matrix Dim.',
             'max_SW' = 'Max SW') %>% 
  tab_header(title = md('**Model Selection Table (BIC) for IPW of Treatment (Surgery) Model**')) %>% 
  tab_options(column_labels.font.size = 20,
              heading.title.font.size = 36,
              heading.subtitle.font.size = 30,
              heading.title.font.weight = 'bold',
              heading.subtitle.font.weight = 'bold',
              column_labels.font.weight = 'bold',
              row_group.font.weight = 'bold') 

gt_bic
cluster_gtsave(gt_bic, 'figures/tables/BIC_table.png')