library(tidyverse)
library(arrow)
library(glue)
library(lubridate)
library(splines)
source('scripts/helpers.R')

### Directory where EHR data is stored
ehr_dir <- '/n/haneuse_ehr_l3/V1.0'
data_dir <- '/n/haneuse_ehr_l3/V1.0/clean_datasets'

### Load in Pre-Saved Data
max_follow_up <- 84
n_trials <- 96
trial_range <- 1:84
# trial_range <- 
#   list('Trials 2005 - 2012' = 1:96,
#        'Trials 2007 - 2012' = 25:96,
#        'Trials 2005 - 2011' = 1:84)


df_trials <- read_parquet(glue('{data_dir}/bariatric_tte/all_trials_combined.parquet'))

### Various eligibility criteria
criteria <- c('elig_age', 'elig_diabetes', 'elig_surgery', 'elig_bmi_missing', 'elig_bmi_range', 'elig_cvd', 'elig_pregnancy', 'elig_cancer')
pre_op <- c('elig_smoking', 'elig_bmi_change')
exclusions <- c('elig_exclusions')

formulae <- 
  list('Linear Main Effects' = 
         surgery ~ 
         ### Categorical
         race_black + gender + smoking_status + 
         hypertension + dyslipidemia + past_1yr_insulin +
         
         ### Continuous
         baseline_age + baseline_bmi + baseline_a1c + 
         bmi_change_1yr + a1c_change_1yr + trial_id,
       
       
       'Madenci Proxy' = 
         surgery ~ 
         ### Categorical
         race_black + gender + past_1yr_insulin + 
         hypertension + dyslipidemia + smoking_status + 
         
         ### Continuous
         ns(trial_id, knots = quantile(trial_id, 0.5), Boundary.knots = quantile(trial_id, c(0.1, 0.9))) + 
         ns(baseline_age, knots = quantile(baseline_age, 0.5), Boundary.knots = quantile(baseline_age, c(0.1, 0.9))) +
         ns(baseline_bmi, knots = quantile(baseline_bmi, 0.5), Boundary.knots = quantile(baseline_bmi, c(0.1, 0.9))) + 
         ns(baseline_a1c, knots = quantile(baseline_a1c, 0.5), Boundary.knots = quantile(baseline_a1c, c(0.1, 0.9))) +
         ns(bmi_change_1yr, knots = quantile(bmi_change_1yr, 0.5), Boundary.knots = quantile(bmi_change_1yr, c(0.1, 0.9))) +
         ns(a1c_change_1yr, knots = quantile(a1c_change_1yr, 0.5), Boundary.knots = quantile(a1c_change_1yr, c(0.1, 0.9)))
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
  
  ### Filter to Eligible Subjects
  df_analysis <- 
    df_trials %>% 
    filter(trial_id %in% trial_range) %>% 
    filter_at(.vars = c(elig_criteria), ~.x) %>% 
    filter(!is.na(baseline_a1c)) ### Remove Missing Baseline a1c for now
  
  for(i in 1:length(formulae)) {
    ### IPW Model for Baseline Confounding
    cat('IPW Model Eligibility Criteria', elig, names(formulae)[i], '\n') 
    ipw_model <-
      speedglm::speedglm(formulae[[i]] ,
                         family = binomial(link = 'logit'),
                         data = df_analysis)
    
    df_analysis$p_surgery <- predict(ipw_model, newdata = df_analysis, type = 'response')
    
    df_analysis <- 
      df_analysis %>% 
      mutate('w_treatment' = case_when(surgery == 0 ~ 1/(1-p_surgery),
                                       surgery == 1 ~ 1/p_surgery),
             'sw_treatment' = case_when(surgery == 0 ~ mean(1-surgery)/(1-p_surgery),
                                        surgery == 1 ~ mean(surgery)/(p_surgery))) %>% 
      ### Trim Weights 
      group_by(surgery) %>%
      mutate('w_treatment' = winsorize(w_treatment, q = c(0, 0.99)),
             'sw_treatment' = winsorize(sw_treatment, q = c(0, 0.99))) %>%
      ungroup()
    
    print(max(df_analysis$sw_treatment))
    
    
    ### Expanded Dataset for Outcome Model
    ### Compute Time to Event Outcpmes
    df_times <-   
      df_analysis %>% 
      mutate('cvd_time' = pmax(1, ceiling(cvd_time))) %>% 
      select(trial_id, subject_id, cvd_time, cvd_event) %>% 
      group_by(trial_id, subject_id) %>% 
      reframe('time' = 1:cvd_time,
              'outcome' = as.numeric(cvd_event == 1 & time == cvd_time)) %>% 
      ungroup() 
    
    cat('Outcome Model Eligibility Criteria', elig, names(formulae)[i], '\n') 
    df_expanded <- 
      df_analysis %>% 
      select(subject_id, trial_id, surgery, w_treatment, sw_treatment) %>% 
      inner_join(df_times, by = c('subject_id', 'trial_id')) 
    
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
    cat('Compute Cumulative Incidence\n') 
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
             'formulae' = names(formulae)[i],
             'outcome_model' = '(Madenci et al., 2024) Knot Placement',
             'formula' = 'CVD ~ surgery * spline(time, knots = (2, 3.5, 5), boundary.knots = (0.5, 6))') %>% 
      ungroup()
    
    df_ci_all <- 
      df_ci_all %>%
      bind_rows(df_ci)
  }
}

write_parquet(df_ci_all, 'data/results/cvd_cumulative_incidence_by_formula.parquet')
