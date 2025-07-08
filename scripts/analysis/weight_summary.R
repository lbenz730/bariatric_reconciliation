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

trial_range <- 1:84

df_trials <- read_parquet(glue('{data_dir}/bariatric_tte/all_trials_combined.parquet'))

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
weight_summary <- NULL
all_weights <- NULL
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
      mutate('w_treatment_trunc' = winsorize(w_treatment, q = c(0, 0.99)),
             'sw_treatment_trunc' = winsorize(sw_treatment, q = c(0, 0.99))) %>%
      ungroup()
    
    ### Weight Summary
    df_summary <- 
      bind_rows(
        df_analysis %>% 
          group_by(surgery = as.character(surgery)) %>% 
          summarise('q_90' = quantile(sw_treatment_trunc, 0.9),
                    'q_95' = quantile(sw_treatment_trunc, 0.95),
                    'q_99' = quantile(sw_treatment_trunc, 0.99),
                    'min' = min(sw_treatment_trunc),
                    'max' = max(sw_treatment_trunc),
                    'max_untruncated' = max(sw_treatment)),
        
        df_analysis %>% 
          summarise('surgery' = 'all', 
                    'q_90' = quantile(sw_treatment_trunc, 0.9),
                    'q_95' = quantile(sw_treatment_trunc, 0.95),
                    'q_99' = quantile(sw_treatment_trunc, 0.99),
                    'min' = min(sw_treatment_trunc),
                    'max' = max(sw_treatment_trunc),
                    'max_untruncated' = max(sw_treatment))
      ) %>% 
      mutate('elig_criteria' = case_when(elig == 1 ~ 'Trial #1A', 
                                         elig == 2 ~ 'Trial #1B',
                                         elig == 3 ~ 'Trial #1C'),
             'ipw_formula' = names(formulae)[[i]])
    
    weight_summary <- 
      weight_summary %>% 
      bind_rows(df_summary)
    
    df_weights <- 
      df_analysis %>% 
      select(subject_id, trial_id, surgery, sw_treatment, sw_treatment_trunc) %>% 
      
      mutate('elig_criteria' = case_when(elig == 1 ~ 'Trial #1A', 
                                         elig == 2 ~ 'Trial #1B',
                                         elig == 3 ~ 'Trial #1C'),
             'ipw_formula' = names(formulae)[[i]])
    
    all_weights <- 
      all_weights %>% 
      bind_rows(df_weights)
    
    
  }
}

all_weights <- 
  all_weights %>%
  mutate('ipw_formula' = case_when(ipw_formula == 'Linear Main Effects' ~ 'Propensity Score Model #1',
                                   ipw_formula == 'Madenci Proxy' ~ 'Propensity Score Model #2'))


### Plots
ggplot(data = all_weights, aes(x = sw_treatment_trunc, y = ..density..)) + 
  facet_wrap(ifelse(surgery == 1, 'Surgery', 'No Surgery') ~ elig_criteria, scales = 'free') + 
  geom_histogram(aes(fill = ipw_formula), position = 'identity', color = 'black', alpha = 0.4) + 
  labs(x = 'Inverse Probability of Treatment Weight',
       y = 'Density',
       fill = '',
       title = 'Distribution of Inverse Probability of Treatment Weights',
       subtitle = 'Stabilized, Truncated at 99th Percentile')
ggsave('figures/results/iptw_distribution.png', height = 9/1.2, width = 16/1.2)


weight_summary <- 
  weight_summary %>% 
  mutate('ipw_formula' = case_when(ipw_formula == 'Linear Main Effects' ~ 'Propensity Score Model #1',
                                   ipw_formula == 'Madenci Proxy' ~ 'Propensity Score Model #2'))
       
write_csv(weight_summary, 'figures/tables/iptw_summary.csv')
