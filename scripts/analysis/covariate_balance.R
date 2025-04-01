library(tidyverse)
library(arrow)
library(glue)
library(lubridate)
library(splines)
library(cobalt)
source('scripts/helpers.R')

### Directory where EHR data is stored
ehr_dir <- '/n/haneuse_ehr_l3/V1.0'
data_dir <- '/n/haneuse_ehr_l3/V1.0/clean_datasets'

### Load in Pre-Saved Data

trial_range <- 1:84

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
df_balance_all <- NULL
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
  
  covariates <- 
    df_analysis %>% 
    select(baseline_age, gender, race_black, 
           smoking_status, baseline_bmi, baseline_a1c, 
           bmi_change_1yr, a1c_change_1yr, 
           past_1yr_insulin, hypertension, dyslipidemia) 
  
  
  ### Unadjusted balance
  unadj_balance <- 
    bal.tab(covariates, 
            treat = 'surgery', 
            data = df_analysis, 
            s.d.denom = "pooled")
  
  df_balance_all <- 
    df_balance_all %>% 
    bind_rows(
      tibble('covariate' = rownames(unadj_balance$Balance),
             'smd' = unadj_balance$Balance$Diff.Un,
             'elig_criteria' = elig, 
             'ipw' = 'Unadjusted')
    )
  
  
  
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
    
    ### Adjusted Balance
    adj_balance <- 
      bal.tab(covariates, 
              treat = 'surgery', 
              data = df_analysis, 
              weights = df_analysis$sw_treatment,
              s.d.denom = "pooled")
    
    df_balance_all <- 
      df_balance_all %>% 
      bind_rows(
        tibble('covariate' = rownames(adj_balance$Balance),
               'smd' = adj_balance$Balance$Diff.Adj,
               'elig_criteria' = elig, 
               'ipw' = paste('Propensity Score Model #', i))
      )
    
  }
}

### Save Results
write_csv(df_balance_all, 'data/results/covariate_balance.csv') 


### Plot Results
df_balance_all <- read_csv('data/results/covariate_balance.csv')
df_plot <- 
  df_balance_all %>% 
  mutate('elig_criteria' = case_when(elig_criteria == 1 ~ 'Eligibility Set #1', 
                                     elig_criteria == 2 ~ 'Eligibility Sets #1 + #2',
                                     elig_criteria == 3 ~ 'Eligibility Sets #1 + #2 + #3'),
         'covariate' = case_when(covariate == 'a1c_change_1yr' ~ '1-Year A1c Change',
                                 covariate == 'baseline_a1c' ~ 'Baseline A1c',
                                 covariate == 'baseline_age' ~ 'Baseline Age',
                                 covariate == 'baseline_bmi' ~ 'Baseline BMI',
                                 covariate == 'bmi_change_1yr' ~ '1-Year BMI Change',
                                 covariate == 'dyslipidemia' ~ 'Dyslipidemia',
                                 covariate == 'gender_M' ~ 'Sex',
                                 covariate == 'hypertension' ~ 'Hypertension',
                                 covariate == 'past_1yr_insulin' ~ 'Insulin Usage Past 12 Months',
                                 covariate == 'race_black' ~ 'Race (African-American Indicator)',
                                 covariate == 'smoking_status_current' ~ 'Smoking Status (Current)',
                                 covariate == 'smoking_status_former' ~ 'Smoking Status (Former)',
                                 covariate == 'smoking_status_never' ~ 'Smoking Status (Never)',
                                 covariate == 'smoking_status_no_self_report' ~ 'Smoking Status (No Self Report)')) %>% 
  mutate('ipw' = fct_relevel(ipw, 'Unadjusted'))



ggplot(df_plot, aes(x = abs(smd), y = covariate)) + 
  facet_wrap(~elig_criteria) + 
  geom_point(aes(color = ipw), size = 3) + 
  labs(x = 'Absolute Standardized Mean Difference',
       y = '',
       color = '') + 
  theme(axis.text = element_text(size = 16),
        strip.text = element_text(size = 16),
        axis.title = element_text(size = 20))

ggsave('figures/results/covariate_balance.png', height = 9, width = 16) 
