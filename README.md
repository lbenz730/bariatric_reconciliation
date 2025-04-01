# Bariatric Surgery Target Trial Emulation

Code for the paper:

Haneuse, S., Benz, L., Smith, V., Arterburn, D., and Macejewski, M.L. "Long-term cardiovascular outcomes following bariatric surgery: Reconciling seemingly conflicting evidence." 

## R Scripts (`scripts/`)

__helpers.R__: Useful helper functions


### Data Processing (`scripts/data`) 
* __cvd_outcomes.R__: Build dataset of CVD diagnoses
* __obesity_comorbidities.R__: Build dataset of obesity comorbidities (e.g., hypertension and dyslipidemia)
* __exclusion_criteria.R__: Additional exclusion criteria processing (Eligibility Set \# 3)
* __weight_cleaning_part2.R__: Some additional pre-processing of body weight measurements
* __build_bariatric_CVD_trial_datasets.R__: Script to build a sequence of target trials
* __combine_trials.R__: Combine trials into one big file for analysis

### Analysis (`scripts/analysis`) 
* __fit_models.R__: Script to fit pooled logistic regression models for confounding and outcome to compute adjusted cumulative incidence curves
* __bootstrap_CI.R__: Bootstrap confidence intervals for cumulative incidence curves (Figure 1, Table 2)
* __plot_CI.R__: Plot cumulative incidence curves
* __population_summary.R__: Script to make summary tables of the population (Table 1)
* __spline_BIC.R__: Analysis of BIC for various IPW models using NCS on various continuous covariates (Propensity Model \#3)
* __patient_flowchart.R__: Numbers for patient flow chart (eFigure 1)
* __covariate_balance.R__: Use `cobalt` to make covariate balance plot (eFigure 2)

## Jobs (`jobs/`)
* __build_cvd_trials.sh__: Build target trials in sequence
* __bootstrap_CI.sh__: Bootstrap confidence intervals for cumulative incidence curves
