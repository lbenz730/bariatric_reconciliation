library(tidyverse)
library(arrow)
library(glue)
library(gt)
source('scripts/helpers.R')

### Load Data
df_ci <- 
  read_parquet('data/results/cvd_cumulative_incidence_by_formula.parquet') %>% 
  rename('ipw_formula' = formulae) %>% 
  mutate('ipw_formula' = 
           fct_relevel(ipw_formula, c('Linear Main Effects', 'Madenci Proxy'))) %>% 
  mutate('elig_criteria' = case_when(elig_criteria == 'No Pre-Op Restrictions' ~ 'Trial #1A',
                                     elig_criteria == 'Pre-Op Restrictions' ~ 'Trial #1B',
                                     T ~ 'Trial #1C')) %>% 
  mutate('ipw_formula' = case_when(ipw_formula == 'Linear Main Effects' ~ 'Propensity Score Model #1',
                                   ipw_formula == 'Madenci Proxy' ~ 'Propensity Score Model #2'))

df_bootstraps <- 
  map_dfr(dir('data/bootstraps/model_comp', full.names = T), read_parquet) %>% 
  mutate('ipw_formula' = fct_relevel(ipw_formula,
                                     'Linear Main Effects', 
                                     'Madenci Proxy')) %>% 
  mutate('elig_criteria' = case_when(elig_criteria == 'No Pre-Op Restrictions' ~ 'Trial #1A',
                                     elig_criteria == 'Pre-Op Restrictions' ~ 'Trial #1B',
                                     T ~ 'Trial #1C')) %>% 
  mutate('ipw_formula' = case_when(ipw_formula == 'Linear Main Effects' ~ 'Propensity Score Model #1',
                                   ipw_formula == 'Madenci Proxy' ~ 'Propensity Score Model #2'))

### Point Estimate Plot
ggplot(df_ci, aes(x = time/12, y = cum_incidence)) + 
  facet_wrap(~elig_criteria) + 
  geom_line(aes(color = as.factor(surgery), lty = ipw_formula), lwd = 1) + 
  scale_x_continuous(breaks = seq(0, 7)) + 
  scale_y_continuous(labels = scales::percent) + 
  scale_color_discrete(labels = c('No Surgery', 'Surgery')) + 
  labs(x = 'Time Since Baseline (Years)',
       y = 'Cumulative Incidence of CVD',
       color = ' Treatment Group',
       linetype = 'Probability of Treatment Model for IPW',
       title = 'Adjusted Cumulative Incidence of CVD',
       subtitle = 'Comparison by Study Eligibility Criteria') + 
  guides(linetype = guide_legend(nrow = 3)) + 
  theme(legend.box = 'vertical',
        legend.title = element_text(size = 14)) 

ggsave('figures/results/adjusted_CI_PLR_point_est.png', height = 9/1.2, width = 16/1.2)


### Compute confidence intervals
df_summary <- 
  df_bootstraps %>% 
  group_by(ipw_formula, elig_criteria, outcome_model, outcome_formula, surgery, time) %>% 
  summarise('mean_ci' = mean(cum_incidence),
            'lower_ci' = 2 * mean(cum_incidence) - quantile(cum_incidence, 0.975),
            'upper_ci' = 2 * mean(cum_incidence) - quantile(cum_incidence, 0.025)) %>% 
  ungroup() 


ggplot(df_summary, aes(x = time/12, y = mean_ci)) + 
  facet_grid(ipw_formula~elig_criteria, labeller = label_wrap_gen(width = 30)) + 
  geom_ribbon(aes(ymin = lower_ci, ymax = upper_ci, fill = as.factor(surgery)), alpha = 0.2) +
  geom_line(data = df_ci, aes(y = cum_incidence, color = as.factor(surgery)), lwd = 1.5) +
  labs(x = 'Time Since Baseline (Years)',
       y = 'Cumulative Incidence of Cardiovascular Events',
       color = '',
       fill = '') + 
  scale_x_continuous(breaks = seq(0, 7)) + 
  scale_y_continuous(labels = scales::percent) + 
  scale_color_discrete(labels = c('No Surgery', 'Surgery')) + 
  scale_fill_discrete(labels = c('No Surgery', 'Surgery')) + 
  theme(strip.text = element_text(size = 14),
        axis.text = element_text(size = 14),
        panel.spacing=unit(0.7,"lines"))

ggsave('figures/results/adjusted_CI_PLR.png', height = 9, width = 16)


### Risk Table 
df_table <-
  df_bootstraps %>% 
  group_by(boot_id) %>% 
  filter(!any(is.na(ipw_formula))) %>% 
  group_by(ipw_formula, elig_criteria, time, outcome_model, outcome_formula, boot_id) %>% 
  summarise('risk_ratio' = cum_incidence[surgery == 1]/cum_incidence[surgery == 0],
            'risk_difference' = cum_incidence[surgery == 1] - cum_incidence[surgery == 0]) %>% 
  group_by(ipw_formula, elig_criteria, time, outcome_model, outcome_formula) %>% 
  summarise('lower_rr' = 2 * mean(risk_ratio, na.rm = T) - quantile(risk_ratio, 0.975, na.rm = T),
            'upper_rr' = 2 * mean(risk_ratio, na.rm = T) - quantile(risk_ratio, 0.025, na.rm = T),
            'lower_rd' = 2 * mean(risk_difference, na.rm = T) - quantile(risk_difference, 0.975, na.rm = T),
            'upper_rd' = 2 * mean(risk_difference, na.rm = T) - quantile(risk_difference, 0.025, na.rm = T)) %>% 
  ungroup() %>% 
  inner_join(
    df_ci %>% 
      group_by(ipw_formula, elig_criteria, time, outcome_model) %>% 
      summarise('mean_rr' = cum_incidence[surgery == 1]/cum_incidence[surgery == 0],
                'mean_rd' = cum_incidence[surgery == 1] - cum_incidence[surgery == 0])
  ) %>% 
  filter(time %in% c(12, 36, 60, 84)) %>% 
  mutate('rr_ci' = paste0(sprintf('%0.2f', mean_rr), ' (', sprintf('%0.2f', lower_rr), ', ', sprintf('%0.2f', upper_rr), ')'),
         'rd_ci' = paste0(sprintf('%0.1f', mean_rd * 100), '% (', sprintf('%0.1f', lower_rd * 100), '% , ', sprintf('%0.1f', upper_rd * 100), '%)')) %>% 
  select(elig_criteria, time, ipw_formula, rr_ci, rd_ci) %>% 
  pivot_wider(names_from = 'time',
              values_from = c('rr_ci', 'rd_ci')) %>% 
  select(ipw_formula, elig_criteria, ends_with('_12'), ends_with('_36'), ends_with('_60'), ends_with('_84')) %>% 
  group_by(elig_criteria)

gt_risk <- 
  gt(df_table) %>% 
  cols_align('center') %>% 
  tab_spanner(columns = ends_with('_12'), '1-Year') %>%
  tab_spanner(columns = ends_with('_36'), '3-Year') %>%
  tab_spanner(columns = ends_with('_60'), '5-Year') %>%
  tab_spanner(columns = ends_with('_84'), '7-Year') %>% 
  tab_style(
    style = list(
      cell_borders(
        sides = "right",
        color = "black",
        weight = px(3)
      )
    ),
    locations = list(
      cells_body(
        columns = c(ipw_formula, contains('rd')))
    )
  ) %>% 
  tab_header(title = md('**Risk Ratios (RR) and Risk Difference (RD) For Surgery vs. No-Surgery**'),
             subtitle = md('**Breakdown by IPW Formula and Eligibility Criteria**')) %>% 
  tab_source_note('Cells shaded red denote non-statistically significant relationships at the 0.05 significance level') %>%  
  tab_options(column_labels.font.size = 20,
              heading.title.font.size = 36,
              heading.subtitle.font.size = 30,
              heading.title.font.weight = 'bold',
              heading.subtitle.font.weight = 'bold',
              column_labels.font.weight = 'bold',
              row_group.font.weight = 'bold') %>% 
  cols_label('elig_criteria' = 'Eligibility Criteria',
             'ipw_formula' = '',
             'rr_ci_12' = 'RR',
             'rd_ci_12' = 'RD',
             'rr_ci_36' = 'RR',
             'rd_ci_36' = 'RD',
             'rr_ci_60' = 'RR',
             'rd_ci_60' = 'RD',
             'rr_ci_84' = 'RR',
             'rd_ci_84' = 'RD')


for(i in 1:nrow(df_table)) {
  for(j in which(grepl('rd', names(df_table)) | grepl('rr', names(df_table)))) {
    
    if(j %% 2 == 1) {
      vec <- as.numeric(unlist(str_extract_all(df_table[[j]][i], '\\d\\.\\d+')))
      if(vec[3] >= 1) {
        gt_risk <- 
          gt_risk %>% 
          tab_style(style = cell_fill(color = 'pink'),
                    locations = cells_body(columns = j, 
                                           rows = i))
      }
    } else {
      vec <- as.numeric(unlist(str_extract_all(gsub('-0\\.0', '-0.1', df_table[[j]][i]), '-?\\d\\.\\d+')))
      if(vec[3] >= 0) {
        gt_risk <- 
          gt_risk %>% 
          tab_style(style = cell_fill(color = 'pink'),
                    locations = cells_body(columns = j, 
                                           rows = i))
      }
    }
  }
}

gt_risk
cluster_gtsave(gt_risk, 'figures/tables/risk_table.html', keep_html = T, zoom = 3, vwidth = 1800)
