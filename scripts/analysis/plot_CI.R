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
           fct_relevel(ipw_formula, c('Linear Main Effects', 
                                      'Madenci Proxy (NCS on Continuous Covariates, No Smoking or 1-Year A1c Change)',
                                      'Best BIC Model (Splines on BMI/Age, No A1c 1yr Change)'))) %>% 
  mutate('elig_criteria' = case_when(elig_criteria == 'No Pre-Op Restrictions' ~ 'Eligibility Set #1',
                                     elig_criteria == 'Pre-Op Restrictions' ~ 'Eligibility Sets #1 + #2',
                                     T ~ 'Eligibility Sets #1 + #2 + #3')) %>% 
  mutate('ipw_formula' = case_when(ipw_formula == 'Linear Main Effects' ~ 'Propensity Score Model #1',
                                   ipw_formula == 'Madenci Proxy (NCS on Continuous Covariates, No Smoking or 1-Year A1c Change)' ~ 'Propensity Score Model #2',
                                   T ~ 'Propensity Score Model #3'))
                                   
                                   

df_bootstraps <- 
  map_dfr(dir('data/bootstraps/model_comp', full.names = T), read_parquet) %>% 
  mutate('ipw_formula' = fct_relevel(ipw_formula,
                                     'Linear Main Effects', 
                                     'Madenci Proxy (NCS on Continuous Covariates, No Smoking or 1-Year A1c Change)',
                                     'Best BIC Model (Splines on BMI/Age, No A1c 1yr Change)')) %>% 
  mutate('elig_criteria' = case_when(elig_criteria == 'No Pre-Op Restrictions' ~ 'Eligibility Set #1',
                                     elig_criteria == 'Pre-Op Restrictions' ~ 'Eligibility Sets #1 + #2',
                                     T ~ 'Eligibility Sets #1 + #2 + #3')) %>% 
  mutate('ipw_formula' = case_when(ipw_formula == 'Linear Main Effects' ~ 'Propensity Score Model #1',
                                   ipw_formula == 'Madenci Proxy (NCS on Continuous Covariates, No Smoking or 1-Year A1c Change)' ~ 'Propensity Score Model #2',
                                   T ~ 'Propensity Score Model #3'))

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
            'lower_ci' = mean(cum_incidence) + qnorm(0.025) * sd(cum_incidence),
            'upper_ci' = mean(cum_incidence) + qnorm(0.975) * sd(cum_incidence)) %>% 
  ungroup() 




ggplot(df_summary, aes(x = time/12, y = mean_ci)) + 
  facet_grid(elig_criteria~ipw_formula, labeller = label_wrap_gen(width = 30)) + 
  geom_ribbon(aes(ymin = lower_ci, ymax = upper_ci, fill = as.factor(surgery)), alpha = 0.2) +
  geom_line(data = df_ci, aes(y = cum_incidence, color = as.factor(surgery)), lwd = 1.5) +
  labs(x = 'Time Since Baseline (Years)',
       y = 'Cumulative Incidence of Cardiovascular Events',
       color = '',
       fill = ''
  ) + 
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
  group_by(ipw_formula, elig_criteria, time, outcome_model, outcome_formula, boot_id) %>% 
  summarise('risk_ratio' = cum_incidence[surgery == 1]/cum_incidence[surgery == 0],
            'risk_difference' = cum_incidence[surgery == 1] - cum_incidence[surgery == 0]) %>% 
  group_by(ipw_formula, elig_criteria, time, outcome_model, outcome_formula) %>% 
  summarise('lower_rr' = mean(risk_ratio) + qnorm(0.025) * sd(risk_ratio),
            'upper_rr' = mean(risk_ratio) + qnorm(0.975) * sd(risk_ratio),
            'lower_rd' = mean(risk_difference) + qnorm(0.025) * sd(risk_difference),
            'upper_rd' = mean(risk_difference) + qnorm(0.975) * sd(risk_difference)) %>% 
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
  group_by(ipw_formula)

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
        columns = c(elig_criteria, contains('rd')))
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
# 
# 
# ### For DURABLE
# durable_ci <- read_parquet('data/results/cvd_cumulative_incidence_DURABLE.parquet')
# durable_bootstraps <- map_dfr(dir('data/bootstraps/DURABLE/', full.names = T), read_parquet)
# 
# ggplot(durable_ci, aes(x = time/12, y = cum_incidence)) + 
#   facet_grid(eligibility ~ weights, labeller = label_wrap_gen(width = 40, multi_line = T)) + 
#   geom_line(aes(color = as.factor(surgery))) + 
#   scale_x_continuous(breaks = seq(0, 7)) + 
#   scale_y_continuous(labels = scales::percent) + 
#   scale_color_discrete(labels = c('No Surgery', 'Surgery')) + 
#   labs(x = 'Time Since Baseline (Years)',
#        y = 'Cumulative Incidence of CVD',
#        color = ' Treatment Group',
#        title = 'Adjusted Cumulative Incidence of CVD',
#        subtitle = 'DURABLE Target Trial(s)') + 
#   guides(linetype = guide_legend(nrow = 3)) + 
#   theme(legend.box = 'vertical',
#         legend.title = element_text(size = 14))
# 
# ggsave('figures/results/DURABLE_adjusted_CI_PLR_point_est.png', height = 9, width = 16)
# 
# 
# 
# 
# 
# 
# durable_summary <- 
#   durable_bootstraps %>% 
#   group_by(eligibility, weights, surgery, time) %>% 
#   summarise('mean_ci' = mean(cum_incidence),
#             'lower_ci' = mean(cum_incidence) + qnorm(0.025) * sd(cum_incidence),
#             'upper_ci' = mean(cum_incidence) + qnorm(0.975) * sd(cum_incidence)) %>% 
#   ungroup()
# 
# ggplot(durable_summary, aes(x = time/12, y = mean_ci)) + 
#   facet_grid(weights ~ eligibility, labeller = label_wrap_gen(width = 25)) + 
#   geom_ribbon(aes(ymin = lower_ci, ymax = upper_ci, fill = as.factor(surgery)), alpha = 0.2) +
#   geom_line(data = durable_ci, aes(y = cum_incidence, color = as.factor(surgery))) +
#   labs(x = 'Time Since Baseline (Years)',
#        y = 'Cumulative Incidence of CVD',
#        color = '',
#        fill = '',
#        title = 'Adjusted Cumulative Incidence of CVD',
#        subtitle = 'DURABLE Target Trial(s)') + 
#   scale_x_continuous(breaks = seq(0, 7)) + 
#   scale_y_continuous(labels = scales::percent) + 
#   scale_color_discrete(labels = c('No Surgery', 'Surgery')) + 
#   scale_fill_discrete(labels = c('No Surgery', 'Surgery'))
# 
# ggsave('figures/results/DURABLE_adjusted_CI_PLR.png', height = 9, width = 16)
# 
# 
# ci_table <- 
#   durable_ci %>% 
#   group_by(eligibility, weights, time) %>% 
#   summarise('rr' = cum_incidence[surgery == 1]/cum_incidence[surgery == 0],
#             'rd' = cum_incidence[surgery == 1] - cum_incidence[surgery == 0]) %>% 
#   filter(time %in% c(12, 36, 60, 84)) %>% 
#   select(eligibility, time, weights, rr, rd) %>% 
#   pivot_wider(names_from = 'time',
#               values_from = c('rr', 'rd')) %>% 
#   select(eligibility, weights, ends_with('_12'), ends_with('_36'), ends_with('_60'), ends_with('_84')) %>% 
#   group_by(eligibility) 
# 
# gt_durable_alt <- 
#   gt(ci_table) %>% 
#   fmt_number(contains('rr'), decimals = 3) %>% 
#   fmt_percent(contains('rd'), decimals = 2) %>% 
#   tab_spanner(columns = ends_with('_12'), '1-Year') %>% 
#   tab_spanner(columns = ends_with('_36'), '3-Year') %>% 
#   tab_spanner(columns = ends_with('_60'), '5-Year') %>% 
#   tab_spanner(columns = ends_with('_84'), '7-Year') %>% 
#   tab_style(
#     style = list(
#       cell_borders(
#         sides = "right",
#         color = "black",
#         weight = px(3)
#       )
#     ),
#     locations = list(
#       cells_body(
#         columns = c(weights, contains('rd')))
#     )
#   ) %>% 
#   tab_header(title = md('**Risk Ratios (RR) and Risk Difference (RD) For Surgery vs. No-Surgery**'),
#              subtitle = md('**Breakdown by IPW Type and Eligibility Criteria**')) %>% 
#   tab_options(column_labels.font.size = 20,
#               heading.title.font.size = 36,
#               heading.subtitle.font.size = 30,
#               heading.title.font.weight = 'bold',
#               heading.subtitle.font.weight = 'bold',
#               column_labels.font.weight = 'bold',
#               row_group.font.weight = 'bold') %>% 
#   cols_label('weights' = 'IPW',
#              'rr_12' = 'RR',
#              'rd_12' = 'RD',
#              'rr_36' = 'RR',
#              'rd_36' = 'RD',
#              'rr_60' = 'RR',
#              'rd_60' = 'RD',
#              'rr_84' = 'RR',
#              'rd_84' = 'RD')
# cluster_gtsave(gt_durable_alt, 'figures/tables/DURABLE_risk_table_no_conf_int.png', vwidth = 1200, vheight = 600)
# 
# 
# 
# durable_table <-
#   durable_bootstraps %>% 
#   group_by(eligibility, weights, time, boot_id) %>% 
#   summarise('risk_ratio' = cum_incidence[surgery == 1]/cum_incidence[surgery == 0],
#             'risk_difference' = cum_incidence[surgery == 1] - cum_incidence[surgery == 0]) %>% 
#   group_by(eligibility, weights, time) %>% 
#   summarise('lower_rr' = mean(risk_ratio) + qnorm(0.025) * sd(risk_ratio),
#             'upper_rr' = mean(risk_ratio) + qnorm(0.975) * sd(risk_ratio),
#             'lower_rd' = mean(risk_difference) + qnorm(0.025) * sd(risk_difference),
#             'upper_rd' = mean(risk_difference) + qnorm(0.975) * sd(risk_difference)) %>% 
#   ungroup() %>% 
#   inner_join(
#     durable_ci %>% 
#       group_by(eligibility, weights, time) %>% 
#       summarise('mean_rr' = cum_incidence[surgery == 1]/cum_incidence[surgery == 0],
#                 'mean_rd' = cum_incidence[surgery == 1] - cum_incidence[surgery == 0])
#   ) %>% 
#   filter(time %in% c(12, 36, 60, 84)) %>% 
#   mutate('rr_ci' = paste0(sprintf('%0.3f', mean_rr), ' (', sprintf('%0.3f', lower_rr), ', ', sprintf('%0.3f', upper_rr), ')'),
#          'rd_ci' = paste0(sprintf('%0.3f', mean_rd * 100), '% (', sprintf('%0.3f', lower_rd * 100), '% , ', sprintf('%0.3f', upper_rd * 100), '%)')) %>% 
#   select(eligibility, time, weights, rr_ci, rd_ci) %>% 
#   pivot_wider(names_from = 'time',
#               values_from = c('rr_ci', 'rd_ci')) %>% 
#   select(eligibility, weights, ends_with('_12'), ends_with('_36'), ends_with('_60'), ends_with('_84')) %>% 
#   group_by(eligibility)
# 
# 
# 
# 
# gt_risk_durable <- 
#   gt(durable_table) %>% 
#   cols_align('center') %>% 
#   tab_spanner(columns = ends_with('_12'), '1-Year') %>% 
#   tab_spanner(columns = ends_with('_36'), '3-Year') %>% 
#   tab_spanner(columns = ends_with('_60'), '5-Year') %>% 
#   tab_spanner(columns = ends_with('_84'), '7-Year') %>% 
#   tab_style(
#     style = list(
#       cell_borders(
#         sides = "right",
#         color = "black",
#         weight = px(3)
#       )
#     ),
#     locations = list(
#       cells_body(
#         columns = c(weights, contains('rd')))
#     )
#   ) %>% 
#   tab_header(title = md('**Risk Ratios (RR) and Risk Difference (RD) For Surgery vs. No-Surgery**'),
#              subtitle = md('**Breakdown by IPW Type and Eligibility Criteria**')) %>% 
#   tab_source_note('Cells shaded red denote non-statistically significant relationships at the 0.05 significance level') %>%  
#   tab_options(column_labels.font.size = 20,
#               heading.title.font.size = 36,
#               heading.subtitle.font.size = 30,
#               heading.title.font.weight = 'bold',
#               heading.subtitle.font.weight = 'bold',
#               column_labels.font.weight = 'bold',
#               row_group.font.weight = 'bold') %>% 
#   cols_label('weights' = 'IPW',
#              'rr_ci_12' = 'RR',
#              'rd_ci_12' = 'RD',
#              'rr_ci_36' = 'RR',
#              'rd_ci_36' = 'RD',
#              'rr_ci_60' = 'RR',
#              'rd_ci_60' = 'RD',
#              'rr_ci_84' = 'RR',
#              'rd_ci_84' = 'RD')
# 
# 
# for(i in 1:nrow(durable_table)) {
#   for(j in which(grepl('rd', names(durable_table)) | grepl('rr', names(durable_table)))) {
#     
#     if(j %% 2 == 1) {
#       vec <- as.numeric(unlist(str_extract_all(durable_table[[j]][i], '\\d\\.\\d+')))
#       if(vec[3] >= 1) {
#         gt_risk_durable <- 
#           gt_risk_durable %>% 
#           tab_style(style = cell_fill(color = 'pink'),
#                     locations = cells_body(columns = j, 
#                                            rows = i))
#       }
#     } else {
#       vec <- as.numeric(unlist(str_extract_all(durable_table[[j]][i], '-?\\d\\.\\d+')))
#       if(vec[3] >= 0) {
#         gt_risk_durable <- 
#           gt_risk_durable %>% 
#           tab_style(style = cell_fill(color = 'pink'),
#                     locations = cells_body(columns = j, 
#                                            rows = i))
#       }
#     }
#   }
# }
# 
# gt_risk_durable
# cluster_gtsave(gt_risk_durable, 'figures/tables/DURABLE_risk_table.png', vwidth = 1800, zoom = 3)
# 
# 
# 
# ### For Sebastien's Talks
# df1 <- 
#   durable_ci %>% 
#   filter(weights == 'IPW: Confounding') %>% 
#   filter(eligibility == 'Eligibility: DURABLE (No Baseline A1c Required)')
# 
# df2 <- 
#   durable_ci %>% 
#   filter(weights == 'IPW: Confounding') %>% 
#   filter(eligibility == 'Eligibility: DURABLE (Baseline A1c Required)')
# 
# durable_summary <- 
#   durable_bootstraps %>% 
#   group_by(weights, eligibility, surgery, time) %>% 
#   summarise('mean_ci' = mean(cum_incidence),
#             'lower_ci' = mean(cum_incidence) + qnorm(0.025) * sd(cum_incidence),
#             'upper_ci' = mean(cum_incidence) + qnorm(0.975) * sd(cum_incidence)) %>% 
#   ungroup()
# 
# 
# durable_summary %>% 
#   filter(weights == 'IPW: Confounding') %>% 
#   filter(eligibility == 'Eligibility: DURABLE (No Baseline A1c Required)') %>% 
#   ggplot(aes(x = time/12, y = mean_ci)) + 
#   geom_ribbon(aes(ymin = lower_ci, ymax = upper_ci, fill = as.factor(surgery)), alpha = 0.2) +
#   geom_line(data = df1, aes(y = cum_incidence, color = as.factor(surgery))) +
#   labs(x = 'Time Since Baseline (Years)',
#        y = 'Cumulative Incidence of CVD',
#        color = '',
#        fill = '',
#        title = 'Adjusted Cumulative Incidence of CVD',
#        subtitle = 'Eligibility: DURABLE (No Baseline A1c Required)') + 
#   scale_x_continuous(breaks = seq(0, 7)) + 
#   scale_y_continuous(labels = scales::percent) + 
#   scale_color_discrete(labels = c('No Surgery', 'Surgery')) + 
#   scale_fill_discrete(labels = c('No Surgery', 'Surgery'))
# ggsave('figures/results/adjusted_CI_DURABLE_confounding_only.png', height = 9/1.2, width = 16/1.2)
# 
# durable_summary %>% 
#   filter(weights == 'IPW: Confounding') %>% 
#   filter(eligibility == 'Eligibility: DURABLE') %>% 
#   ggplot(aes(x = time/12, y = mean_ci)) + 
#   geom_ribbon(aes(ymin = lower_ci, ymax = upper_ci, fill = as.factor(surgery)), alpha = 0.2) +
#   geom_line(data = df2, aes(y = cum_incidence, color = as.factor(surgery))) +
#   labs(x = 'Time Since Baseline (Years)',
#        y = 'Cumulative Incidence of CVD',
#        color = '',
#        fill = '',
#        title = 'Adjusted Cumulative Incidence of CVD',
#        subtitle = 'Eligibility: DURABLE (Baseline A1c Required)') + 
#   scale_x_continuous(breaks = seq(0, 7)) + 
#   scale_y_continuous(labels = scales::percent) + 
#   scale_color_discrete(labels = c('No Surgery', 'Surgery')) + 
#   scale_fill_discrete(labels = c('No Surgery', 'Surgery'))
# ggsave('figures/results/adjusted_CI_DURABLE_confounding_only_A1c_required.png', height = 9/1.2, width = 16/1.2)
