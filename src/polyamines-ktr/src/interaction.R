library(jstable)
library(tidyverse)


#### Prepare data ####
elspa_fu3 = elspa_fu %>%
  select(time_gf_fu4, event_gf_fu4, time_mort_fu4, event_mort_fu4,
         log2_AcPUT_24h, log2_N8AcSPD_24h, log2_N1AcSPD_24h,
         log2_PUT_24h, log2_SPD_24h, log2_SPM_24h, log2_CAD_24h,
         Age_visit_1, Geslacht, BMI1, SBP1, Roken_1,
         Dummy_Proteinurie, eGFRgrp, Calcineurineremmer_1,
         Dummy_Antihypertensives, Statine_1, log2_SumOfvezel,
         log2_SumOfeiwittot, log2_SumOfalcohol)


#### Create stratified variables ####
# Age
median_age = median(elspa_fu3$Age_visit_1, na.rm = TRUE)
elspa_fu3$Age_visit_1_strat = cut(elspa_fu3$Age_visit_1,
                                   c(-Inf, median_age, Inf),
                                   labels = c(paste0('< ', round(median_age, 1)),
                                             paste0('> ', round(median_age, 1))))

# BMI
median_bmi = median(elspa_fu3$BMI1, na.rm = TRUE)
elspa_fu3$BMI1_strat = cut(elspa_fu3$BMI1,
                           c(-Inf, median_bmi, Inf),
                           labels = c(paste0('< ', round(median_bmi, 1)),
                                     paste0('> ', round(median_bmi, 1))))

# SBP
median_sbp = median(elspa_fu3$SBP1, na.rm = TRUE)
elspa_fu3$SBP1_strat = cut(elspa_fu3$SBP1,
                           c(-Inf, median_sbp, Inf),
                           labels = c(paste0('< ', round(median_sbp, 1)),
                                     paste0('> ', round(median_sbp, 1))))

# Fiber intake
median_fiber = median(elspa_fu3$log2_SumOfvezel, na.rm = TRUE)
elspa_fu3$log2_SumOfvezel_strat = cut(elspa_fu3$log2_SumOfvezel,
                                       c(-Inf, median_fiber, Inf),
                                       labels = c(paste0('< ', round(median_fiber, 1)),
                                                 paste0('> ', round(median_fiber, 1))))

# Protein intake
median_protein = median(elspa_fu3$log2_SumOfeiwittot, na.rm = TRUE)
elspa_fu3$log2_SumOfeiwittot_strat = cut(elspa_fu3$log2_SumOfeiwittot,
                                          c(-Inf, median_protein, Inf),
                                          labels = c(paste0('< ', round(median_protein, 1)),
                                                    paste0('> ', round(median_protein, 1))))

# Alcohol intake
median_alcohol = median(elspa_fu3$log2_SumOfalcohol, na.rm = TRUE)
elspa_fu3$log2_SumOfalcohol_strat = cut(elspa_fu3$log2_SumOfalcohol,
                                         c(-Inf, median_alcohol, Inf),
                                         labels = c(paste0('< ', round(median_alcohol, 1)),
                                                   paste0('> ', round(median_alcohol, 1))))


#### Define subgroups and covariates ####
subgroups = c('Age_visit_1_strat', 'Geslacht', 'BMI1_strat', 'SBP1_strat',
              'Roken_1', 'Dummy_Proteinurie', 'eGFRgrp',
              'Calcineurineremmer_1', 'Dummy_Antihypertensives', 'Statine_1',
              'log2_SumOfvezel_strat', 'log2_SumOfeiwittot_strat',
              'log2_SumOfalcohol_strat')

covs = c('Age_visit_1_strat', 'Geslacht', 'BMI1_strat', 'SBP1_strat',
         'Roken_1', 'Dummy_Proteinurie', 'eGFRgrp',
         'Calcineurineremmer_1', 'Dummy_Antihypertensives', 'Statine_1',
         'log2_SumOfvezel_strat', 'log2_SumOfeiwittot_strat',
         'log2_SumOfalcohol_strat')


#### Graft failure interaction analysis ####
# N1-acetylspermidine
els_int_N1AcSPD = TableSubgroupMultiCox(
  formula = Surv(time_gf_fu4, event_gf_fu4) ~ log2_N1AcSPD_24h,
  var_subgroups = subgroups,
  var_cov = covs,
  data = elspa_fu3
)

# N8-acetylspermidine
els_int_N8AcSPD = TableSubgroupMultiCox(
  formula = Surv(time_gf_fu4, event_gf_fu4) ~ log2_N8AcSPD_24h,
  var_subgroups = subgroups,
  var_cov = covs,
  data = elspa_fu3
)


#### Mortality interaction analysis ####
# N-acetylputrescine
els_int_AcPUT = TableSubgroupMultiCox(
  formula = Surv(time_mort_fu4, event_mort_fu4) ~ log2_AcPUT_24h,
  var_subgroups = subgroups,
  var_cov = covs,
  data = elspa_fu3
)

# N1-acetylspermidine
els_int_N1AcSPD_mort = TableSubgroupMultiCox(
  formula = Surv(time_mort_fu4, event_mort_fu4) ~ log2_N1AcSPD_24h,
  var_subgroups = subgroups,
  var_cov = covs,
  data = elspa_fu3
)

# N8-acetylspermidine
els_int_N8AcSPD_mort = TableSubgroupMultiCox(
  formula = Surv(time_mort_fu4, event_mort_fu4) ~ log2_N8AcSPD_24h,
  var_subgroups = subgroups,
  var_cov = covs,
  data = elspa_fu3
)

# Spermidine
els_int_SPD_mort = TableSubgroupMultiCox(
  formula = Surv(time_mort_fu4, event_mort_fu4) ~ log2_SPD_24h,
  var_subgroups = subgroups,
  var_cov = covs,
  data = elspa_fu3
)


#### Table 3S: Graft failure interaction ####
# Extract N1-acetylspermidine data (rows 40-42)
data_subset = els_int_N1AcSPD[40:42, ]

table_3s = data.frame(
  Polyamines = rep('N1-acetylspermidine', nrow(data_subset)),
  Variable = data_subset$Variable,
  Count = data_subset$Count,
  Percent = data_subset$Percent,
  HR_CI = sprintf('%.2f (%.2f, %.2f)', 
                  as.numeric(data_subset$'Point Estimate'),
                  as.numeric(data_subset$Lower),
                  as.numeric(data_subset$Upper)),
  P_for_interaction = data_subset$'P for interaction'
)

print(table_3s)
write.csv(table_3s, 'Table_3S_Graft_Failure_Interaction.csv', row.names = FALSE)


#### Table 5S: Mortality interaction ####
# N-acetylputrescine (rows 14-16)
data_acput = els_int_AcPUT[14:16, ]

table_5s_acput = data.frame(
  Polyamines = rep('N-acetylputrescine', nrow(data_acput)),
  Variable = data_acput$Variable,
  Count = data_acput$Count,
  Percent = data_acput$Percent,
  HR_CI = sprintf('%.2f (%.2f, %.2f)', 
                  as.numeric(data_acput$'Point Estimate'),
                  as.numeric(data_acput$Lower),
                  as.numeric(data_acput$Upper)),
  P_for_interaction = data_acput$'P for interaction'
)

# Spermidine (rows 5-7, 17-19)
data_spd = els_int_SPD_mort[c(5:7, 17:19), ]

table_5s_spd = data.frame(
  Polyamines = rep('Spermidine', nrow(data_spd)),
  Variable = data_spd$Variable,
  Count = data_spd$Count,
  Percent = data_spd$Percent,
  HR_CI = sprintf('%.2f (%.2f, %.2f)', 
                  as.numeric(data_spd$'Point Estimate'),
                  as.numeric(data_spd$Lower),
                  as.numeric(data_spd$Upper)),
  P_for_interaction = data_spd$'P for interaction'
)

# Combine
table_5s = rbind(table_5s_acput, table_5s_spd)

print(table_5s)
write.csv(table_5s, 'Table_5S_Mortality_Interaction.csv', row.names = FALSE)
