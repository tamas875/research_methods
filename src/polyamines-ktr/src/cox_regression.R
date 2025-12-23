library(survival)


#### Create eGFR groups ####
elspa_fu$eGFRgrp = cut(elspa_fu$CKD_EPI_1, breaks = c(-Inf, 30, 60, Inf),
                       labels = c('low', 'medium', 'high'))


#### Define variables ####
pa_vars = c('log2_AcPUT_24h', 'log2_N1AcSPD_24h', 'log2_N8AcSPD_24h',
            'log2_CAD_24h', 'log2_PUT_24h', 'log2_SPM_24h', 'log2_SPD_24h')
pa_names = c('N-acetylputrescine', 'N1-acetylspermidine', 'N8-acetylspermidine',
             'Cadaverine', 'Putrescine', 'Spermine', 'Spermidine')

adj_demo = c('Age_visit_1', 'Geslacht', 'BMI1')
adj_life = c('SBP1', 'Roken_1')
adj_clin = c('log2_U_Tot_Eiwit_24h_totaal_1', 'eGFRgrp')
adj_med = c('Calcineurineremmer_1', 'Dummy_Antihypertensives', 'Statine_1')
adj_diet = c('log2_SumOfvezel', 'log2_SumOfeiwittot', 'log2_SumOfalcohol')

adj_sets = list(
  NULL,
  adj_demo,
  c(adj_demo, adj_life),
  c(adj_demo, adj_life, adj_clin),
  c(adj_demo, adj_life, adj_clin, adj_med),
  c(adj_demo, adj_life, adj_clin, adj_med, adj_diet)
)

n_vars = length(pa_vars)
n_models = length(adj_sets)


#### Graft failure analysis ####
gf_hr = matrix(NA, nrow = n_models, ncol = n_vars)
gf_p = matrix(NA, nrow = n_models, ncol = n_vars)
rownames(gf_hr) = paste0('Model ', 1:n_models)
rownames(gf_p) = paste0('Model ', 1:n_models)
colnames(gf_hr) = pa_names
colnames(gf_p) = pa_names

for (i in seq_along(pa_vars)) {
  var = pa_vars[i]
  
  for (j in seq_along(adj_sets)) {
    # Build formula
    if (is.null(adj_sets[[j]])) {
      formula = as.formula(paste0('Surv(time_gf_fu4, event_gf_fu4) ~ ', var))
    } else {
      formula = as.formula(paste0('Surv(time_gf_fu4, event_gf_fu4) ~ ', var, ' + ',
                                  paste(adj_sets[[j]], collapse = ' + ')))
    }
    
    # Fit model
    model = coxph(formula, data = elspa_fu)
    hr = exp(coef(model))[1]
    ci = exp(confint(model))[1, ]
    p_val = summary(model)$coefficients[1, 'Pr(>|z|)']
    
    gf_hr[j, i] = sprintf('%.2f [%.2f-%.2f]', hr, ci[1], ci[2])
    gf_p[j, i] = p_val
  }
}

# Add significance stars
gf_hr_star = gf_hr
for (j in 1:nrow(gf_p)) {
  p_adj = p.adjust(gf_p[j, ], method = 'BH')
  for (i in 1:ncol(gf_p)) {
    if (!is.na(p_adj[i])) {
      if (p_adj[i] < 0.001) {
        gf_hr_star[j, i] = paste0(gf_hr[j, i], '***')
      } else if (p_adj[i] < 0.01) {
        gf_hr_star[j, i] = paste0(gf_hr[j, i], '**')
      } else if (p_adj[i] < 0.05) {
        gf_hr_star[j, i] = paste0(gf_hr[j, i], '*')
      }
    }
  }
}

print(gf_hr_star)
table(elspa_fu$event_gf_fu4)
summary(elspa_fu$time_gf_fu4)
write.csv(gf_hr_star, 'cox_graft_failure.csv', row.names = TRUE)


#### Mortality analysis ####
mt_hr = matrix(NA, nrow = n_models, ncol = n_vars)
mt_p = matrix(NA, nrow = n_models, ncol = n_vars)
rownames(mt_hr) = paste0('Model ', 1:n_models)
rownames(mt_p) = paste0('Model ', 1:n_models)
colnames(mt_hr) = pa_names
colnames(mt_p) = pa_names

for (i in seq_along(pa_vars)) {
  var = pa_vars[i]
  
  for (j in seq_along(adj_sets)) {
    # Build formula
    if (is.null(adj_sets[[j]])) {
      formula = as.formula(paste0('Surv(time_mort_fu4, event_mort_fu4) ~ ', var))
    } else {
      formula = as.formula(paste0('Surv(time_mort_fu4, event_mort_fu4) ~ ', var, ' + ',
                                  paste(adj_sets[[j]], collapse = ' + ')))
    }
    
    # Fit model
    model = coxph(formula, data = elspa_fu)
    hr = exp(coef(model))[1]
    ci = exp(confint(model))[1, ]
    p_val = summary(model)$coefficients[1, 'Pr(>|z|)']
    
    mt_hr[j, i] = sprintf('%.2f [%.2f-%.2f]', hr, ci[1], ci[2])
    mt_p[j, i] = p_val
  }
}

# Add significance stars
mt_hr_star = mt_hr
for (j in 1:nrow(mt_p)) {
  p_adj = p.adjust(mt_p[j, ], method = 'BH')
  for (i in 1:ncol(mt_p)) {
    if (!is.na(p_adj[i])) {
      if (p_adj[i] < 0.001) {
        mt_hr_star[j, i] = paste0(mt_hr[j, i], '***')
      } else if (p_adj[i] < 0.01) {
        mt_hr_star[j, i] = paste0(mt_hr[j, i], '**')
      } else if (p_adj[i] < 0.05) {
        mt_hr_star[j, i] = paste0(mt_hr[j, i], '*')
      }
    }
  }
}

print(mt_hr_star)
table(elspa_fu$event_mort_fu4)
summary(elspa_fu$time_mort_fu4)
write.csv(mt_hr_star, 'cox_mortality.csv', row.names = TRUE)


#### Extract Model 6 results for plotting ####
plot_gf = data.frame(Polyamines = pa_names, HR = numeric(n_vars), 
                     CI_lower = numeric(n_vars), CI_upper = numeric(n_vars), 
                     P_value = numeric(n_vars))
plot_mt = data.frame(Polyamines = pa_names, HR = numeric(n_vars), 
                     CI_lower = numeric(n_vars), CI_upper = numeric(n_vars), 
                     P_value = numeric(n_vars))

for (i in seq_along(pa_vars)) {
  var = pa_vars[i]
  
  # Graft failure Model 6
  formula_gf = as.formula(paste0('Surv(time_gf_fu4, event_gf_fu4) ~ ', var, ' + ',
                                 paste(adj_sets[[6]], collapse = ' + ')))
  model_gf = coxph(formula_gf, data = elspa_fu)
  hr_gf = exp(coef(model_gf))[1]
  ci_gf = exp(confint(model_gf))[1, ]
  p_gf = summary(model_gf)$coefficients[1, 'Pr(>|z|)']
  
  plot_gf$HR[i] = hr_gf
  plot_gf$CI_lower[i] = ci_gf[1]
  plot_gf$CI_upper[i] = ci_gf[2]
  plot_gf$P_value[i] = p_gf
  
  # Mortality Model 6
  formula_mt = as.formula(paste0('Surv(time_mort_fu4, event_mort_fu4) ~ ', var, ' + ',
                                 paste(adj_sets[[6]], collapse = ' + ')))
  model_mt = coxph(formula_mt, data = elspa_fu)
  hr_mt = exp(coef(model_mt))[1]
  ci_mt = exp(confint(model_mt))[1, ]
  p_mt = summary(model_mt)$coefficients[1, 'Pr(>|z|)']
  
  plot_mt$HR[i] = hr_mt
  plot_mt$CI_lower[i] = ci_mt[1]
  plot_mt$CI_upper[i] = ci_mt[2]
  plot_mt$P_value[i] = p_mt
}

plot_gf$P_adj = p.adjust(plot_gf$P_value, method = 'BH')
plot_gf$HR_CI = sprintf('%.2f [%.2f-%.2f]', plot_gf$HR, plot_gf$CI_lower, plot_gf$CI_upper)
plot_gf$P_fmt = ifelse(plot_gf$P_adj < 0.001, '< 0.001', sprintf('%.3f', plot_gf$P_adj))

plot_mt$P_adj = p.adjust(plot_mt$P_value, method = 'BH')
plot_mt$HR_CI = sprintf('%.2f [%.2f-%.2f]', plot_mt$HR, plot_mt$CI_lower, plot_mt$CI_upper)
plot_mt$P_fmt = ifelse(plot_mt$P_adj < 0.001, '< 0.001', sprintf('%.3f', plot_mt$P_adj))

print(plot_gf)
print(plot_mt)
