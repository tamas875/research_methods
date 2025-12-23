library(survival)


#### Median imputation ####
for (col in names(elspa)) {
  if (is.numeric(elspa[[col]])) {
    median_val = median(elspa[[col]], na.rm = TRUE)
    elspa[[col]][is.na(elspa[[col]])] = median_val
  } else if (is.character(elspa[[col]]) || is.factor(elspa[[col]])) {
    mode_val = names(sort(table(elspa[[col]]), decreasing = TRUE))[1]
    elspa[[col]][is.na(elspa[[col]])] = mode_val
  }
}

elspa_median = elspa


#### Convert to log2 scale ####
names(elspa_median)
log2_var = names(elspa_median)[c(7, 26:40, 59:65)]
for (var in log2_var) {
  elspa_median[[paste0('log2_', var)]] = log2(elspa_median[[var]] + 0.001)
}
summary(elspa_median)

#### Create eGFR groups ####
elspa_median$eGFRgrp = cut(elspa_median$CKD_EPI_1, 
                            breaks = c(-Inf, 30, 60, Inf),
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
    model = coxph(formula, data = elspa_median)
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
write.csv(gf_hr_star, 'cox_median_graft_failure.csv', row.names = TRUE)


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
    model = coxph(formula, data = elspa_median)
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
write.csv(mt_hr_star, 'cox_median_mortality.csv', row.names = TRUE)
