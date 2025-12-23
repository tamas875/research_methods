library(lm.beta)
library(car)
library(ggplot2)


#### Convert to log2 scale ####
names(elspa_fu)
log2_var = names(elspa_fu)[c(7, 26:40, 59:65)]
for (var in log2_var) {
  elspa_fu[[paste0('log2_', var)]] = log2(elspa_fu[[var]] + 0.001)
}


#### Define variables ####
names(elspa_fu)
predictors = names(elspa_fu)[69:83]
outcomes = c('log2_AcPUT_24h', 'log2_N1AcSPD_24h', 'log2_N8AcSPD_24h',
             'log2_CAD_24h', 'log2_PUT_24h', 'log2_SPM_24h', 'log2_SPD_24h')
adjusters = c('Age_visit_1', 'Geslacht', 'BMI1', 'CKD_EPI_1')
outcome_names = gsub('log2_|_24h', '', outcomes)


#### Univariable regression ####
n_pred = length(predictors)
n_out = length(outcomes)
univ_beta = matrix(NA, nrow = n_pred, ncol = n_out)
univ_p = matrix(NA, nrow = n_pred, ncol = n_out)
rownames(univ_beta) = predictors
rownames(univ_p) = predictors
colnames(univ_beta) = outcome_names
colnames(univ_p) = outcome_names

for (i in seq_along(outcomes)) {
  outcome = outcomes[i]
  
  for (j in seq_along(predictors)) {
    predictor = predictors[j]
    
    # Standardize data
    is_cat = is.factor(elspa_fu[[predictor]]) || is.character(elspa_fu[[predictor]])
    if (is_cat) {
      formula = as.formula(paste(outcome, '~', predictor))
      model = lm(formula, data = elspa_fu)
    } else {
      data_std = data.frame(outcome_var = scale(elspa_fu[[outcome]])[, 1],
                            predictor_var = scale(elspa_fu[[predictor]])[, 1])
      formula = as.formula('outcome_var ~ predictor_var')
      model = lm(formula, data = data_std)
    }
    
    # Extract results
    model_summary = summary(model)
    coef_table = model_summary$coefficients

    if (is_cat) {
      coef = coef_table[2, 1]
      se = coef_table[2, 2]
      p_val = coef_table[2, 4]
    } else {
      coef = coef_table['predictor_var', 1]
      se = coef_table['predictor_var', 2]
      p_val = coef_table['predictor_var', 4]
    }
    
    ci_low = coef - 1.96 * se
    ci_high = coef + 1.96 * se
    
    univ_beta[j, i] = sprintf('%.2f (%.2f, %.2f)', coef, ci_low, ci_high)
    univ_p[j, i] = p_val
  }
}


#### Multivariable regression ####
multi_beta = matrix(NA, nrow = n_pred, ncol = n_out)
multi_p = matrix(NA, nrow = n_pred, ncol = n_out)
vif_mat = matrix(NA, nrow = n_pred, ncol = n_out)
rownames(multi_beta) = predictors
rownames(multi_p) = predictors
rownames(vif_mat) = predictors
colnames(multi_beta) = outcome_names
colnames(multi_p) = outcome_names
colnames(vif_mat) = outcome_names

for (i in seq_along(outcomes)) {
  outcome = outcomes[i]
  
  for (j in seq_along(predictors)) {
    predictor = predictors[j]
    
    # Standardize data
    is_cat = is.factor(elspa_fu[[predictor]]) || is.character(elspa_fu[[predictor]])
    data_list = list()
    data_list$outcome_var = scale(elspa_fu[[outcome]])[, 1]
    if (is_cat) {
      data_list$predictor_var = elspa_fu[[predictor]]
    } else {
      data_list$predictor_var = scale(elspa_fu[[predictor]])[, 1]
    }
    
    # Add adjusters
    for (adj in adjusters) {
      if (is.numeric(elspa_fu[[adj]]) && !is.factor(elspa_fu[[adj]])) {
        data_list[[adj]] = scale(elspa_fu[[adj]])[, 1]
      } else {
        data_list[[adj]] = elspa_fu[[adj]]
      }
    }
    data_std = as.data.frame(data_list)
    
    # Fit model
    formula = as.formula(paste('outcome_var ~', paste(c('predictor_var', adjusters), collapse = ' + ')))
    model = lm(formula, data = data_std)
    
    # Extract results
    model_summary = summary(model)
    coef_table = model_summary$coefficients
    
    # Get coefficient and CI
    if (is_cat) {
      coef = coef_table[2, 1]
      se = coef_table[2, 2]
      p_val = coef_table[2, 4]
    } else {
      coef = coef_table['predictor_var', 1]
      se = coef_table['predictor_var', 2]
      p_val = coef_table['predictor_var', 4]
    }
    
    ci_low = coef - 1.96 * se
    ci_high = coef + 1.96 * se
    
    multi_beta[j, i] = sprintf('%.2f (%.2f, %.2f)', coef, ci_low, ci_high)
    multi_p[j, i] = p_val
    
    # Calculate VIF
    if (!is_cat) {
      vif_val = car::vif(model)
      vif_mat[j, i] = round(vif_val['predictor_var'], 2)
    }
  }
}


#### Add significance stars ####
# Adjust p-values and add stars to univariable results
univ_beta_star = univ_beta
for (j in 1:ncol(univ_p)) {
  p_adj = p.adjust(univ_p[, j], method = 'BH')
  
  for (i in 1:nrow(univ_p)) {
    if (!is.na(p_adj[i])) {
      if (p_adj[i] < 0.001) {
        univ_beta_star[i, j] = paste0(univ_beta[i, j], '***')
      } else if (p_adj[i] < 0.01) {
        univ_beta_star[i, j] = paste0(univ_beta[i, j], '**')
      } else if (p_adj[i] < 0.05) {
        univ_beta_star[i, j] = paste0(univ_beta[i, j], '*')
      }
    }
  }
}

# Adjust p-values and add stars to multivariable results
multi_beta_star = multi_beta
for (j in 1:ncol(multi_p)) {
  p_adj = p.adjust(multi_p[, j], method = 'BH')
  
  for (i in 1:nrow(multi_p)) {
    if (!is.na(p_adj[i])) {
      if (p_adj[i] < 0.001) {
        multi_beta_star[i, j] = paste0(multi_beta[i, j], '***')
      } else if (p_adj[i] < 0.01) {
        multi_beta_star[i, j] = paste0(multi_beta[i, j], '**')
      } else if (p_adj[i] < 0.05) {
        multi_beta_star[i, j] = paste0(multi_beta[i, j], '*')
      }
    }
  }
}
print(multi_beta_star)
print(vif_mat)
write.csv(multi_beta_star, 'multi_beta_significant.csv', row.names = TRUE)

