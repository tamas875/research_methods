library(dplyr)
library(tidyr)
library(ggplot2)
library(rstatix)
library(ggpubr)
library(svglite)
library(patchwork)
library(rsvg)


#### Define polyamine variables ####
pa_vars = c('AcPUT_24h', 'N1AcSPD_24h', 'N8AcSPD_24h',
            'CAD_24h', 'PUT_24h', 'SPM_24h', 'SPD_24h', 'totalPA_24h')

metabolites = c('Total acetylated polyamines',
                'N-acetylputrescine',
                'N1-acetylspermidine', 
                'N8-acetylspermidine',
                'Cadaverine',
                'Putrescine',
                'Spermine',
                'Spermidine',
                'Total polyamines')


#### Prepare datasets ####
donor$totalAcetylPA_24h = donor$AcPUT_24h + donor$N1AcSPD_24h + donor$N8AcSPD_24h
elspa$totalAcetylPA_24h = elspa$AcPUT_24h + elspa$N1AcSPD_24h + elspa$N8AcSPD_24h
pa_vars = c('totalAcetylPA_24h', pa_vars)

donor_data = donor %>%
  select(Subjectnr, all_of(pa_vars)) %>%
  mutate(group = 'Donors', gf_group = 'Donors')

elspa_data = elspa %>%
  mutate(group = ifelse(event_mort_fu4 == 1, 'Recipients (Death)', 'Recipients (Alive)'),
         gf_group = ifelse(event_gf_fu4 == 1, 'Recipients (Graft Failure)', 'Recipients (Functioning Graft)')) %>%
  select(Subjectnr, all_of(pa_vars), group, gf_group)

all_data = rbind(donor_data, elspa_data)


#### Table 2: Donors vs Recipients ####
table2_data = all_data %>% 
  mutate(group = ifelse(group == 'Donors', 'Donors', 'Recipients'))

table2 = data.frame(Metabolite = metabolites)

for (group_name in c('Donors', 'Recipients')) {
  group_data = table2_data[table2_data$group == group_name, ]
  
  mean_vals = sapply(pa_vars, function(x) mean(group_data[[x]], na.rm = TRUE))
  sd_vals = sapply(pa_vars, function(x) sd(group_data[[x]], na.rm = TRUE))
  median_vals = sapply(pa_vars, function(x) median(group_data[[x]], na.rm = TRUE))
  q1_vals = sapply(pa_vars, function(x) quantile(group_data[[x]], 0.25, na.rm = TRUE))
  q3_vals = sapply(pa_vars, function(x) quantile(group_data[[x]], 0.75, na.rm = TRUE))
  
  table2[[paste0(group_name, '_Mean_SD')]] = sprintf('%.2f ± %.2f', mean_vals, sd_vals)
  table2[[paste0(group_name, '_Median_IQR')]] = sprintf('%.2f [%.2f-%.2f]', median_vals, q1_vals, q3_vals)
}

mw_p = c()
for (var in pa_vars) {
  p_val = wilcox.test(as.formula(paste(var, '~ group')), data = table2_data)$p.value
  mw_p = append(mw_p, ifelse(p_val < 0.001, '<0.001', sprintf('%.3f', p_val)))
}
table2$Mann_Whitney_P = mw_p

print(table2)
write.csv(table2, 'Table2_Donor_vs_Recipients.csv', row.names = FALSE)


#### Tables 3 & 4: Summary statistics ####
table3 = data.frame(Metabolite = metabolites)
table4 = data.frame(Metabolite = metabolites)

for (group_name in c('Donors', 'Recipients (Alive)', 'Recipients (Death)')) {
  group_data = all_data[all_data$group == group_name, ]
  
  mean_vals = sapply(pa_vars, function(x) mean(group_data[[x]], na.rm = TRUE))
  sd_vals = sapply(pa_vars, function(x) sd(group_data[[x]], na.rm = TRUE))
  median_vals = sapply(pa_vars, function(x) median(group_data[[x]], na.rm = TRUE))
  q1_vals = sapply(pa_vars, function(x) quantile(group_data[[x]], 0.25, na.rm = TRUE))
  q3_vals = sapply(pa_vars, function(x) quantile(group_data[[x]], 0.75, na.rm = TRUE))
  
  col_prefix = gsub(' ', '_', group_name)
  table3[[paste0(col_prefix, '_Mean_SD')]] = sprintf('%.2f ± %.2f', mean_vals, sd_vals)
  table3[[paste0(col_prefix, '_Median_IQR')]] = sprintf('%.2f [%.2f-%.2f]', median_vals, q1_vals, q3_vals)
}

for (group_name in c('Donors', 'Recipients (Functioning Graft)', 'Recipients (Graft Failure)')) {
  group_data = all_data[all_data$gf_group == group_name, ]
  
  mean_vals = sapply(pa_vars, function(x) mean(group_data[[x]], na.rm = TRUE))
  sd_vals = sapply(pa_vars, function(x) sd(group_data[[x]], na.rm = TRUE))
  median_vals = sapply(pa_vars, function(x) median(group_data[[x]], na.rm = TRUE))
  q1_vals = sapply(pa_vars, function(x) quantile(group_data[[x]], 0.25, na.rm = TRUE))
  q3_vals = sapply(pa_vars, function(x) quantile(group_data[[x]], 0.75, na.rm = TRUE))
  
  col_prefix = gsub(' ', '_', group_name)
  table4[[paste0(col_prefix, '_Mean_SD')]] = sprintf('%.2f ± %.2f', mean_vals, sd_vals)
  table4[[paste0(col_prefix, '_Median_IQR')]] = sprintf('%.2f [%.2f-%.2f]', median_vals, q1_vals, q3_vals)
}

print(table3)
print(table4)


#### Tables 5 & 6: Statistical tests ####
table5 = data.frame(Metabolite = metabolites)
table6 = data.frame(Metabolite = metabolites)

# Table 5: Mortality groups
kw_p = c()
for (var in pa_vars) {
  kw_result = kruskal.test(as.formula(paste(var, '~ group')), data = all_data)
  kw_p = append(kw_p, ifelse(kw_result$p.value < 0.001, '<0.001', sprintf('%.3f', kw_result$p.value)))
}
table5$Overall_P = kw_p

comparisons = list(c('Donors', 'Recipients (Alive)'),
                   c('Donors', 'Recipients (Death)'),
                   c('Recipients (Alive)', 'Recipients (Death)'))

for (comp in comparisons) {
  col_name = paste0(gsub(' |\\(|\\)', '', comp), collapse = '_vs_')
  dunn_p = c()
  
  for (var in pa_vars) {
    dunn_result = all_data %>%
      dunn_test(as.formula(paste(var, '~ group')), p.adjust.method = 'bonferroni')
    
    idx = which((dunn_result$group1 == comp[1] & dunn_result$group2 == comp[2]) |
                (dunn_result$group1 == comp[2] & dunn_result$group2 == comp[1]))
    
    if (length(idx) > 0) {
      p_val = dunn_result$p.adj[idx]
      dunn_p = append(dunn_p, ifelse(p_val < 0.001, '<0.001', 
                                     ifelse(p_val >= 0.05, 'NS', sprintf('%.3f', p_val))))
    } else {
      dunn_p = append(dunn_p, 'NS')
    }
  }
  table5[[col_name]] = dunn_p
}

# Table 6: Graft failure groups
kw_p = c()
for (var in pa_vars) {
  kw_result = kruskal.test(as.formula(paste(var, '~ gf_group')), data = all_data)
  kw_p = append(kw_p, ifelse(kw_result$p.value < 0.001, '<0.001', sprintf('%.3f', kw_result$p.value)))
}
table6$Overall_P = kw_p

comparisons = list(c('Donors', 'Recipients (Functioning Graft)'),
                   c('Donors', 'Recipients (Graft Failure)'),
                   c('Recipients (Functioning Graft)', 'Recipients (Graft Failure)'))

for (comp in comparisons) {
  col_name = paste0(gsub(' |\\(|\\)', '', comp), collapse = '_vs_')
  dunn_p = c()
  
  for (var in pa_vars) {
    dunn_result = all_data %>%
      dunn_test(as.formula(paste(var, '~ gf_group')), p.adjust.method = 'bonferroni')
    
    idx = which((dunn_result$group1 == comp[1] & dunn_result$group2 == comp[2]) |
                (dunn_result$group1 == comp[2] & dunn_result$group2 == comp[1]))
    
    if (length(idx) > 0) {
      p_val = dunn_result$p.adj[idx]
      dunn_p = append(dunn_p, ifelse(p_val < 0.001, '<0.001', 
                                     ifelse(p_val >= 0.05, 'NS', sprintf('%.3f', p_val))))
    } else {
      dunn_p = append(dunn_p, 'NS')
    }
  }
  table6[[col_name]] = dunn_p
}

print(table5)
print(table6)


#### Generate violin plots ####
plot_data_mort = all_data
plot_data_mort$group = factor(plot_data_mort$group, 
                              levels = c('Donors', 'Recipients (Alive)', 'Recipients (Death)'))
plot_data_gf = all_data
plot_data_gf$gf_group = factor(plot_data_gf$gf_group, 
                               levels = c('Donors', 'Recipients (Functioning Graft)', 'Recipients (Graft Failure)'))

plots_mort = list()
plots_gf = list()

for (i in seq_along(pa_vars)) {
  var = pa_vars[i]
  met_name = metabolites[i]
  
  # Mortality plot
  stats_mort = table5[i, ]
  y_max_mort = max(plot_data_mort[[var]], na.rm = TRUE)
  y_range_mort = diff(range(plot_data_mort[[var]], na.rm = TRUE))
  
  p_mort = ggplot(plot_data_mort, aes(x = group, y = .data[[var]], fill = group)) +
    geom_violin(trim = FALSE, alpha = 0.6) +
    geom_boxplot(width = 0.15, fill = 'white', alpha = 0.4, outlier.size = 0.5) +
    scale_fill_manual(values = c('Donors' = '#2C699A', 
                                  'Recipients (Alive)' = '#54B435', 
                                  'Recipients (Death)' = '#E74646')) +
    scale_x_discrete(labels = c('Donors', 'Recipients\n(Alive)', 'Recipients\n(Death)')) +
    labs(x = '', y = 'Excretion (μmol/24h)') +
    theme_classic() +
    theme(plot.title = element_text(size = 12, face = 'bold', hjust = 0.5),
          axis.text = element_text(size = 10, color = 'black'),
          axis.title.y = element_text(size = 10),
          legend.position = 'none')
  
  # Graft failure plot
  stats_gf = table6[i, ]
  y_max_gf = max(plot_data_gf[[var]], na.rm = TRUE)
  y_range_gf = diff(range(plot_data_gf[[var]], na.rm = TRUE))
  
  p_gf = ggplot(plot_data_gf, aes(x = gf_group, y = .data[[var]], fill = gf_group)) +
    geom_violin(trim = FALSE, alpha = 0.6) +
    geom_boxplot(width = 0.15, fill = 'white', alpha = 0.4, outlier.size = 0.5) +
    scale_fill_manual(values = c('Donors' = '#2C699A', 
                                  'Recipients (Functioning Graft)' = '#8B4789', 
                                  'Recipients (Graft Failure)' = '#E07B39')) +
    scale_x_discrete(labels = c('Donors', 'Recipients\n(Functioning Graft)', 'Recipients\n(Graft Failure)')) +
    labs(x = '', y = 'Excretion (μmol/24h)') +
    theme_classic() +
    theme(plot.title = element_text(size = 12, face = 'bold', hjust = 0.5),
          axis.text = element_text(size = 10, color = 'black'),
          axis.title.y = element_text(size = 10),
          legend.position = 'none')
  
  # Add titles
  if (met_name == 'N1-acetylspermidine') {
    p_mort = p_mort + ggtitle(expression(bold(paste('N'^'1', '-acetylspermidine'))))
    p_gf = p_gf + ggtitle(expression(bold(paste('N'^'1', '-acetylspermidine'))))
  } else if (met_name == 'N8-acetylspermidine') {
    p_mort = p_mort + ggtitle(expression(bold(paste('N'^'8', '-acetylspermidine'))))
    p_gf = p_gf + ggtitle(expression(bold(paste('N'^'8', '-acetylspermidine'))))
  } else {
    p_mort = p_mort + ggtitle(met_name)
    p_gf = p_gf + ggtitle(met_name)
  }
  
  # Add Kruskal-Wallis annotations
  kw_label_mort = paste('Kruskal-Wallis, p', 
                        ifelse(stats_mort$Overall_P == '<0.001', '< 0.001', paste('=', stats_mort$Overall_P)))
  p_mort = p_mort + annotate('text', x = 1, y = y_max_mort * 1.33, label = kw_label_mort, size = 3)
  
  kw_label_gf = paste('Kruskal-Wallis, p', 
                      ifelse(stats_gf$Overall_P == '<0.001', '< 0.001', paste('=', stats_gf$Overall_P)))
  p_gf = p_gf + annotate('text', x = 1, y = y_max_gf * 1.33, label = kw_label_gf, size = 3)
  
  # Add pairwise comparisons for mortality
  dunn_results_mort = all_data %>%
    dunn_test(as.formula(paste(var, '~ group')), p.adjust.method = 'bonferroni')
  dunn_results_mort$p.adj.signif = cut(dunn_results_mort$p.adj, 
                                       breaks = c(-Inf, 0.001, 0.01, 0.05, Inf), 
                                       labels = c('***', '**', '*', 'ns'))
  dunn_sig_mort = dunn_results_mort[dunn_results_mort$p.adj.signif != 'ns', ]
  
  if (nrow(dunn_sig_mort) > 0) {
    dunn_sig_mort$y.position = seq(y_max_mort + y_range_mort * 0.05, y_max_mort + y_range_mort * 0.20, 
                                   length.out = nrow(dunn_sig_mort))
    p_mort = p_mort + stat_pvalue_manual(dunn_sig_mort, label = 'p.adj.signif', 
                                         tip.length = 0.02, bracket.nudge.y = 0)
  }
  
  # Add pairwise comparisons for graft failure
  dunn_results_gf = all_data %>%
    dunn_test(as.formula(paste(var, '~ gf_group')), p.adjust.method = 'bonferroni')
  dunn_results_gf$p.adj.signif = cut(dunn_results_gf$p.adj, 
                                     breaks = c(-Inf, 0.001, 0.01, 0.05, Inf), 
                                     labels = c('***', '**', '*', 'ns'))
  dunn_sig_gf = dunn_results_gf[dunn_results_gf$p.adj.signif != 'ns', ]
  
  if (nrow(dunn_sig_gf) > 0) {
    dunn_sig_gf$y.position = seq(y_max_gf + y_range_gf * 0.05, y_max_gf + y_range_gf * 0.20, 
                                 length.out = nrow(dunn_sig_gf))
    p_gf = p_gf + stat_pvalue_manual(dunn_sig_gf, label = 'p.adj.signif', 
                                     tip.length = 0.02, bracket.nudge.y = 0)
  }
  
  plots_mort[[i]] = p_mort
  plots_gf[[i]] = p_gf
}

# Create legends
mort_colors = c('Donors' = '#2C699A', 
                'Recipients (Alive)' = '#54B435', 
                'Recipients (Death)' = '#E74646')
gf_colors = c('Donors' = '#2C699A', 
              'Recipients (Functioning Graft)' = '#8B4789', 
              'Recipients (Graft Failure)' = '#E07B39')

legend_data_mort = plot_data_mort %>%
  group_by(group) %>%
  summarise(n = n(), .groups = 'drop')

legend_data_gf = plot_data_gf %>%
  group_by(gf_group) %>%
  summarise(n = n(), .groups = 'drop')

# Create legend plots
legend_plot_mort = ggplot(legend_data_mort, aes(x = 1, y = 1, fill = group)) +
  geom_tile(alpha = 0.6) +
  scale_fill_manual(values = mort_colors,
                    labels = paste0(legend_data_mort$group, ' (n=', legend_data_mort$n, ')')) +
  guides(fill = guide_legend(title = NULL)) +
  theme_void() +
  theme(legend.text = element_text(size = 10),
        legend.key.size = unit(0.5, 'cm'))

legend_plot_gf = ggplot(legend_data_gf, aes(x = 1, y = 1, fill = gf_group)) +
  geom_tile(alpha = 0.6) +
  scale_fill_manual(values = gf_colors,
                    labels = paste0(legend_data_gf$gf_group, ' (n=', legend_data_gf$n, ')')) +
  guides(fill = guide_legend(title = NULL)) +
  theme_void() +
  theme(legend.text = element_text(size = 10),
        legend.key.size = unit(0.5, 'cm'))

legend_mort = get_legend(legend_plot_mort)
legend_gf = get_legend(legend_plot_gf)

# Combine plots
plots_grid_mort = wrap_plots(plots_mort[c(9, 1, 2, 3, 4, 5, 6, 7, 8)], ncol = 3)
final_plot_mort = (plot_spacer() | as_ggplot(legend_mort)) / plots_grid_mort +
  plot_layout(heights = c(0.1, 0.9), widths = c(0.7, 0.3))

plots_grid_gf = wrap_plots(plots_gf[c(9, 1, 2, 3, 4, 5, 6, 7, 8)], ncol = 3)
final_plot_gf = (plot_spacer() | as_ggplot(legend_gf)) / plots_grid_gf +
  plot_layout(heights = c(0.1, 0.9), widths = c(0.7, 0.3))

# Save plots
ggsave('vio_plot1.svg', final_plot_mort, width = 15, height = 15, units = 'in', dpi = 300, device = svglite)
rsvg_pdf('vio_plot1.svg', 'vio_plot1.pdf')

ggsave('vio_plot2.svg', final_plot_gf, width = 15, height = 15, units = 'in', dpi = 300, device = svglite)
rsvg_pdf('vio_plot2.svg', 'vio_plot2.pdf')
