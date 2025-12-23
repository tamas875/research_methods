library(splines)
library(rms)
library(ggplot2)
library(ggpubr)
library(svglite)
library(patchwork)
library(rsvg)


#### Prepare data ####
d_gf = elspa_fu[, c('time_gf_fu4', 'event_gf_fu4', 'log2_N8AcSPD_24h',
                    'log2_N1AcSPD_24h', 'log2_SPD_24h', 'Age_visit_1',
                    'Geslacht', 'BMI1')]
ddist_gf = datadist(d_gf)
options(datadist = 'ddist_gf')

d_mort = elspa_fu[, c('time_mort_fu4', 'event_mort_fu4', 'log2_AcPUT_24h',
                      'log2_N8AcSPD_24h', 'log2_N1AcSPD_24h', 'log2_SPD_24h',
                      'Age_visit_1', 'Geslacht', 'BMI1')]
ddist_mort = datadist(d_mort)
options(datadist = 'ddist_mort')


#### Find optimal knots for graft failure ####
# N1-acetylspermidine
knots_range = 3:7
aic_gf1 = c()
for (k in knots_range) {
  fit = cph(Surv(time_gf_fu4, event_gf_fu4) ~ 
              rcs(log2_N1AcSPD_24h, k) + Age_visit_1 + Geslacht + BMI1,
            data = d_gf, x = TRUE)
  aic_gf1 = append(aic_gf1, AIC(fit))
}
best_k_gf1 = knots_range[which.min(aic_gf1)]

# N8-acetylspermidine
aic_gf2 = c()
for (k in knots_range) {
  fit = cph(Surv(time_gf_fu4, event_gf_fu4) ~ 
              rcs(log2_N8AcSPD_24h, k) + Age_visit_1 + Geslacht + BMI1,
            data = d_gf, x = TRUE)
  aic_gf2 = append(aic_gf2, AIC(fit))
}
best_k_gf2 = knots_range[which.min(aic_gf2)]

# Spermidine
aic_gf3 = c()
for (k in knots_range) {
  fit = cph(Surv(time_gf_fu4, event_gf_fu4) ~ 
              rcs(log2_SPD_24h, k) + Age_visit_1 + Geslacht + BMI1,
            data = d_gf, x = TRUE)
  aic_gf3 = append(aic_gf3, AIC(fit))
}
best_k_gf3 = knots_range[which.min(aic_gf3)]


#### Find optimal knots for mortality ####
# N1-acetylspermidine
aic_mort1 = c()
for (k in knots_range) {
  fit = cph(Surv(time_mort_fu4, event_mort_fu4) ~ 
              rcs(log2_N1AcSPD_24h, k) + Age_visit_1 + Geslacht + BMI1,
            data = d_mort, x = TRUE)
  aic_mort1 = append(aic_mort1, AIC(fit))
}
best_k_mort1 = knots_range[which.min(aic_mort1)]

# N8-acetylspermidine
aic_mort2 = c()
for (k in knots_range) {
  fit = cph(Surv(time_mort_fu4, event_mort_fu4) ~ 
              rcs(log2_N8AcSPD_24h, k) + Age_visit_1 + Geslacht + BMI1,
            data = d_mort, x = TRUE)
  aic_mort2 = append(aic_mort2, AIC(fit))
}
best_k_mort2 = knots_range[which.min(aic_mort2)]

# Spermidine
aic_mort3 = c()
for (k in knots_range) {
  fit = cph(Surv(time_mort_fu4, event_mort_fu4) ~ 
              rcs(log2_SPD_24h, k) + Age_visit_1 + Geslacht + BMI1,
            data = d_mort, x = TRUE)
  aic_mort3 = append(aic_mort3, AIC(fit))
}
best_k_mort3 = knots_range[which.min(aic_mort3)]

# N-acetylputrescine
aic_mort4 = c()
for (k in knots_range) {
  fit = cph(Surv(time_mort_fu4, event_mort_fu4) ~ 
              rcs(log2_AcPUT_24h, k) + Age_visit_1 + Geslacht + BMI1,
            data = d_mort, x = TRUE)
  aic_mort4 = append(aic_mort4, AIC(fit))
}
best_k_mort4 = knots_range[which.min(aic_mort4)]


#### Graft failure models ####
fit_gf1 = cph(Surv(time_gf_fu4, event_gf_fu4) ~
                rcs(log2_N1AcSPD_24h, best_k_gf1) + Age_visit_1 + Geslacht + BMI1,
              data = d_gf, x = TRUE)
pred_gf1 = Predict(fit_gf1, log2_N1AcSPD_24h, ref.zero = TRUE, fun = exp)

fit_gf2 = cph(Surv(time_gf_fu4, event_gf_fu4) ~
                rcs(log2_N8AcSPD_24h, best_k_gf2) + Age_visit_1 + Geslacht + BMI1,
              data = d_gf, x = TRUE)
pred_gf2 = Predict(fit_gf2, log2_N8AcSPD_24h, ref.zero = TRUE, fun = exp)

fit_gf3 = cph(Surv(time_gf_fu4, event_gf_fu4) ~
                rcs(log2_SPD_24h, best_k_gf3) + Age_visit_1 + Geslacht + BMI1,
              data = d_gf, x = TRUE)
pred_gf3 = Predict(fit_gf3, log2_SPD_24h, ref.zero = TRUE, fun = exp)


#### Mortality models ####
fit_mort1 = cph(Surv(time_mort_fu4, event_mort_fu4) ~
                  rcs(log2_N1AcSPD_24h, best_k_mort1) + Age_visit_1 + Geslacht + BMI1,
                data = d_mort, x = TRUE)
pred_mort1 = Predict(fit_mort1, log2_N1AcSPD_24h, ref.zero = TRUE, fun = exp)

fit_mort2 = cph(Surv(time_mort_fu4, event_mort_fu4) ~
                  rcs(log2_N8AcSPD_24h, best_k_mort2) + Age_visit_1 + Geslacht + BMI1,
                data = d_mort, x = TRUE)
pred_mort2 = Predict(fit_mort2, log2_N8AcSPD_24h, ref.zero = TRUE, fun = exp)

fit_mort3 = cph(Surv(time_mort_fu4, event_mort_fu4) ~
                  rcs(log2_SPD_24h, best_k_mort3) + Age_visit_1 + Geslacht + BMI1,
                data = d_mort, x = TRUE)
pred_mort3 = Predict(fit_mort3, log2_SPD_24h, ref.zero = TRUE, fun = exp)

fit_mort4 = cph(Surv(time_mort_fu4, event_mort_fu4) ~
                  rcs(log2_AcPUT_24h, best_k_mort4) + Age_visit_1 + Geslacht + BMI1,
                data = d_mort, x = TRUE)
pred_mort4 = Predict(fit_mort4, log2_AcPUT_24h, ref.zero = TRUE, fun = exp)


#### Create plots ####
# Plot A: N1-acetylspermidine for graft failure
plot_a = ggplot(pred_gf1) +
  geom_line(aes(x = log2_N1AcSPD_24h, y = yhat), color = '#009E73', alpha = 0.8) +
  geom_ribbon(aes(x = log2_N1AcSPD_24h, ymin = lower, ymax = upper), 
              alpha = 0.3, fill = '#009E73') +
  geom_hline(yintercept = 1, color = 'grey20', linetype = 2) +
  geom_histogram(data = d_gf, aes(x = log2_N1AcSPD_24h,
                                  y = after_stat(count) / max(after_stat(count)) * max(pred_gf1$yhat) * 1.5),
                breaks = seq(min(pred_gf1$log2_N1AcSPD_24h), max(pred_gf1$log2_N1AcSPD_24h), length.out = 20),
                fill = 'grey', alpha = 0.3, color = 'grey30') +
  labs(x = expression(paste('N'^'1', '-acetylspermidine excretion (log2 μmol/24h)')),
       y = 'HR for graft failure', title = 'A') +
  theme_pubr() +
  theme(plot.title = element_text(face = 'bold', size = 16),
        plot.subtitle = element_blank(),
        plot.caption = element_blank()) +
  scale_y_continuous(sec.axis = sec_axis(~ . / (max(pred_gf1$yhat) * 1.5) * 
                                           max(table(cut(d_gf$log2_N1AcSPD_24h, breaks = 30)), na.rm = TRUE),
                                        name = 'Frequency'))
plot_a


# Plot B: N8-acetylspermidine for graft failure
plot_b = ggplot(pred_gf2) +
  geom_line(aes(x = log2_N8AcSPD_24h, y = yhat), color = '#009E73', alpha = 0.8) +
  geom_ribbon(aes(x = log2_N8AcSPD_24h, ymin = lower, ymax = upper), 
              alpha = 0.3, fill = '#009E73') +
  geom_hline(yintercept = 1, color = 'grey20', linetype = 2) +
  geom_histogram(data = d_gf, aes(x = log2_N8AcSPD_24h,
                                  y = after_stat(count) / max(after_stat(count)) * max(pred_gf2$yhat) * 1.5),
                breaks = seq(min(pred_gf2$log2_N8AcSPD_24h), max(pred_gf2$log2_N8AcSPD_24h), length.out = 20),
                fill = 'grey', alpha = 0.3, color = 'grey30') +
  labs(x = expression(paste('N'^'8', '-acetylspermidine excretion (log2 μmol/24h)')),
       y = 'HR for graft failure', title = 'B') +
  theme_pubr() +
  theme(plot.title = element_text(face = 'bold', size = 16),
        plot.subtitle = element_blank(),
        plot.caption = element_blank()) +
  scale_y_continuous(sec.axis = sec_axis(~ . / (max(pred_gf2$yhat) * 1.5) * 
                                           max(table(cut(d_gf$log2_N8AcSPD_24h, breaks = 30)), na.rm = TRUE),
                                        name = 'Frequency'))
plot_b

# Plot C: Spermidine for graft failure
plot_c = ggplot(pred_gf3) +
  geom_line(aes(x = log2_SPD_24h, y = yhat), color = '#009E73', alpha = 0.8) +
  geom_ribbon(aes(x = log2_SPD_24h, ymin = lower, ymax = upper), 
              alpha = 0.3, fill = '#009E73') +
  geom_hline(yintercept = 1, color = 'grey20', linetype = 2) +
  geom_histogram(data = d_gf, aes(x = log2_SPD_24h,
                                  y = after_stat(count) / max(after_stat(count)) * max(pred_gf3$yhat) * 1.5),
                breaks = seq(min(pred_gf3$log2_SPD_24h), max(pred_gf3$log2_SPD_24h), length.out = 20),
                fill = 'grey', alpha = 0.3, color = 'grey30') +
  labs(x = 'Spermidine excretion (log2 μmol/24h)',
       y = 'HR for graft failure', title = 'C') +
  theme_pubr() +
  theme(plot.title = element_text(face = 'bold', size = 16),
        plot.subtitle = element_blank(),
        plot.caption = element_blank()) +
  scale_y_continuous(sec.axis = sec_axis(~ . / (max(pred_gf3$yhat) * 1.5) * 
                                           max(table(cut(d_gf$log2_SPD_24h, breaks = 30)), na.rm = TRUE),
                                        name = 'Frequency'))
plot_c

# Plot D: N1-acetylspermidine for mortality
plot_d = ggplot(pred_mort1) +
  geom_line(aes(x = log2_N1AcSPD_24h, y = yhat), color = '#d55e00', alpha = 0.8) +
  geom_ribbon(aes(x = log2_N1AcSPD_24h, ymin = lower, ymax = upper), 
              alpha = 0.3, fill = '#d55e00') +
  geom_hline(yintercept = 1, color = 'grey20', linetype = 2) +
  geom_histogram(data = d_mort, aes(x = log2_N1AcSPD_24h,
                                    y = after_stat(count) / max(after_stat(count)) * max(pred_mort1$yhat) * 1.5),
                breaks = seq(min(pred_mort1$log2_N1AcSPD_24h), max(pred_mort1$log2_N1AcSPD_24h), length.out = 20),
                fill = 'grey', alpha = 0.3, color = 'grey30') +
  labs(x = expression(paste('N'^'1', '-acetylspermidine excretion (log2 μmol/24h)')),
       y = 'HR for mortality', title = 'D') +
  theme_pubr() +
  theme(plot.title = element_text(face = 'bold', size = 16),
        plot.subtitle = element_blank(),
        plot.caption = element_blank()) +
  scale_y_continuous(sec.axis = sec_axis(~ . / (max(pred_mort1$yhat) * 1.5) * 
                                           max(table(cut(d_mort$log2_N1AcSPD_24h, breaks = 30)), na.rm = TRUE),
                                        name = 'Frequency'))
plot_d

# Plot E: N8-acetylspermidine for mortality
plot_e = ggplot(pred_mort2) +
  geom_line(aes(x = log2_N8AcSPD_24h, y = yhat), color = '#d55e00', alpha = 0.8) +
  geom_ribbon(aes(x = log2_N8AcSPD_24h, ymin = lower, ymax = upper), 
              alpha = 0.3, fill = '#d55e00') +
  geom_hline(yintercept = 1, color = 'grey20', linetype = 2) +
  geom_histogram(data = d_mort, aes(x = log2_N8AcSPD_24h,
                                    y = after_stat(count) / max(after_stat(count)) * max(pred_mort2$yhat) * 1.5),
                breaks = seq(min(pred_mort2$log2_N8AcSPD_24h), max(pred_mort2$log2_N8AcSPD_24h), length.out = 20),
                fill = 'grey', alpha = 0.3, color = 'grey30') +
  labs(x = expression(paste('N'^'8', '-acetylspermidine excretion (log2 μmol/24h)')),
       y = 'HR for mortality', title = 'E') +
  theme_pubr() +
  theme(plot.title = element_text(face = 'bold', size = 16),
        plot.subtitle = element_blank(),
        plot.caption = element_blank()) +
  scale_y_continuous(sec.axis = sec_axis(~ . / (max(pred_mort2$yhat) * 1.5) * 
                                           max(table(cut(d_mort$log2_N8AcSPD_24h, breaks = 30)), na.rm = TRUE),
                                        name = 'Frequency'))
plot_e

# Plot F: Spermidine for mortality
plot_f = ggplot(pred_mort3) +
  geom_line(aes(x = log2_SPD_24h, y = yhat), color = '#d55e00', alpha = 0.8) +
  geom_ribbon(aes(x = log2_SPD_24h, ymin = lower, ymax = upper), 
              alpha = 0.3, fill = '#d55e00') +
  geom_hline(yintercept = 1, color = 'grey20', linetype = 2) +
  geom_histogram(data = d_mort, aes(x = log2_SPD_24h,
                                    y = after_stat(count) / max(after_stat(count)) * max(pred_mort3$yhat) * 1.5),
                breaks = seq(min(pred_mort3$log2_SPD_24h), max(pred_mort3$log2_SPD_24h), length.out = 20),
                fill = 'grey', alpha = 0.3, color = 'grey30') +
  labs(x = 'Spermidine excretion (log2 μmol/24h)',
       y = 'HR for mortality', title = 'F') +
  theme_pubr() +
  theme(plot.title = element_text(face = 'bold', size = 16),
        plot.subtitle = element_blank(),
        plot.caption = element_blank()) +
  scale_y_continuous(sec.axis = sec_axis(~ . / (max(pred_mort3$yhat) * 1.5) * 
                                           max(table(cut(d_mort$log2_SPD_24h, breaks = 30)), na.rm = TRUE),
                                        name = 'Frequency'))
plot_f

# Plot G: N-acetylputrescine for mortality
plot_g = ggplot(pred_mort4) +
  geom_line(aes(x = log2_AcPUT_24h, y = yhat), color = '#d55e00', alpha = 0.8) +
  geom_ribbon(aes(x = log2_AcPUT_24h, ymin = lower, ymax = upper), 
              alpha = 0.3, fill = '#d55e00') +
  geom_hline(yintercept = 1, color = 'grey20', linetype = 2) +
  geom_histogram(data = d_mort, aes(x = log2_AcPUT_24h,
                                    y = after_stat(count) / max(after_stat(count)) * max(pred_mort4$yhat) * 1.5),
                breaks = seq(min(pred_mort4$log2_AcPUT_24h), max(pred_mort4$log2_AcPUT_24h), length.out = 20),
                fill = 'grey', alpha = 0.3, color = 'grey30') +
  labs(x = 'N-acetylputrescine excretion (log2 μmol/24h)',
       y = 'HR for mortality', title = 'G') +
  theme_pubr() +
  theme(plot.title = element_text(face = 'bold', size = 16),
        plot.subtitle = element_blank(),
        plot.caption = element_blank()) +
  scale_y_continuous(sec.axis = sec_axis(~ . / (max(pred_mort4$yhat) * 1.5) * 
                                           max(table(cut(d_mort$log2_AcPUT_24h, breaks = 30)), na.rm = TRUE),
                                        name = 'Frequency'))
plot_g

#### Combine and save plots ####
empty_plot = ggplot() + theme_void()
combined_plot = wrap_plots(
  plot_a, plot_d,
  plot_b, plot_e,
  plot_c, plot_f,
  empty_plot, plot_g,
  ncol = 2, widths = c(1, 1), heights = rep(1, 4)
)

ggsave('combined_splines_plot.svg', plot = combined_plot,
       width = 15, height = 24, units = 'in', dpi = 300, device = svglite)
rsvg_pdf('combined_splines_plot.svg', 'combined_splines_plot.pdf')
