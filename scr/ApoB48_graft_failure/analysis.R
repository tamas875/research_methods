#Set up data for analysis

library(foreign)

d.aiko <- read.spss('~/Data/transplant_lines/aiko_cohort/Aiko cohort.sav', use.value.labels = TRUE, to.data.frame = TRUE)
d.apob48 <- read.csv('~/Data/transplant_lines/aiko_cohort/apob48.csv')

colnames(d.apob48) <- c('ki_num', 'b_apob48_fu1')

d.aiko <- merge(d.aiko, d.apob48, by = 'ki_num')
d.aiko <- d.aiko[which(is.na(d.aiko['CKD_EPI']) == FALSE & is.na(d.aiko['graft_failure_2012']) == FALSE & is.na(d.aiko['b_apob48_fu1']) == FALSE),]
rm(d.apob48)

#Survival analysis

library(survival)
library(CoxDisp)
library(DescribeData)

#Make tertiles of apob48
if (FALSE) {
  #Thus is outdated code
  split.data <- split(d.aiko, d.aiko$sexe_pat)
  tertiles.apob48 <- lapply(split.data, function(group) {
    group$tertiles <- cut(group$b_apob48_fu1, breaks = quantile(group$b_apob48_fu1, probs = c(0, 1/3, 2/3, 1)), labels = c("T1", "T2", "T3"))
    return(group)
    d.aiko <- do.call(rbind, tertiles.apob48)
    rm(split.data, tertiles.apob48)
  })
}
if (TRUE) {
  d.aiko$tertiles <- cut(d.aiko$b_apob48_fu1, breaks = c(-Inf, quantile(d.aiko$b_apob48_fu1, probs = 1/3), quantile(d.aiko$b_apob48_fu1, probs = 2/3), Inf), labels = c("Low tertile", "Medium tertile", "High tertile"))
}

d.aiko$c_apob48_log2 <- log2(d.aiko$b_apob48_fu1)

d.aiko$egfr_strat <- cut(d.aiko$CKD_EPI, 
                         c(-Inf, median(d.aiko$CKD_EPI, na.rm = TRUE), Inf), 
                         labels = c(paste0('< ', 
                                           round(median(d.aiko$CKD_EPI), 1)), 
                                    paste0('> ', 
                                           round(median(d.aiko$CKD_EPI), 1))))
d.aiko$hdl_strat <- cut(d.aiko$con_hdlc, 
                        c(-Inf, median(d.aiko$con_hdlc, na.rm = TRUE), Inf), 
                        labels = c(paste0('< ', 
                                          round(median(d.aiko$con_hdlc), 1)), 
                                   paste0('> ', 
                                          round(median(d.aiko$con_hdlc), 1))))
d.aiko$tgl_strat <- cut(d.aiko$con_tg, 
                        c(-Inf, median(d.aiko$con_tg, na.rm = TRUE), Inf), 
                        labels = c(paste0('< ', 
                                          round(median(d.aiko$con_tg), 1)), 
                                   paste0('> ', 
                                          round(median(d.aiko$con_tg), 1))))
d.aiko$waist_strat <- cut(d.aiko$gem_tail, 
                          c(-Inf, median(d.aiko$gem_tail, na.rm = TRUE), Inf), 
                          labels = c(paste0('< ', 
                                            round(median(d.aiko$gem_tail, na.rm = TRUE), 1)), 
                                     paste0('> ', 
                                            round(median(d.aiko$gem_tail, na.rm = TRUE), 1))))
d.aiko$tc_strat <- cut(d.aiko$con_tc, 
                       c(-Inf, median(d.aiko$con_tc, na.rm = TRUE), Inf), 
                       labels = c(paste0('< ', 
                                         round(median(d.aiko$con_tc), 1)), 
                                  paste0('> ', 
                                         round(median(d.aiko$con_tc), 1))))
d.aiko$apoa1_strat <- cut(d.aiko$con_apoa, 
                          c(-Inf, median(d.aiko$con_apoa, na.rm = TRUE), Inf), 
                          labels = c(paste0('< ', 
                                            round(median(d.aiko$con_apoa), 1)), 
                                     paste0('> ', 
                                            round(median(d.aiko$con_apoa), 1))))
d.aiko$apob_strat <- cut(d.aiko$con_apob, 
                         c(-Inf, median(d.aiko$con_apob, na.rm = TRUE), Inf), 
                         labels = c(paste0('< ', 
                                           round(median(d.aiko$con_apob), 1)), 
                                    paste0('> ', 
                                           round(median(d.aiko$con_apob), 1))))
d.aiko$time_tx_strat <- cut(d.aiko$tijd_bas, c(-Inf, 
                                               median(d.aiko$tijd_bas, na.rm = TRUE), Inf),
                            labels = c(paste0('< ', 
                                              round(median(d.aiko$tijd_bas), 1)), 
                                       paste0('> ', 
                                              round(median(d.aiko$tijd_bas), 1))))
d.aiko$hba1c_strat <- cut(d.aiko$per_hba1, c(-Inf, 
                                             median(d.aiko$per_hba1, na.rm = TRUE), Inf),
                          labels = c(paste0('< ', 
                                            round(median(d.aiko$per_hba1, na.rm = TRUE), 1)), 
                                     paste0('> ', 
                                            round(median(d.aiko$per_hba1, na.rm = TRUE), 1))))
d.aiko$diabetes_fu1 <- as.factor(ifelse(d.aiko$fvg_diab  == 'no                                                          ', 'No', 'Yes'))
d.aiko$statin_use_fu1 <- as.factor(ifelse(d.aiko$type_sta == 'geen', 'No', 'Yes'))
d.aiko$prot_urea <- ifelse(d.aiko$exc_pro == 0, 'No', 'Yes')
levels(d.aiko$sexe_pat)[which(levels(d.aiko$sexe_pat) == 'vrouw')] <- 'Female'
levels(d.aiko$sexe_pat)[which(levels(d.aiko$sexe_pat) == 'man')] <- 'Male' 

levels(d.aiko$di_sig)[which(levels(d.aiko$di_sig) == 'geen roken')] <- 'No'
levels(d.aiko$di_sig)[which(levels(d.aiko$di_sig) == 'wel roken')] <- 'Yes'

levels(d.aiko$type_cni)[which(levels(d.aiko$type_cni) == 'geen')] <- 'None'
levels(d.aiko$type_cni)[which(levels(
  d.aiko$type_cni) == 'Ciclosporine (N')] <- 'Cyclosporine'
levels(d.aiko$type_cni)[which(levels(
  d.aiko$type_cni) == 'Tacrolimus (Pro')] <- 'Tacrolimus'

levels(d.aiko$type_plr)[which(levels(d.aiko$type_plr) == 'geen')] <- 'No'
levels(d.aiko$type_plr)[which(levels(
  d.aiko$type_plr) == 'Azathioprine (I')] <- 'Azathioprine'
levels(d.aiko$type_plr)[which(levels(
  d.aiko$type_plr) == 'Mycofenolzuur (')] <- 'Mycophenolic acid'


baseline.vars <- c('age_pat', 'sexe_pat', 'tijd_bas',
                   'con_krea', 'exc_pro', 'CKD_EPI', 'tijd_dia', 'num_mmto',  
                   'bmi_poli', 'gem_tail', 'di_sig', 'con_gluc', 'gem_sbp', 'gem_dbp',
                   'con_crp', 'diabetes_fu1', 'con_insu', 'per_hba1', 'iri_homa',
                   'statin_use_fu1', 'con_tc', 'con_tg', 'con_hdlc', 'con_apoa', 'con_apob', 
                   'dag_pred', 'type_cni', 'type_plr')

baseline.table <- DescribeData(variables = baseline.vars, 
                               normal = c('con_tc', 'gem_tail', 'gem_dbp',
                                          'con_apob'), 
                               group = 'tertiles',
                               df = d.aiko)

surv.obj <- Surv(d.aiko$followup_gf_dood_2012, as.numeric(d.aiko$graft_failure_2012))

cox.1 <- coxph(surv.obj ~ c_apob48_log2, data = d.aiko)
cox.2 <- coxph(surv.obj ~ c_apob48_log2 + sexe_pat + age_pat, data = d.aiko)
cox.3 <- coxph(surv.obj ~ c_apob48_log2 + sexe_pat + age_pat + 
                 prot_urea + type_cni + type_plr, data = d.aiko)
cox.4 <- coxph(surv.obj ~ c_apob48_log2 + sexe_pat + age_pat + 
                 prot_urea + type_cni + type_plr + con_tg, 
               data = d.aiko)
cox.5 <- coxph(surv.obj ~ c_apob48_log2 + sexe_pat + age_pat + 
                 prot_urea + type_cni + type_plr + con_tg + gem_tail, 
               data = d.aiko)
# Reviewer 1 comments 6 and 8, additional adjustments in the Cox model
cox.6 <- coxph(surv.obj ~ c_apob48_log2 + sexe_pat + age_pat + 
                 prot_urea + type_cni + type_plr + con_tg + gem_tail + 
                 con_crp + crf_sv + type_sta + dag_pred, 
               data = d.aiko)
cox.7 <- coxph(surv.obj ~ c_apob48_log2 + sexe_pat + age_pat + 
                 prot_urea + type_cni + type_plr + con_tg + gem_tail + 
                 con_crp + crf_sv + type_sta + dag_pred + groente + fruit, 
               data = d.aiko)


cox.list <- list(cox.1, cox.2, cox.3, cox.4, cox.5, cox.6, cox.7)
model.names <- c('model 1', 'model 2', 'model 3', 'model 4', 'model 5', 'model 6', 'model 7')

cox.summary <- CoxSummary(models = cox.list, var = 'c_apob48_log2', summ = TRUE, styled = FALSE)
cox.summary.s <- CoxSummary(models = cox.list, var = 'c_apob48_log2', model.names = model.names, summ = TRUE, styled = TRUE)

label.models <- c('Crude analysis',
                  'Model 1 + age and sex',
                  'Model 2 + proteinuria and \nimmunosuppressive medication use',
                  'Model 3 + triglycerides',
                  'Model 4 + waist circumference',
                  'Model 5 + CRP, pre-transplant diagnosis,\nstatin use and prednisolone dose',
                  'Model 6 + daily fruit and vegetable intake')
library(grid)
#Plot HR and export to PDF file
forest.plot <- ForestPlotCox(cox.summary, 
                             label.vars = label.models,
                             cex = 1.0,
                             cex.xlab = 0.9,
                             cex.axis = 0.9,
                             box.size = 0.125,
                             graph.width = unit(3.5, 'cm'),
                             xticks = seq(0.5,3.5, by = 1.0),
                             box.col = 'black',# '#212427',
                             zero.col = 'gray',
                             lines.col = 'gray',
                             graph.pos = 3)
pdf('~/Research/Apob48_GF_AIKO/manuscript/submission/CKJ/figure_3.pdf', width = 8.5, height = 5.5)
forest.plot
dev.off()

print(cox.summary.s)

healthy.controls <- c(3.0, 3.5, 3.4, 4.5, 3.6, 5.4, 4.8, 3.6, 3.2, 3.8, 
                      5.8, 4.9, 5.1, 5.0, 4.7, 5.3, 3.2, 5.0, 5.1)

d.interactions <- InteractCox(independent.var = 'c_apob48_log2',
                              tests.interaction = list(
                                c('sexe_pat', 
                                  levels(d.aiko$sexe_pat)[1],
                                  levels(d.aiko$sexe_pat)[2]),
                                c('egfr_strat', 
                                  levels(d.aiko$egfr_strat)[1],
                                  levels(d.aiko$egfr_strat)[2]),
                                c('prot_urea', 'No', 'Yes'),
                                c('time_tx_strat', 
                                  levels(d.aiko$time_tx_strat)[1],
                                  levels(d.aiko$time_tx_strat)[2]),
                                c('waist_strat', 
                                  levels(d.aiko$waist_strat)[1],
                                  levels(d.aiko$waist_strat)[2]),
                                c('di_sig', 
                                  levels(d.aiko$di_sig)[1],
                                  levels(d.aiko$di_sig)[2]),
                                c('hba1c_strat',
                                  levels(d.aiko$hba1c_strat)[1],
                                  levels(d.aiko$hba1c_strat[2])),
                                c('tc_strat', 
                                  levels(d.aiko$tc_strat)[1], 
                                  levels(d.aiko$tc_strat)[2]),
                                c('hdl_strat', 
                                  levels(d.aiko$hdl_strat)[1], 
                                  levels(d.aiko$hdl_strat)[2]),
                                c('tgl_strat', 
                                  levels(d.aiko$tgl_strat)[1], 
                                  levels(d.aiko$tgl_strat)[2]),
                                c('apoa1_strat', 
                                  levels(d.aiko$apoa1_strat)[1], 
                                  levels(d.aiko$apoa1_strat)[2]),
                                c('apob_strat', 
                                  levels(d.aiko$apob_strat)[1], 
                                  levels(d.aiko$apob_strat)[2])
                              ),
                              surv.obj = surv.obj,
                              df = d.aiko,
                              adjust.for = '',
                              row.labels = c('Sex',
                                             'eGFR', 
                                             'Proteinuria', 
                                             'Years since TX',
                                             'Waist circumference',
                                             'Smoking',
                                             'HbA1c',
                                             'Total cholesterol',
                                             'HDL cholesterol',
                                             'Triglycerides',
                                             'Apo A-I',
                                             'Apo B'))

cox.apob.1 <- coxph(surv.obj ~ con_apob, data = d.aiko)
cox.apob.summary <- CoxSummary(models = list(cox.apob.1), var = 'con_apob', model.names = c('con_apob'))

# Revision for CKJ

quantile(d.aiko[which(d.aiko$type_cni == 'Tacrolimus'),'b_apob48_fu1'])
quantile(d.aiko[which(d.aiko$type_cni == 'Cyclosporine'),'b_apob48_fu1'])
quantile(d.aiko[which(d.aiko$type_cni == 'None'),'b_apob48_fu1'])

fit.1.cni <- coxph(surv.obj ~ c_apob48_log2 + type_cni, data = d.aiko)
fit.2.cni <- coxph(surv.obj ~ c_apob48_log2 + type_cni + c_apob48_log2:type_cni, data = d.aiko)
anova(fit.1.cni, fit.2.cni)

round((length(which(!is.na(d.aiko$groente))) / nrow(d.aiko)) * 100, 1)
round((length(which(!is.na(d.aiko$fruit))) / nrow(d.aiko)) * 100, 1)

baseline.vars <- c('fruit', 'groente')
DescribeData(variables = baseline.vars, normal = c(), group = 'tertiles', df = d.aiko)

cox.cvd.1 <- cox.1 <- coxph(Surv(d.aiko$T_CVD, as.numeric(d.aiko$di_CVD)) ~ c_apob48_log2, data = d.aiko)
summary(cox.cvd.1)

aiko.raw <- read.spss('~/Data/transplant_lines/aiko_cohort/Aiko cohort.sav', use.value.labels = TRUE, to.data.frame = TRUE)
# Use the NITRA registry to look up which KTR used which induction therapy
nitra.raw = read.spss('~/Data/nitra/nitra_basic_data.sav', use.value.labels = TRUE, to.data.frame = TRUE)
nitra = nitra.raw[which(nitra.raw$UMCGNR %in% d.aiko$UMCGNR),c('UMCGNR', 'GEBDAT', 'GESLACHT', 'KFUP_IS_IN', 'AANTTX')]
aiko.raw = merge(nitra, aiko.raw, by = 'UMCGNR')
aiko.raw = unique(aiko.raw[which(aiko.raw$AANTTX == aiko.raw$num_ntx & aiko.raw$GEBDAT == aiko.raw$dat_gebp),])
aiko.raw$KFUP_IS_IN = as.factor(as.character(aiko.raw$KFUP_IS_IN))
aiko.induction.therapy <- merge(aiko.raw[c('UMCGNR', 'dat_gebp', 'KFUP_IS_IN')], d.aiko, by = c('UMCGNR', 'dat_gebp'))
rm(aiko.raw, nitra.raw, nitra)
summary(aiko.induction.therapy$KFUP_IS_IN)
round((nrow(aiko.induction.therapy) / nrow(d.aiko)) * 100, 0)

d.aiko$c_apob48_tg_log2 <- log2((d.aiko$b_apob48_fu1/d.aiko$con_tg))
cox.suppl.1 <- coxph(surv.obj ~ c_apob48_tg_log2, data = d.aiko)
cox.suppl.2 <- coxph(surv.obj ~ c_apob48_tg_log2 + sexe_pat + age_pat, data = d.aiko)
cox.suppl.3 <- coxph(surv.obj ~ c_apob48_tg_log2 + sexe_pat + age_pat + 
                 prot_urea + type_cni + type_plr, data = d.aiko)
cox.suppl.4 <- coxph(surv.obj ~ c_apob48_tg_log2 + sexe_pat + age_pat + 
                 prot_urea + type_cni + type_plr + con_tg, 
               data = d.aiko)
cox.suppl.5 <- coxph(surv.obj ~ c_apob48_tg_log2 + sexe_pat + age_pat + 
                 prot_urea + type_cni + type_plr + con_tg + gem_tail, 
               data = d.aiko)
cox.suppl.6 <- coxph(surv.obj ~ c_apob48_tg_log2 + sexe_pat + age_pat + 
                 prot_urea + type_cni + type_plr + con_tg + gem_tail + 
                 con_crp + crf_sv + statin_use_fu1 + dag_pred, 
               data = d.aiko)
cox.suppl.7 <- coxph(surv.obj ~ c_apob48_tg_log2 + sexe_pat + age_pat + 
                 prot_urea + type_cni + type_plr + con_tg + gem_tail + 
                 con_crp + crf_sv + statin_use_fu1 + dag_pred + groente + fruit, 
               data = d.aiko)
cox.suppl.list <- list(cox.suppl.1, cox.suppl.2, cox.suppl.3, cox.suppl.4, cox.suppl.5, cox.suppl.6, cox.suppl.7)

cox.suppl.summary <- CoxSummary(models = cox.suppl.list, var = 'c_apob48_tg_log2', summ = TRUE, styled = FALSE)

CoxSummary(models = list(cox.suppl.7), var = 'c_apob48_tg_log2', summ = TRUE, styled = FALSE)

library(splines)
library(rms)

d.splines <- d.aiko[,c("c_apob48_tg_log2", "followup_gf_dood_2012", "graft_failure_2012", "sexe_pat", "age_pat", "CKD_EPI")]

d.splines <- na.omit(d.splines)
ddist <- datadist(d.splines)
refvalue <- median(d.splines$c_apob48_tg_log2)
ddist$limits$c_apob48_tg_log2[2]<-refvalue
options(datadist="ddist")
#Test how many knots fit best for the splines
for (i in 3:7) {
  fit <- rms::cph(Surv(d.splines$followup_gf_dood_2012, as.numeric(d.splines$graft_failure_2012)) ~ rcs(c_apob48_tg_log2,i) +
                    age_pat + as.numeric(sexe_pat), 
                  data = d.splines, 
                  x=TRUE)
  tmp = extractAIC(fit)
  #print(tmp)
  if(i==3) {AIC = tmp[2];nk = 3}
  if(tmp[2] < AIC) {AIC = tmp[2]; nk = i}
  #print(nk)
}

fit <- cph(Surv(d.splines$followup_gf_dood_2012, as.numeric(d.splines$graft_failure_2012)) ~ rcs(c_apob48_tg_log2,4) + age_pat + as.numeric(sexe_pat), data = d.splines, x=TRUE)
pred.hr <- Predict(fit,c_apob48_tg_log2,ref.zero=TRUE,fun=exp)

#Plot HR and export to PDF file
pdf('~/Research/Apob48_GF_AIKO/manuscript/submission/CKJ/suppl_fig_2.pdf', width = 7.5, height = 6)
PlotHr(pred.hr, d.splines, ind.var = 'c_apob48_tg_log2', 
       #breaks = seq(0.5, 5.0, by = 0.25),
       axis.1.at = seq(0.0, 5.0, by = 1.0),
       axis.1.labels = seq(0.0, 5.0, by = 1.0),
       axis.4.at = seq(0.0, 120, by = 20),
       h.wd = 1.0,
       h.border = 'darkgray',
       abline.col = 'darkgray',
       lim.x = c(0.0, 5.0),
       line.wd = 1.5,
       polygon.alpha = 0.15,
       axes.wd = 0.75,
       box = TRUE,
       cex = 0.9,
       cex.labels = 1.0,
       y.lab = 'aHR for graft failure', 
       x.lab = 'ApoB48/triglycerides ratio', 
       round.y = 0, 
       round.x = 1
)
dev.off()

# Previously another variable used, related to death censoring but this is not the right variable.
summary(glm((d.aiko$stat_pat_f8_2012 == 'overleden') ~ d.aiko$c_apob48_log2, family = 'binomial'))
round(exp(coef(glm((d.aiko$stat_pat_f8_2012 == 'overleden') ~ d.aiko$c_apob48_log2, family = 'binomial'))), 2)
round(exp(confint.default(glm((d.aiko$stat_pat_f8_2012 == 'overleden') ~ d.aiko$c_apob48_log2, family = 'binomial'))), 2)

d.aiko$vg_mi[which(d.aiko$vg_mi == '                                                            ')] <- NA
d.aiko$history_mi <- as.factor(ifelse(d.aiko$vg_mi == 'n                                                           ' 
                            | d.aiko$vg_mi == 'no                                                          ',
                            "no", "yes"))

d.aiko$vg_ticva[which(d.aiko$vg_ticva == '                        ')] <- NA
d.aiko$history_tia_cva <- as.factor(ifelse(d.aiko$vg_ticva == 'no                      ' ,
                                      "no", "yes"))

baseline.vars <- c('history_mi', 'history_tia_cva')
DescribeData(variables = baseline.vars, normal = c(), group = 'tertiles', df = d.aiko)

summary(d.aiko$di_rapa)
round((length(which(d.aiko$di_rapa == 'Sirolimus'))/nrow(d.aiko)) * 100, 1)

levels(d.aiko$type_sta)[which(levels(d.aiko$type_sta) == 'Atorvastine (Lipitor')] <- 'Atorvastin'
levels(d.aiko$type_sta)[which(levels(d.aiko$type_sta) == 'Fluvastatine (Lescol')] <- 'Fluvastatin'
levels(d.aiko$type_sta)[which(levels(d.aiko$type_sta) == 'Pravastatine (Selekt')] <- 'Pravastatin'
levels(d.aiko$type_sta)[which(levels(d.aiko$type_sta) == 'Simvastatine (Zocor)')] <- 'Simvastatin'
levels(d.aiko$type_sta)[which(levels(d.aiko$type_sta) == 'geen')] <- 'No statin use'

baseline.vars.statines <- c('age_pat', 'sexe_pat', 'b_apob48_fu1', 'con_tc', 'con_tg', 'con_hdlc', 'con_ldlc',
                            'con_apoa', 'con_apob')
baseline.table.statines <- DescribeData(variables = baseline.vars.statines, 
                               normal = c('con_tc', 'con_apob'), 
                               group = 'type_sta',
                               df = d.aiko)
if (FALSE){
  write.csv(baseline.table.statines, '~/Research/Apob48_GF_AIKO/manuscript/submission/CKJ/table_statins.csv')
}

# Summary on acute rejections
d.aiko$acute_rej = d.aiko$type_rej != 'geen'
summary(d.aiko$acute_rej)

