# Data analysis for investigating the association between normal
# triglycerides and incident hypertension in the PREVEND cohort.

data <- dataPrepPrevend('~/medical_data/PREVEND/data')

library(dplyr)
d.prevend <- filter(
  data,
  data$c_met_syndr_fu1 < 3 &
    data$b_tgl_fu1 < 1.69 &
    data$q_last_meal_fu1 != "> 1 cracker" &
    data$q_last_drink_fu1 != "> 1 cup of tea" &
    data$pe_sbp_fu1 < 130 &
    data$pe_dbp_fu1 < 80 &
    data$q_phar_antilip_fu1 != "Yes" &
    data$c_diabetes_fu1 != "Yes" &
    data$HYP_1 != "Yes" &
    is.na(data$event_hyp_fu1) == 0 &
    is.na(data$surv_yrs_hyp_fu1) == 0 &
    is.na(data$b_tgl_fu1) == 0
)

length(d.prevend$b_tgl_fu1)
summary(d.prevend$surv_yrs_hyp_fu1)
summary(d.prevend$event_hyp_fu1)

d.prevend$c_tgl_tertiles_fu1 <-
  tertile(d.prevend$b_tgl_mgdl_fu1, group = d.prevend$sex)

boxplot(
  d.prevend$b_tgl_mgdl_fu1 ~ d.prevend$event_hyp_fu1,
  outline = FALSE,
  ylab = "Triglycerides (mg/dL)",
  xlab = "",
  ylim = c(20, 155),
  las = 1,
  names = c("No hypertension", "Hypertension")
)

#Baseline characteristics table

bsl.var <-
  c(
    'b_tgl_mgdl_fu1',
    'event_hyp_fu1',
    'age_fu2',
    'sex',
    'pe_bmi_cat_fu1',
    'b_glucose_fu1',
    'pe_sbp_fu1',
    'pe_dbp_fu1',
    'b_total_chol_mgdl_fu1',
    'b_hdl_mgdl_fu1',
    'c_ldl_mgdl_fu1',
    'b_crp_fu1',
    'c_egfr_fu1',
    'q_alcohol_fu1',
    'q_curr_smoke_fu1'
  )
bsl.table <-
  DescribeData(
    variables = bsl.var,
    normal = bsl.var.nrmal.distr,
    group = 'c_tgl_tertiles_fu1',
    data = d.prevend,
    nmiss = FALSE
  )
View(bsl.table)

write.csv(bsl.table,
          "~/Desktop/baseline_table_tgl_hyp.csv",
          row.names = TRUE)


# Survival analysis

library(survival)
library(survminer)

surv_object_analysis <-
  Surv(time = d.prevend$surv_yrs_hyp_fu1,
       event = as.numeric(d.prevend$event_hyp_fu1))

d.prevend$c_tgl_log2_fu1 <- log2(d.prevend$b_tgl_mgdl_fu1)

d.prevend$pe_bmi_cat_fu1 <- droplevels(d.prevend$pe_bmi_cat_fu1)

model1 <-
  coxph(surv_object_analysis ~ c_tgl_log2_fu1, data = d.prevend)
model2 <-
  coxph(surv_object_analysis ~ c_tgl_log2_fu1 + age_fu1 + sex, data = d.prevend)
model3 <-
  coxph(surv_object_analysis ~ c_tgl_log2_fu1 + age_fu1 + sex + pe_bmi_cat_fu1,
        data = d.prevend)
model4 <-
  coxph(
    surv_object_analysis ~ c_tgl_log2_fu1 + age_fu1 + sex + pe_bmi_cat_fu1 + c_egfr_fu1,
    data = d.prevend
  )
model5 <-
  coxph(
    surv_object_analysis ~ c_tgl_log2_fu1 + age_fu1 + sex + pe_bmi_cat_fu1 + c_egfr_fu1 + c_homa_ir_fu1,
    data = d.prevend
  )
model6 <-
  coxph(
    surv_object_analysis ~ c_tgl_log2_fu1 + age_fu1 + sex + pe_bmi_cat_fu1 + c_egfr_fu1 + c_homa_ir_fu1 + q_alcohol_fu1 + q_curr_smoke_fu1,
    data = d.prevend
  )

models <- list(model1, model2, model3, model4, model5, model6)

d.cox <-
  CoxSummary(
    models = models,
    independent.var = 'c_tgl_log2_fu1',
    summ = TRUE,
    styled = FALSE,
    data = d.prevend
  )

d.cox
#write.csv(df_cox, "~/Desktop/cox_regression_table_hyp.csv", row.names = TRUE)

#Forest plot representation for the models

dark.blue <- '#284387'
  
hrzl.lines = list(
  '1' = gpar(
    lty = 1,
    lwd = 2,
    col = 'black'
  ),
  '2' = gpar(
    lty = 1,
    lwd = 2,
    col = 'black'
  ),
  '3' = gpar(lty = 1, lwd = 0.5),
  '4' = gpar(lty = 1, lwd = 0.5),
  '5' = gpar(lty = 1, lwd = 0.5),
  '6' = gpar(lty = 1, lwd = 0.5),
  '7' = gpar(lty = 1, lwd = 0.5)
)

ForestPlotCox(
  d.cox,
  label.vars = list(
    "Crude",
    "Model 1 + Age and Sex",
    "Model 2 + BMI",
    'Model 3 + eGFR',
    'Model 4 + HOMA-IR',
    'Model 5 + Smoking and alcohol use'
  ),
  hrzl.lines = hrzl.lines,
  box.size = 0.15,
  box.col = dark.blue,
  cex = 1.2,
  graph.pos = 3
)


# Continous hazard ratio

library(splines)

d.prevend.splines <-
  d.prevend[, c("b_tgl_mgdl_fu1",
                "age_fu1",
                "sex",
                "surv_yrs_hyp_fu1",
                "event_hyp_fu1")]

ddist <- datadist(d.prevend.splines)
refvalue <- median(d.prevend$b_tgl_mgdl_fu1)
ddist$limits$b_tgl_mgdl_fu1[2] <- refvalue
options(datadist = "ddist")

fit <-
  cph(
    Surv(
      d.prevend.splines$surv_yrs_hyp_fu1,
      as.numeric(d.prevend.splines$event_hyp_fu1)
    ) ~ log2(b_tgl_mgdl_fu1) + age_fu1 + sex,
    data = d.prevend.splines
  )
pred_HR <- Predict(fit, b_tgl_mgdl_fu1, ref.zero = TRUE, fun = exp)

background_col <- "#284387"
line.wd <- 3
axis.wd <- 3
  
  pdf("~/Desktop/tgl_hr_hyp.pdf")
  
  par(mar = c(5, 4, 4, 4) + 0.3)
  par(xpd = NA)
  
  ylim.bot <- min(pred_HR$lower)
  ylim.top <- max(pred_HR$upper)
  
  #Density plot
  dens <-
    density(d.prevend.splines[, colnames(d.prevend.splines) == "b_tgl_mgdl_fu1"]) #Calculate density
  plot(
    dens$x,
    dens$y,
    col = ggplot2::alpha(background_col, 0.5),
    type = "l",
    xlab = "",
    ylab = "",
    xaxt = "n",
    yaxt = "n",
    lwd = 3,
    bty = 'n'
  )
  polygon(
    dens$x,
    dens$y,
    col = ggplot2::alpha(background_col, 0.25),
    border = ggplot2::alpha(background_col, 0.5)
  )
  #axis(side = 4, at = pretty(range(dens$y))[-length(pretty(range(dens$y)))])
  #mtext("Fraction of population (Density)", side = 4, line = 4)
  
  par(new = TRUE)
  
  plot(
    pred_HR[, colnames(pred_HR) == "b_tgl_mgdl_fu1"],
    pred_HR[, colnames(pred_HR) == "yhat"],
    xlab = "Triglycerides (mg/dL)",
    ylab = paste0("aHR of hypertension"),
    type = "l",
    ylim = c(ylim.bot, ylim.top),
    col = "red",
    lwd = axis.wd,
    bty = 'n',
    font.lab = 2,
    axes = FALSE
  )
  axis(1, c(30, 60, 90, 120, 150), lwd = axis.wd, font = 2)
  axis(2,
       c(0.5, 0.7, 0.9, 1.1, 1.3, 1.5, 1.7),
       lwd = axis.wd,
       font = 2)
  lines(pred_HR[, colnames(pred_HR) == "b_tgl_mgdl_fu1"],
        pred_HR[, colnames(pred_HR) == "lower"],
        lty = 2,
        lwy = 1.5,
        lwd = axis.wd)
  lines(pred_HR[, colnames(pred_HR) == "b_tgl_mgdl_fu1"],
        pred_HR[, colnames(pred_HR) == "upper"],
        lty = 2,
        lwy = 1.5,
        lwd = axis.wd)
  lines(
    x = range(pred_HR[, colnames(pred_HR) == "b_tgl_mgdl_fu1"]),
    y = c(1, 1),
    lty = 3,
    col = "grey40",
    lwd = axis.wd
  )
  points(as.numeric(refvalue), 1, pch = 16, cex = 1.5)
  text(refvalue + 15, 0.93, paste0("ref = ", round(refvalue, 1)), font = 2)
  
  legend(
    "topleft",
    lty = c(1, 2),
    lwd = c(3, 3),
    col = c("red", "black"),
    c("estimation", "95% CI"),
    bty = "n",
    cex = 0.9,
    text.font = 2
  )
  
  dev.off()
  
# Mediation analysis
  
library(mediation)
  
#Blood pressure
d.prevend$c_sbp_log2_fu1 <- log2(d.prevend$pe_sbp_fu1)
  
summary(
  glm(
    event_hyp_fu1 ~ c_tgl_log2_fu1 + age_fu1 + sex,
    family = "binomial",
    data = d.prevend
    )
)

model_m <- lm(c_sbp_log2_fu1 ~ c_tgl_log2_fu1, data = d.prevend)
summary(model_m)

model_y <-
  glm(
    event_hyp_fu1 ~ c_sbp_log2_fu1 + c_tgl_log2_fu1 + age_fu1 + sex,
    family = "binomial",
    data = d.prevend
  )
summary(model_y)
  
mediation <- mediation::mediate(
  model.m = model_m,
  model.y = model_y,
  treat = "c_tgl_log2_fu1",
  mediator = 'c_sbp_log2_fu1',
  sims = 500
)
summary(mediation)
  
#Revision JAHA
  
print(paste0('Percent White: ', (round(
  length(which(d.prevend$race == 'Caucasian')) / nrow(d.prevend) * 100, 1
  )), '%'))
  
library(FSA)
  
DT.age <- dunnTest(as.numeric(d.prevend$age_fu1) ~ d.prevend$c_tgl_tertiles_fu1, method = 'bh')
DT.age
  
DT.glucose <- dunnTest(as.numeric(d.prevend$b_glucose_fu1) ~ d.prevend$c_tgl_tertiles_fu1, method = 'bh')
DT.glucose
  
DT.sbp <- dunnTest(as.numeric(d.prevend$pe_sbp_fu1) ~ d.prevend$c_tgl_tertiles_fu1, method = 'bh')
DT.sbp
  
DT.dbp <- dunnTest(as.numeric(d.prevend$pe_dbp_fu1) ~ d.prevend$c_tgl_tertiles_fu1, method = 'bh')
DT.dbp
  
DT.total.chol <- dunnTest(as.numeric(d.prevend$b_total_chol_fu1) ~ d.prevend$c_tgl_tertiles_fu1, method = 'bh')
DT.total.chol
  
DT.hdl <- dunnTest(as.numeric(d.prevend$b_hdl_fu1) ~ d.prevend$c_tgl_tertiles_fu1, method = 'bh')
DT.hdl
  
DT.ldl <- dunnTest(as.numeric(d.prevend$c_ldl_fu1) ~ d.prevend$c_tgl_tertiles_fu1, method = 'bh')
DT.ldl
  
DT.crp <- dunnTest(as.numeric(d.prevend$b_crp_fu1) ~ d.prevend$c_tgl_tertiles_fu1, method = 'bh')
DT.crp

DT.egfr <- dunnTest(as.numeric(d.prevend$c_egfr_fu1) ~ d.prevend$c_tgl_tertiles_fu1, method = 'bh')
DT.egfr

library(survival)
  
intr.model.1 <- coxph(surv_object_analysis ~ d.prevend$b_tgl_fu1 + d.prevend$sex)
intr.model.2 <- coxph(surv_object_analysis ~ d.prevend$b_tgl_fu1*d.prevend$sex)
summary(intr.model.2)
anova(intr.model.1, intr.model.2)

d.prevend$c_u_sodium_fu1 <- (d.prevend$u_sodium_vol1_fu1 + d.prevend$u_sodium_vol2_fu1) / 2

sodium.model.1 <- coxph(surv_object_analysis ~ log2(d.prevend$c_u_sodium_fu1))
summary(sodium.model.1)

d.prevend$c_u_uric_acid_fu1 <- (d.prevend$u_uric_acid_vol1_fu1 + d.prevend$u_uric_acid_vol2_fu1) / 2
cor(d.prevend$c_u_uric_acid_fu1, d.prevend$b_tgl_fu1, use = 'complete.obs')
  
  
