library(dplyr)
library(psych)

#### Merge the data ####
else_data <- rename(else_data, UMCG_nummer = "UMCG.")
pa_sub <- pa[255:877, c(1, 15:21, 23)]
donor_sub <- pa[5:251, c(1, 15:21, 23)]

colnames(pa_sub) <- as.character(pa_sub[1, ])
print(colnames(pa_sub)[9])
print(colnames(pa_sub)[1])
colnames(pa_sub)[9] <- "n1_n8"
colnames(pa_sub)[1] <- "Subjectnr"
pa_sub <- pa_sub[-1, ]
pa_sub <- pa_sub %>%
  mutate(across(everything(), ~as.numeric(as.character(.))))
pa_sub <- na.omit(pa_sub)

colnames(donor_sub) <- as.character(donor_sub[1, ])
print(colnames(donor_sub)[9])
print(colnames(donor_sub)[1])
colnames(donor_sub)[9] <- "n1_n8"
colnames(donor_sub)[1] <- "Subjectnr"
colnames(donor_sub)[2] <- "AcPUT"
colnames(donor_sub)[3] <- "N8AcSPD"
colnames(donor_sub)[4] <- "N1AcSPD"


donor_sub <- donor_sub[-1, ]
donor_sub <- donor_sub %>%
  mutate(across(everything(), ~as.numeric(as.character(.))))
donor_sub <- na.omit(donor_sub)



elspa <- else_data %>%
  inner_join(outcm, by = "UMCG_nummer") %>%
  inner_join(pa_sub, by = "Subjectnr") %>%
  select(
    Subjectnr,
    Age_visit_1,
    Geslacht,
    BMI1,
    SerumKreat_1,
    CKD_EPI_1,
    U_Tot_Eiwit_24h_totaal_1,
    Dummy_Proteinurie,
    U_Kalium_24h_totaal_1,
    PRE_EMPTIVE,
    follow_up1,
    SBP1,
    Roken_1,
    Dummy_Diabetes,
    Statine_1,
    Dummy_Antihypertensives,
    Hb_1,
    HbA1c1,
    hsCRP_1,
    Serum_Albumine_1,
    Serum_cholesterol_1,
    Serum_HDL_Cholesterol_1,
    Prednisolon_1,
    Calcineurineremmer_1,
    Proliferatieremmer_1,
    gdag_brood,
    gdag_peulvruchten,
    gdag_groente,
    gdag_fruit,
    gdag_vis,
    gdag_melkprod,
    gdag_koffie,
    gdag_thee,
    gdag_vleesprod,
    SumOfkCal,
    SumOfvezel,
    SumOfeiwittot,
    SumOfvettot,
    SumOfkhtot,
    SumOfalcohol,
    Vol_24hU_1,
    AcPUT,
    N8AcSPD,
    N1AcSPD,
    PUT,
    CAD,
    SPD,
    SPM,
    time_gf_fu4,
    event_gf_fu4,
    time_mort2_fu4,
    event_mort2_fu4,
    event_Miscell2_fu4,
    event_cvd2_fu4_v2,
    event_infect2_fu4,
    event_Malignant2_fu4,
    time_mort_fu4,
    event_mort_fu4
  )
summary(elspa)

donor <- donor_data %>%
  inner_join(donor_sub, by = "Subjectnr") %>%
  select(
    Subjectnr,
    Age_fu1,
    Gender,
    BMI_fu1,
    eGFR_fu1,
    U_24h_Vol_fu1,
    AcPUT,
    N8AcSPD,
    N1AcSPD,
    PUT,
    CAD,
    SPD,
    SPM
  )

#### Calculate PA excretion ####

excretion_variables <- c("AcPUT", "N8AcSPD", "N1AcSPD",
                         "PUT", "CAD", "SPD", "SPM")

elspa <- elspa %>%
  mutate(across(all_of(excretion_variables),
                list(h24 = ~.x * Vol_24hU_1 / 1000),
                .names = "{.col}_24h")) %>%
  mutate(totalPA_24h = rowSums(across(ends_with("_24h")), na.rm = TRUE))

donor <- donor %>%
  mutate(across(all_of(excretion_variables),
                list(h24 = ~.x * U_24h_Vol_fu1 / 1000),
                .names = "{.col}_24h")) %>%
  mutate(totalPA_24h = rowSums(across(ends_with("_24h")), na.rm = TRUE))


#### Convert Prednisolon_1 to numeric (prednisolone dose) ####
elspa$Prednisolon_1 <- as.numeric(elspa$Prednisolon_1)

#### Active smoking ####
elspa$Roken_1 <- factor(ifelse(elspa$Roken_1 %in% c("Nooit", "Ooit"), 0,
                               ifelse(elspa$Roken_1 == "Ja", 1, NA)),
                        levels = c(0, 1), labels = c("No", "Yes"))
table(elspa$Roken_1)


#### Statine use ####
elspa$Statine_1 <- factor(ifelse(elspa$Statine_1 %in% "Geen", 0, 1),
                          levels = c(0, 1), labels = c("No", "Yes"))

table(elspa$Statine_1)

#### Finish ####

