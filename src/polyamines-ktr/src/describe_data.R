library(nortest)
library(CBCgrps)


# Group totalPA based on the terciles
elspa_des <- elspa
quantiles_else <- quantile(elspa_des$totalPA_24h,
                           probs = c(0, 1 / 3, 2 / 3, 1),
                           na.rm = TRUE)



elspa_des$totalPAgrp <- cut(elspa_des$totalPA_24h,
                            breaks = quantiles_else,
                            include.lowest = TRUE,
                            labels = c("T1", "T2", "T3"))


#### Baseline characters ####
summary(elspa_des)
names(elspa_des)

elspa_des1 <- elspa_des %>%
  select(totalPAgrp, totalPA_24h,
         Age_visit_1 : SumOfalcohol)
elspa_des1 <- as.data.frame(elspa_des1)


# Distribution - Q-Q plots
lapply(names(elspa_des1)[sapply(elspa_des1, is.numeric)], function(var) {
  # Q-Q plot
  qqnorm(elspa_des1[[var]],
         main = paste("Q-Q Plot of", var),
         col = "darkblue")
  qqline(elspa_des1[[var]],
         col = "red",
         lwd = 2)
})

table1 <- multigrps(df = elspa_des1,
                    "totalPAgrp",
                    norm.rd = 1,
                    sk.rd = 1,
                    cat.rd = 0)
print(table1)

write.csv(table1,"PA_BaselineCharacter.csv")

#### Finish ####