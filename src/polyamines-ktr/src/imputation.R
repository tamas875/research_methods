
library(missForest)
set.seed(123)
names(elspa)
elspa_fu <- elspa
summary(elspa_fu)

elspa_imp <- missForest(elspa_fu, ntree = 100, verbose = TRUE, maxiter = 100)
elspa_fu <- elspa_imp$ximp
elspa_fu <- as.data.frame(elspa_fu)
summary(elspa_fu)
summary(elspa)

#### Finish ####
