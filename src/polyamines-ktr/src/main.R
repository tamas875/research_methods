library(rstudioapi)
library(foreign)
library(haven)
library(readxl)


# Get script directory 
if (!(exists('current_wd'))) {
  current_wd = paste0(sub('/[^/]*$', '', getSourceEditorContext()$path), '/')
}

data_dirs <- read.csv(file.path(current_wd, "data_dirs.csv"))

# Import all data
else_data <- read.spss(data_dirs[1, ], to.data.frame = TRUE)
outcm <- read.spss(data_dirs[2, ], to.data.frame = TRUE)
pa <- read_excel(data_dirs[3, ])
donor_data <- read.spss(data_dirs[4, ], to.data.frame = TRUE)

# Preprocess the data
source(file.path(current_wd, "preprocess_data.R"))

# Describe the data
source(file.path(current_wd, "describe_data.R"))

# Describe the data2
source(file.path(current_wd, "describe_data2.R"))

# Imputation
source(file.path(current_wd, "imputation.R"))

# Run linear regression models
source(file.path(current_wd, "linear_regression.R"))

# Run cox regression models
source(file.path(current_wd, "cox_regression.R"))

# Run forest plots
source(file.path(current_wd, "forest_plot.R"))

# Run spline models
source(file.path(current_wd, "splines.R"))

# Run interaction analyses
source(file.path(current_wd, "interaction.R"))

# Run sensitivity analyses
source(file.path(current_wd, "cox_median.R"))
