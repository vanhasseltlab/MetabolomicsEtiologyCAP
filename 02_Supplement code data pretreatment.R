## Supplement code data pretreatment

# This is code to replicate the data pretreatment steps from the paper:
# "Metabolomic profiling of microbial disease etiology in community-acquired pneumonia"  
# Code developed by Ilona den Hartog, Laura Zwep, and Coen van Hasselt. 

# Code contents:
# - Data preatreament: Applying log-transformation and 
#   applying standardization to correct for heteroscedasticity

## Setup R --------------------------------------------------------------------------------------
# R version: 3.6.3
# Rstudio version: 1.1.463

# Clear global environment
rm(list=ls())

# Set working directory and datafolder
setwd("~/ZonMW/Supplement")
datafolder <- "~/ZonMW/Supplement/"

# Load libraries

## Input ---------------------------------------------------------------------------------------
# Load metabolomics data after quality control
diagnosisdata_clean <- read.csv2(paste0(datafolder, "Data/", "01_diagnosisdata_clean.csv"))

# Define metabolite and metadata range in data after quality control
# Columns before metabolite data are defined as 'begin'
begin = 1
# Columns after metabolite data are defined as 'end'
end = (length(diagnosisdata_clean)-3):(length(diagnosisdata_clean))
# Columns containing the metabolite data are defined as 'metrange'
metrange <- 2:(length(diagnosisdata_clean)-4)

## Data pretreatment ---------------------------------------------------------------------------
# - Metabolite data follows a log-normal distribution as is expected according to biology. 
#   Therfore, we perform a log-transformation on all metabolites.
# - To correct for heteroscedasticity, we standardize per metabolite by 
#   subtracting the metabolite mean and dividing by the metabolite standard deviation 
#   to reach a mean of 0 and a standard deviation of 1. 

# Log2 transformation of metabolite values +1. 
#  Log2 transformation is performed on the metabolite values + 1 to remove zero values, 
#  because 0 values become infinite after log transformation and should therefore be avoided.
logdata <- diagnosisdata_clean
logdata[,metrange] <- apply(logdata[,metrange], MARGIN = 2, FUN = function(x) { 
  log2(x+1)
})

# Autoscaling of log transformed data. 
diagnosisdata_pretreated <- logdata
diagnosisdata_pretreated[,metrange] <- apply(diagnosisdata_pretreated[,metrange], MARGIN=2, FUN=function(x) {
  (x - mean(x)) / sd(x)
})

## Output --------------------------------------------------------------------------------------
# Write csv file
write.csv2(diagnosisdata_pretreated, file = paste(datafolder, "Data/", "02_diagnosisdata_pretreated.csv", sep=""), row.names = FALSE)
# Save Rdata file
save(diagnosisdata_pretreated, file = paste(datafolder, "Data/", "02_diagnosisdata_pretreated.Rdata", sep=""))


