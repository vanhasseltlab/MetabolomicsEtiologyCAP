## Supplement code data cleaning

# This is code to replicate the data cleaning steps from the paper:
# "Metabolomic profiling of microbial disease etiology in community-acquired pneumonia"  
# Code developed by Ilona den Hartog, Laura Zwep, and Coen van Hasselt. 

# Code contents:
# - Data cleaning: Removing missing data after metabolomics quality control

## Setup R --------------------------------------------------------------------------------------
# R version: 3.6.3
# Rstudio version: 1.1.463

# Clear global environment
rm(list=ls())

# Set working directory and datafolder
setwd("~/ZonMW/Supplement")
datafolder <- "~/ZonMW/Supplement/"

# Load libraries
library(readxl)


## Input data ----------------------------------------------------------------------------------
# Load metabolomics data after quality control
diagnosisdata_afterQC <- read_excel(paste0(datafolder, "Data/", "Supplementary material III, Table S4, metabolomics data after quality control.xlsx"), na="NA")

# Define metabolite and metadata range in data after quality control
# Columns before metabolite data are defined as 'begin'
begin = 1
# Columns after metabolite data are defined as 'end'
end = (length(diagnosisdata_afterQC)-3):(length(diagnosisdata_afterQC))
# Columns containing the metabolite data are defined as 'metrange'
metrange <- 2:(length(diagnosisdata_afterQC)-4)

# Define metabolite data as numeric data for analysis
diagnosisdata_afterQC[,metrange] <- apply(diagnosisdata_afterQC[,metrange], 2, function(x) {
  as.numeric(as.character(x))
})

## Data cleaning --------------------------------------------------------------------------------
# - Samples (in rows) with more then 10 missing metabolites are removed. 
#   This removes samples with no measurements in at least one platform.
# - Metabolites (in columns) with missing values are removed.

# If row has multiple NA's (>10), remove row from dataframe
# Subset metabolites from metabolomics data
df_sub <- diagnosisdata_afterQC[-c(begin, end)]
# Calculate the numer of NA's per row
narows <- apply(df_sub, 1, function(X){sum(is.na(X))})
# Subset metabolomics data keeping only rows with les then 10 missing metabolites
df2 <- diagnosisdata_afterQC[narows < 10,]

# If column has one NA, remove column. This removes metabolites with missing values
# Subset metabolites metabolomics data
df2_sub <- df2[-c(begin, end)]
# Calculate the numer of NA's per column
nacols <- apply(df2_sub, 2, function(X){sum(is.na(X))})
# Subset metabolomics data keeping only columns without missing metabolites
keep <- colnames(df2_sub[, nacols==0])
# Combine metadata with clean metabolite data again
diagnosisdata_clean <- cbind(sampleID = df2[,begin], df2[,keep], df2[,end])

## Output data ----------------------------------------------------------------------------------

# Save clean metabolomics dataframe in .csv and Rdata files
write.csv2(diagnosisdata_clean, file = paste(datafolder, "Data/", "01_diagnosisdata_clean.csv", sep=""), row.names = FALSE)
save(diagnosisdata_clean, file = paste(datafolder, "Data/", "01_diagnosisdata_clean.Rdata", sep=""))


