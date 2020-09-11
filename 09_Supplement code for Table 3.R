## Supplement code Table 3

# This is code to replicate Table 3 the paper:
# "Metabolomic profiling of microbial disease etiology in community-acquired pneumonia"  
# Code developed by Ilona den Hartog, Laura Zwep, and Coen van Hasselt. 

# Code contents:
# - Input data
# - Data analysis
# - Table generation
# - Output

## Setup R --------------------------------------------------------------------------------------
# R version: 3.6.3
# Rstudio version: 1.1.463

# Clear global environment
rm(list=ls())

# Set working directory and datafolder
setwd("~/ZonMW/Supplement")
datafolder <- "~/ZonMW/Supplement/"

# Load libraries
library(tidyr)
library(table1)

## Input ---------------------------------------------------------------------------------------

# Load elastic net data
load("~/ZonMW/Supplement/Results/03_Elastic_net_results.Rdata")

## Data analysis -------------------------------------------------------------------------------

### Mean sensitivity and specificity
#Calculate mean sensitivity and specificity for all folds and repeats
outerCV_sens_spec <- data.frame(matrix(data=NA, nrow = 15*100, ncol = 6))
names(outerCV_sens_spec) <- c("repeatnr", "comparison", "fold", "Sensitivity", "Specificity", "BER")
for (r in 1:100){
  start <- 1 + (r-1)*15
  stop <- r*15
  
  df <- repeat_data[[r]]$elnet_data      #ca = TP, cb =FP, cc = FN, cd = TN
  df$Sensitivity <- df$ca/(df$ca+df$cc)
  df$Specificity <- df$cd/(df$cb+df$cd)
  df$BER <- 0.5*( df$cb / (df$ca + df$cb) + df$cc / (df$cc + df$cd))
  df2 <- subset(df, select = c("comparison", "fold", "Sensitivity", "Specificity", "BER"))
  
  outerCV_sens_spec[start:stop,] <- cbind("repeatnr" = rep(r, 15), df2)
}

# Sensitivity
## Mean sensitivity for all folds and repeats (outerCV)
atyp_pneu_SensMean <- mean(outerCV_sens_spec$Sensitivity[outerCV_sens_spec$comparison==1], na.rm = TRUE)
pneu_vir_SensMean <- mean(outerCV_sens_spec$Sensitivity[outerCV_sens_spec$comparison==2], na.rm = TRUE)
atyp_vir_SensMean <- mean(outerCV_sens_spec$Sensitivity[outerCV_sens_spec$comparison==3], na.rm = TRUE)

## Standard deviation sensitivity for all folds and repeats (outerCV)
atyp_pneu_SensSD <- sd(outerCV_sens_spec$Sensitivity[outerCV_sens_spec$comparison==1], na.rm = TRUE)
pneu_vir_SensSD <- sd(outerCV_sens_spec$Sensitivity[outerCV_sens_spec$comparison==2], na.rm = TRUE)
atyp_vir_SensSD <- sd(outerCV_sens_spec$Sensitivity[outerCV_sens_spec$comparison==3], na.rm = TRUE)

#Specitivity
## Mean specificity for all folds and repeats (outerCV)
atyp_pneu_SpecMean <- mean(outerCV_sens_spec$Specificity[outerCV_sens_spec$comparison==1], na.rm = TRUE)
pneu_vir_SpecMean <- mean(outerCV_sens_spec$Specificity[outerCV_sens_spec$comparison==2], na.rm = TRUE)
atyp_vir_SpecMean <- mean(outerCV_sens_spec$Specificity[outerCV_sens_spec$comparison==3], na.rm = TRUE)

## Standard deviation BER for all folds and repeats (outerCV)
atyp_pneu_SpecSD <- sd(outerCV_sens_spec$Specificity[outerCV_sens_spec$comparison==1], na.rm = TRUE)
pneu_vir_SpecSD <- sd(outerCV_sens_spec$Specificity[outerCV_sens_spec$comparison==2], na.rm = TRUE)
atyp_vir_SpecSD <- sd(outerCV_sens_spec$Specificity[outerCV_sens_spec$comparison==3], na.rm = TRUE)


## Table generation --------------------------------------------------------------------------------
factsheet_sens_spec <- data.frame(comparison = c("S. pneumoniae - (atypical + viral)", "Atypical bacterial - (S.pneumoniae + viral)", "Viral - (S. pneumoniae + atypical)"),
                                  Sensitivity = c(paste(signif(atyp_pneu_SensMean, digits = 2), " (", signif(atyp_pneu_SensSD, digits = 2), ")", sep=""),
                                                  paste(signif(pneu_vir_SensMean, digits = 2), " (", signif(pneu_vir_SensSD, digits = 2), ")", sep=""),
                                                  paste(signif(atyp_vir_SensMean, digits = 2), " (", signif(atyp_vir_SensSD, digits = 2), ")", sep="")),
                                  Specitivity = c(paste(signif(atyp_pneu_SpecMean, digits = 2), " (", signif(atyp_pneu_SpecSD, digits = 2), ")", sep=""),
                                                  paste(signif(pneu_vir_SpecMean, digits = 2), " (", signif(pneu_vir_SpecSD, digits = 2), ")", sep=""),
                                                  paste(signif(atyp_vir_SpecMean, digits = 2), " (", signif(atyp_vir_SpecSD, digits = 2), ")", sep="")))
#View(factsheet_sens_spec)

## Output --------------------------------------------------------------------------------------------------------

# Create tables folder
dir.create(paste(datafolder,"Tables", sep=""))

# Save table as .csv file
# Write csv file
write.csv2(factsheet_sens_spec, file = paste(datafolder, "Tables/", "Table3.csv", sep=""), row.names = FALSE)

