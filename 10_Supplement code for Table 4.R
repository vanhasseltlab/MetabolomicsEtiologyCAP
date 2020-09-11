## Supplement code Table 4

# This is code to replicate Table 4 the paper:
# "Metabolomic profiling of microbial disease etiology in community-acquired pneumonia"  
# Code developed by Ilona den Hartog, Laura Zwep, and Coen van Hasselt. 

# Code contents:
# - Input data
# - Data analysis:
#   - Elastic net model without age and sex
#   - Elastic net model with age and sex
#   - Linear model of age and sex
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
elnetdata_without <- repeat_data

# Load elastic net data with covariates age and sex
load("~/ZonMW/Supplement/Results/04_Elastic_net_with_covariates_results.Rdata")
elnetdata_with <- repeat_data

# Load linear regression of covariates age and sex 
load("~/ZonMW/Supplement/Results/05_Linear_regression_of_covariates_results.Rdata")
lineardata <- repeat_data

resultsdata <- list(elnetdata_without, elnetdata_with, lineardata)

## Data analysis: Elastic net model without age and sex ---------------------------------------------

factsheet <- data.frame(Comparison = c("S. pneumoniae - (atypical + viral)", 
                                       "Atypical bacterial - (S. pneumoniae + viral)", 
                                       "Viral - (S. pneumoniae + atypical)"),
                        `Elastic net without age and sex: BER` = NA, 
                        `Elastic net without age and sex: nvars` = NA, 
                        `Elastic net with age and sex: BER` = NA, 
                        `Elastic net with age and sex: nvars` = NA, 
                        `Linear model age and sex`= NA)
k = 2
l = 3
for (a in resultsdata) { 
  #Put all optimal parameters (a, l and BER) from innerCV for all folds and repeats in one dataframe
  innerCV_oppars_BER <- data.frame(matrix(data=NA, nrow = 15*100, ncol = 6))
  names(innerCV_oppars_BER) <- c("repeatnr", "comparison", "fold", "alpha", "lambda", "BER")
  for (r in 1:100){
    start <- 1 + (r-1)*15
    stop <- r*15
    innerCV_oppars_BER[start:stop,] <- cbind("repeatnr" = rep(r, 15), a[[r]]$optimal_parameters)
  }
  
  # Calculate BER for all folds and repeats (per alpha and lambda combination resulting from innerCV)
  outerCV_BER_a_l <- data.frame(matrix(data=NA, nrow = 15*100, ncol = 4))
  names(outerCV_BER_a_l) <- c("repeatnr", "comparison", "fold", "BER")
  for (r in 1:100){
    start <- 1 + (r-1)*15
    stop <- r*15
    
    df <- a[[r]]$elnet_data
    df$BER <- 0.5*( df$cb / (df$ca + df$cb) + df$cc / (df$cc + df$cd))
    df2 <- subset(df, select = c("comparison", "fold", "BER"))
    
    outerCV_BER_a_l[start:stop,] <- cbind("repeatnr" = rep(r, 15), df2)
  }
  outerCV_BER_a_l <- outerCV_BER_a_l %>% 
    mutate(alpha = innerCV_oppars_BER$alpha) %>% 
    mutate(lambda = innerCV_oppars_BER$lambda)
  
  
  # Make dataframe that includes the amount of variables selected:
  outerCV_nvars <- data.frame(matrix(data = NA, nrow = 5*3*100, ncol = 4))
  names(outerCV_nvars) <- c("repeatnr", "comparison", "fold", "nvars")
  for (r in 1:100){
    for (i in 1:3){
      for (f in 1:5){
        start <- 1 + (f-1) + (i-1)*5 + (r-1)*15
        stop <- f + (i-1)*5 + (r-1)*15
        
        nvars <- as.numeric(as.character(sum(a[[r]]$elnet_model[[i]][[f]]$beta[,1] != 0)))
        outerCV_nvars[start:stop,] <- cbind("repeatnr" = r, "comparison" = i, "fold" = f, "nvars" = nvars)
      }
    }
  }
  
  # Make dataframe that includes a, l, BER and nvars outerCV
  outerCV_BER_nvars <- data.frame(cbind(outerCV_BER_a_l, nvars = outerCV_nvars$nvars))
  
  ### Numbers for paper Elasic Net without covariates:
  ## Mean BER for all folds and repeats (outerCV)
  pneu_BERmean <- mean(outerCV_BER_a_l$BER[outerCV_BER_a_l$comparison==1])
  atyp_BERmean <- mean(outerCV_BER_a_l$BER[outerCV_BER_a_l$comparison==2])
  vir_BERmean <- mean(outerCV_BER_a_l$BER[outerCV_BER_a_l$comparison==3])
  
  ## Standard deviation BER for all folds and repeats (outerCV)
  pneu_BERsd <- sd(outerCV_BER_a_l$BER[outerCV_BER_a_l$comparison==1])
  atyp_BERsd <- sd(outerCV_BER_a_l$BER[outerCV_BER_a_l$comparison==2])
  vir_BERsd <- sd(outerCV_BER_a_l$BER[outerCV_BER_a_l$comparison==3])
  
  # Add mean (sd) BER to table
  factsheet[1,k] <- paste(signif(pneu_BERmean, digits = 2), " (", signif(pneu_BERsd, digits = 2), ")", sep="")
  factsheet[2,k] <- paste(signif(atyp_BERmean, digits = 2), " (", signif(atyp_BERsd, digits = 2), ")", sep="")
  factsheet[3,k] <- paste(signif(vir_BERmean, digits = 2), " (", signif(vir_BERsd, digits = 2), ")", sep="")
  
  
  ## Mean number of variables selected for all folds and repeats (outerCV)
  pneu_nvarsmean <- mean(outerCV_nvars$nvars[outerCV_nvars$comparison==1])
  atyp_nvarsmean <- mean(outerCV_nvars$nvars[outerCV_nvars$comparison==2])
  vir_nvarsmean <- mean(outerCV_nvars$nvars[outerCV_nvars$comparison==3])
  
  ## Standard deviation number of variables selected for all folds and repeats (outerCV)
  pneu_nvarssd <- sd(outerCV_nvars$nvars[outerCV_nvars$comparison==1])
  atyp_nvarssd <- sd(outerCV_nvars$nvars[outerCV_nvars$comparison==2])
  vir_nvarssd <- sd(outerCV_nvars$nvars[outerCV_nvars$comparison==3])

  # Add number of variables to table
  factsheet[1,l] <- paste(signif(pneu_nvarsmean, digits = 2), " (", signif(pneu_nvarssd, digits = 2), ")", sep="")
  factsheet[2,l] <- paste(signif(atyp_nvarsmean, digits = 2), " (", signif(atyp_nvarssd, digits = 2), ")", sep="")
  factsheet[3,l] <- paste(signif(vir_nvarsmean, digits = 2), " (", signif(vir_nvarssd, digits = 2), ")", sep="")
  
  k = k + 2
  l = l + 2
}


## Output ------------------------------------------------------------------------------------------
# Save table as .csv file
# Write csv file
write.csv2(factsheet, file = paste(datafolder, "Tables/", "Table4.csv", sep=""), row.names = FALSE)

