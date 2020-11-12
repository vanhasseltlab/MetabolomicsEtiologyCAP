## This is code to replicate the data analysis steps from the paper:
## "Metabolomic profiling of microbial disease etiology in community-acquired pneumonia"  
## Code developed by Ilona den Hartog, Laura Zwep, and Coen van Hasselt. 

# Code contents
# - Data cleaning
# - Data pretreatment
# - Data analysis for the comparisons:
#   - Comparison 1: S. pneumoniae versus atypical bacteria + virusses
#   - Comparison 2: Atypical bacteria versus S. pneumoniae + virusses
#   - Comparison 3: Virusses versus S. pneumoniae + atypical bacteria
# - Statistical data analysis:
#   - T-test + FDR multiple testing correction for comparison of 1 pathogen 
#      group versus the remaining 2 pathogen groups. 
  

## Setup R --------------------------------------------------------------------
# R version: 3.6.3
# Rstudio version: 1.1.463

# Load libraries
library(tidyverse)
library(data.table)
library(dplyr)
library(glmnet)
library(caret)
library(ROCR)
library(ggplot2)
library(ggpubr)
library(readxl)

# Source functions
source("functions/LeidenColoring.R")
source("functions/DataCleaning.R")
source("functions/Ttesting.R")
source("functions/MakeFolds.R")
source("functions/LogisticRegressionFolds.R")
source("functions/ElasticNetCV.R")
source("functions/CovElasticNetCV.R")
source("functions/ElasticNetPerformance.R")
source("functions/CovElasticNetPerformance.R")
source("functions/ElnetROC.R")


## Input data -----------------------------------------------------------------
# Load metabolomics data after quality control
data.reduced <- read.csv2("data/00_data_reduced.csv")

# Define metabolite and metadata range in data after quality control
nonmetrange <- c("sample.id", "pathogen.group", 
                 "age", "sex", "psi.score")
metrange <- setdiff(names(data.reduced), nonmetrange)
cov <- c("age", "sex")


## Data cleaning --------------------------------------------------------------
data.clean <- DataCleaning(data.reduced, metrange, nonmetrange)
#Reset metrange after cleaning
metrange <- setdiff(names(data.clean), nonmetrange)

## Data pretreatment ----------------------------------------------------------
# Log2 transformation of metabolite values +1. 
logdata <- data.clean
logdata[, metrange] <- apply(logdata[, metrange], MARGIN = 2, FUN = function(x) { 
  log2(x + 1)
})

# Autoscaling of log transformed data. 
data.pretreated <- logdata
data.pretreated[, metrange] <- apply(data.pretreated[, metrange], MARGIN = 2, FUN = function(x) {
  (x - mean(x)) / sd(x)
})

## Data per comparison -------------------------------------------------------
dat <- data.pretreated

## ATYPICAL versus all other groups
# Add column specificating if pathogen is atypical or other
dat <- mutate(dat, "atyp.other" = NA)
for (i in 1:nrow(dat)){
  if (dat$pathogen.group[i] == "Viral"){
    dat$atyp.other[i] <- "other"
  } else if (dat$pathogen.group[i] == "S. pneumoniae") { 
    dat$atyp.other[i] <- "other"
  } else {
    dat$atyp.other[i] <- "atyp"
  }
}
dat$atyp.other <- as.factor(dat$atyp.other)

## S.PNEUMONIAE verus all other groups
# Add column specificating if pathogen is s.pneumoniae or other
dat <- mutate(dat, "spneu.other" = NA)
for (i in 1:nrow(dat)){
  if (dat$pathogen.group[i] == "Viral"){
    dat$spneu.other[i] <- "other"
  } else if (dat$pathogen.group[i] == "Atypical") { 
    dat$spneu.other[i] <- "other"
  } else {
    dat$spneu.other[i] <- "spneu"
  }
}
dat$spneu.other <- as.factor(dat$spneu.other)

## VIRAL pathogens versus all other groups
# Add column specificating if pathogen is bacterial or viral
dat <- mutate(dat, "bac.vir" = NA)
for (i in 1:nrow(dat)){
  if (dat$pathogen.group[i] == "Viral"){
    dat$bac.vir[i] <- "Viral"
  } else {
    dat$bac.vir[i] <- "Bacterial"
  }
}
dat$bac.vir <- as.factor(dat$bac.vir)

## Data analysis univariate: Students T-test----------------------------------
# Atypical versus other groups
pval.atyp.other <- Ttesting(metrange, dat$atyp.other, "atyp", "other", dat)
sign.mets.atyp.other <- pval.atyp.other[pval.atyp.other$q.value < 0.05, ]

# S. pneumoniae versus other groups
pval.spneu.other <- Ttesting(metrange, dat$spneu.other, "spneu", "other", dat)
sign.mets.spneu.other <- pval.spneu.other[pval.spneu.other$q.value < 0.05, ]

# Viral versus other groups (bacterial)
pval.bac.vir <- Ttesting(metrange, dat$bac.vir, "Viral", "Bacterial", dat)
sign.mets.bac.vir <- pval.bac.vir[pval.bac.vir$q.value < 0.05, ]

## Performance of single metabolites in cross validated logistic regression models ----
# Logistic regression model of significant metabolites for atyp versus other groups
logit.glycylglycine <- LogisticRegressionFolds(dat, variable = "Glycylglycine", 
                                               outcome = "atyp.other", 
                                               repeats = 100, folds = 5)
logit.sdma <- LogisticRegressionFolds(dat, variable = "SDMA", 
                                      outcome = "atyp.other", 
                                      repeats = 100, folds = 5)
logit.lpi18.1 <- LogisticRegressionFolds(dat, variable = "LPI.18.1.", 
                                         outcome = "atyp.other", 
                                         repeats = 100, folds = 5)
logit.age.sex <- LogisticRegressionFolds(dat, variable = c("age", "sex"), 
                                         outcome = "atyp.other", 
                                         repeats = 100, folds = 5)
# Linear model all significant metabolites combined
logit.sigatyp <- LogisticRegressionFolds(dat, variable = c("Glycylglycine", "SDMA", 
                                                           "LPI.18.1."),
                                    outcome = "atyp.other", 
                                    repeats = 100, folds = 5)
# Linear model all significant metabolites + age
logit.sigatyp.age.sex <- LogisticRegressionFolds(dat, variable = c("Glycylglycine", 
                                                                   "SDMA", 
                                                                   "LPI.18.1.",
                                                                   "age", "sex"), 
                                                 outcome = "atyp.other", 
                                                 repeats = 100, folds = 5)


## Data analysis multivariate: Elastic net -----------------------------------

## Run elastic net analysis WITHOUT covariates: 
 # *** approximately 5H runtime per comparison ***
# # Atypical versus other
# elnet.atyp.other.results <- ElasticNetCV(dat, metrange, outcome = "atyp.other", 
#                                  alphas = seq(0, 1, 0.05) , lambdas = seq(0, 1, 0.05), 
#                                  repeats = 100, folds = 5)
# save(elnet.atyp.other.results, file = "results/elnet.atyp.other.results.R")
# 
# # S.pneumoniae versus other
# elnet.spneu.other.results <- ElasticNetCV(dat, metrange, outcome = "spneu.other", 
#                                          alphas = seq(0, 1, 0.05) , lambdas = seq(0, 1, 0.05), 
#                                          repeats = 100, folds = 5)
# save(elnet.spneu.other.results, file = "results/elnet.spneu.other.results.R")
# 
# # Viral versus bacterial
# elnet.bac.vir.results <- ElasticNetCV(dat, metrange, outcome = "bac.vir", 
#                                          alphas = seq(0, 1, 0.05) , lambdas = seq(0, 1, 0.05), 
#                                          repeats = 100, folds = 5)
# save(elnet.bac.vir.results, file = "results/elnet.bac.vir.results.R")
## OR LOAD RESULTS
load("results/elnet.atyp.other.results.R")
load("results/elnet.spneu.other.results.R")
load("results/elnet.bac.vir.results.R")

## Run elastic net analysis WITH covariates sex and age:
# *** approximately 5H runtime per comparison ***
dat$sex <- as.numeric(dat$sex)
# # Atypical versus other
# cov.elnet.atyp.other.results <- CovElasticNetCV(dat, metrange, cov, 
#                                                 outcome = "atyp.other", 
#                                                 alphas = seq(0, 1, 0.05), 
#                                                 lambdas = seq(0, 1, 0.05),
#                                                 repeats = 100, folds = 5)
# save(cov.elnet.atyp.other.results, file = "results/cov.elnet.atyp.other.results.R")
# # Spneumoniae versus other
# cov.elnet.spneu.other.results <- CovElasticNetCV(dat, metrange, cov,  
#                                                 outcome = "spneu.other", 
#                                                 alphas = seq(0, 1, 0.05), 
#                                                 lambdas = seq(0, 1, 0.05),
#                                                 repeats = 100, folds = 5)
# save(cov.elnet.spneu.other.results, file = "results/cov.elnet.spneu.other.results.R")
# # Viral versus bacterial
# cov.elnet.bac.vir.results <- CovElasticNetCV(dat, metrange, cov, 
#                                                 outcome = "bac.vir", 
#                                                 alphas = seq(0, 1, 0.05) , 
#                                                 lambdas = seq(0, 1, 0.05),
#                                                 repeats = 100, folds = 5)
# save(cov.elnet.bac.vir.results, file = "results/cov.elnet.bac.vir.results.R")
## OR LOAD RESULTS
load("results/cov.elnet.atyp.other.results.R")
load("results/cov.elnet.spneu.other.results.R")
load("results/cov.elnet.bac.vir.results.R")


## Results elastic net model performance within cross validation --------------
#Without covariates
elnet.atyp.other.performance <- ElasticNetPerformance(dat, metrange, outcome = "atyp.other", 
                                                       repeats = 100, folds = 5, 
                                                       elnet.atyp.other.results)
elnet.spneu.other.performance <- ElasticNetPerformance(dat, metrange, outcome = "spneu.other", 
                                                      repeats = 100, folds = 5, 
                                                      elnet.spneu.other.results)
elnet.bac.vir.performance <- ElasticNetPerformance(dat, metrange, outcome = "bac.vir", 
                                                   repeats = 100, folds = 5, 
                                                   elnet.bac.vir.results)
# With covariates sex and age                                    
cov.elnet.atyp.other.performance <- CovElasticNetPerformance(dat, metrange, outcome = "atyp.other", 
                                                      repeats = 100, folds = 5, 
                                                      cov.elnet.atyp.other.results)
cov.elnet.spneu.other.performance <- CovElasticNetPerformance(dat, metrange, outcome = "spneu.other", 
                                                       repeats = 100, folds = 5, 
                                                       cov.elnet.spneu.other.results)
cov.elnet.bac.vir.performance <- CovElasticNetPerformance(dat, metrange, outcome = "bac.vir", 
                                                   repeats = 100, folds = 5, 
                                                   cov.elnet.bac.vir.results)

## ROC curves of cross validated elastic net models ---------------------------
source("scripts/ROC_curve_generation_code.R")

## VIP score plot of atyp-other elastic net regression model-------------------
source("scripts/VIP_plot_generation_code.R")

## Performance table ----------------------------------------------------------
performance.overview <- data.frame(Comparison = character(), Model = character(),
                                   Variables = character(), AUC = character(), 
                                   Sensitivity =  character(), 
                                   Specificity = character(), 
                                   BER = character())
performance.overview <- rbind(performance.overview,
                              logit.glycylglycine[[4]],
                              logit.sdma[[4]],
                              logit.lpi18.1[[4]], 
                              logit.age.sex[[4]], 
                              logit.sigatyp[[4]], 
                              logit.sigatyp.age.sex[[4]], 
                              elnet.atyp.other.performance[[1]], 
                              cov.elnet.atyp.other.performance[[1]], 
                              elnet.spneu.other.performance[[1]], 
                              cov.elnet.spneu.other.performance[[1]], 
                              elnet.bac.vir.performance[[1]], 
                              cov.elnet.bac.vir.performance[[1]])
performance.overview$Comparison <- c(rep("Atypical – (S. pneumoniae + viral)", 8), 
                                     rep("S. pneumoniae – (atypical + viral)", 2), 
                                     rep("Viral – (S. pneumoniae + atypical)", 2))

## Figures hyperparameter selection and BER per number of variables selected -----
#  (supplementary material) 
source("scripts/Hyperparameter_selection_figure_code.R")




