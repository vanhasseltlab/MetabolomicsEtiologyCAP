## Supplement code data analysis: linear regression of covariates

# This is code to replicate the data analysis steps from the paper:
# "Metabolomic profiling of microbial disease etiology in community-acquired pneumonia"  
# Code developed by Ilona den Hartog, Laura Zwep, and Coen van Hasselt. 

# Code contents:
# - Data analysis for the comparisons:
#   - Comparison 1: S. pneumoniae versus atypical bacteria + virusses
#   - Comparison 2: Atypical bacteria versus S. pneumoniae + virusses
#   - Comparison 3: Virusses versus S. pneumoniae + atypical bacteria
# 
# - Performed data analysis:
#   - Linear regression with only the covariates sex and age for comparison of 
#     1 pathogen group versus the remaining 2 pathogen groups. 


## Setup R --------------------------------------------------------------------------------------
# R version: 3.6.3
# Rstudio version: 1.1.463

# Clear global environment
rm(list=ls())

# Set working directory and datafolder
setwd("~/ZonMW/Supplement")
datafolder <- "~/ZonMW/Supplement/"

# Load libraries
library(dplyr)
library(glmnet)
library(caret)

## Input ---------------------------------------------------------------------------------------
# Load metabolomics data after quality control
diagnosisdata <- read.csv2(paste0(datafolder, "Data/", "02_diagnosisdata_pretreated.csv"))

# Define metabolite and metadata range in data after quality control
# Columns before metabolite data are defined as 'begin'
begin = 1
# Columns after metabolite data are defined as 'end'
end = (length(diagnosisdata)-3):(length(diagnosisdata))
# Columns containing the metabolite data are defined as 'metrange'
metrange <- 2:(length(diagnosisdata)-4)

# Define data and group identifiers per comparison
spneu_rest <- mutate(diagnosisdata, test_groups = ifelse(diagnosisdata$pathogen_group == "S.pneumoniae", yes = "S.pneumoniae", no = "Atypical or viral"))
spneu_rest$test_groups <- factor(spneu_rest$test_groups, levels = c("S.pneumoniae", "Atypical or viral"))
atyp_rest <- mutate(diagnosisdata, test_groups = ifelse(diagnosisdata$pathogen_group == "atypical", yes = "Atypical", no = "S.pneumoniae or viral"))
atyp_rest$test_groups <- factor(atyp_rest$test_groups, levels = c("Atypical", "S.pneumoniae or viral"))
vir_bac <- mutate(diagnosisdata, test_groups = ifelse(diagnosisdata$pathogen_group == "viral", yes = "Viral", no = "Bacterial"))
vir_bac$test_groups <- factor(vir_bac$test_groups, levels = c("Viral", "Bacterial"))

## Data analysis: Linear regression of covariates -------------------------------------

## List to store dataframes and lists per repeat
repeat_data <- list()
repeatnr = 0

#For 100 repeats (s in 1:100)
for (s in 1:100){    
  repeatnr = repeatnr + 1
  print(paste("repeat =", repeatnr))
  
  ## List to store linear model
  lin_model <- list()
  ## Dataframe to store prediction data from linear model
  lin_data <- as.data.frame(matrix(data = NA, nrow = 15, ncol = 6))
  names(lin_data) <- c("comparison", "fold", "ca", "cb", "cc", "cd")
  
  # Reset i_count
  i_count = 0
  
  # Set seed
  set.seed(s)
  
  # 5 fold cross validation for comparison of all groups, if checking only one comparison:
  #     i = spneu_rest #Set i_count to 0
  #     i = atyp_rest #Set i_count to 1
  #     i = vir_bac #Set i_count to 2
  for (i in list(spneu_rest, atyp_rest, vir_bac)){
    #counter for group comparison  
    i_count = i_count + 1
    print(paste("comparison = ", i_count))
    
    # Adjust lists to add comparison nr:
    lin_model[[i_count]] <- list()
    
    # Determine foldsize per group (20% of samples in testfold, 5 folds)
    tsize_test <- round(sum(i$test_groups == levels(i$test_groups)[1])*0.2)
    tsize_rest <- round(sum(i$test_groups == levels(i$test_groups)[2])*0.2)
    
    # Make 5 stratified folds and save in folds list
    flds <- list()
    
    # Fold 1 - 4
    rest <- i
    for (n in 1:4){
      # Sample for testset by MNC name per group
      tfold_test <- sample(rest$sampleID[rest$test_groups == levels(i$test_groups)[1]], size = tsize_test)
      tfold_rest <- sample(rest$sampleID[rest$test_groups == levels(i$test_groups)[2]], size = tsize_rest)
      # Save in folds list
      flds[[n]] <- c(tfold_test, tfold_rest)
      # gather remaining samples
      rest <- rest %>% 
        filter(!(sampleID %in% tfold_test | sampleID %in% tfold_rest))
    }
    
    # Fold 5 (save all remaining samples in folds list
    flds[[5]] <- as.integer(rest$sampleID)
    # Name folds
    names(flds) <- c("Fold1", "Fold2", "Fold3", "Fold4", "Fold5")
    
    
    for (f in 1:5){
      # counter for debugging  
      print(paste("f = ", f))
      
      testfold <- i %>%
        filter(as.integer(sampleID) %in% flds[[f]])
      
      trainfold <- i %>%
        filter(!(as.integer(sampleID) %in% flds[[f]]))
      
      
      ## Build generalized linear model to predict testfold
      #Build generalized linear model on trainfold, save variables nessecary for Balanced error rate determination
      lin_model <- glm(test_groups ~ leeftijd + geslacht, data = trainfold, family = "binomial")
      
      # Predict testfold with linear model
      p <- predict(lin_model, newdata = testfold, type = "response")
      predclasses <- factor((p>0.5)+1, levels = c(1,2))
      
      # confustion matrix
      cmatrix <- caret::confusionMatrix(predclasses, as.factor(as.numeric(testfold$test_groups)))
      
      # Store linear model inlist
      lin_model[[i_count]][[f]] <- list(lin_model)
      
      # Store predictions in dataframe: comparison, testfoldnr, nr correctly predicted (ca, cd), nr wrongly predicted (cb, cc)
      op_count <- (i_count-1)*5 + f
      lin_data[op_count,] <- c(comparison = i_count, fold = f, ca = cmatrix$table[1,1], cb = cmatrix$table[2,1], cc = cmatrix$table[1,2], cd = cmatrix$table[2,2])
      
    }
  }
  
  # Save all dataframes and lists in repeat_data list
  repeat_data[[repeatnr]] <- list(lin_model, lin_data)
  names(repeat_data[[repeatnr]]) <- c("linear_model", "linear_model_data")
  
}


## Output --------------------------------------------------------------------------------------

# Create results folder (if not already excisting)
# dir.create(paste(datafolder,"Results", sep=""))

# Save results elastic net regressions
save(repeat_data, file = paste(datafolder, "Results/", "05_Linear_regression_of_covariates_results", sep=""))






