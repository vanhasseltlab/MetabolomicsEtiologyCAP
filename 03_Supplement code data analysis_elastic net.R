## Supplement code data analysis: elastic net

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
#   - Elastic net for comparison of 1 pathogen group versus the remaining 2 pathogen groups. 


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

## Data analysis: Elastic net regression (without covariates). Note: VERY LONG RUNTIME!-------------------------------------

# Define sequence of alpha's and labda's
alphas <- seq(0, 1, 0.05) 
labdas <- seq(0, 1, 0.05)

## List to store dataframes and lists per repeat
repeat_data <- list()
repeatnr = 0

#For 100 repeats (s in 1:100)
for (s in 1:100){    
  repeatnr = repeatnr + 1
  print(paste("repeat =", repeatnr))
  
  ## List to store data of optimization of alpha and lambda inside nested loop
  data_f_a_l_BER <- list()
  ## Dataframe to store data of optimizated of alpha and lambda 
  optimal_parameters <- as.data.frame(matrix(data = NA, nrow = 15, ncol = 5))
  names(optimal_parameters) <- c("comparison", "fold", "alpha", "lambda", "BER")
  ## List to store elastic net model
  elnet_model <- list()
  ## List to store coefficients of elastic net model
  elnet_coef <- list()
  ## Dataframe to store prediction data from elastic net model
  elnet_data <- as.data.frame(matrix(data = NA, nrow = 15, ncol = 6))
  names(elnet_data) <- c("comparison", "fold", "ca", "cb", "cc", "cd")
  ## Dataframe to store real data & predictions of testfold elastic net model
  elnet_pred <- as.data.frame(matrix(data = NA, nrow = 252, ncol = 4))
  names(elnet_pred) <- c("comparison", "fold", "y", "y_hat")
  
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
    data_f_a_l_BER[[i_count]] <- list()
    elnet_model[[i_count]] <- list()
    elnet_coef[[i_count]] <- list()
    
    # Determine foldsize per group (20% of samples in testfold, 5 folds)
    tsize_test <- round(sum(i$test_groups == levels(i$test_groups)[1])*0.2)
    tsize_rest <- round(sum(i$test_groups == levels(i$test_groups)[2])*0.2)
    
    # Make 5 stratified folds and save in folds list
    flds <- list()
    
    # Fold 1 - 4
    rest <- i
    for (n in 1:4){
      # Sample for testset by MNC name per group
      tfold_test <- sample(rest$SampleID[rest$test_groups == levels(i$test_groups)[1]], size = tsize_test)
      tfold_rest <- sample(rest$SampleID[rest$test_groups == levels(i$test_groups)[2]], size = tsize_rest)
      # Save in folds list
      flds[[n]] <- c(tfold_test, tfold_rest)
      # gather remaining samples
      rest <- rest %>% 
        filter(!(SampleID %in% tfold_test | SampleID %in% tfold_rest))
    }
    
    # Fold 5 (save all remaining samples in folds list
    flds[[5]] <- as.integer(rest$SampleID)
    # Name folds
    names(flds) <- c("Fold1", "Fold2", "Fold3", "Fold4", "Fold5")
    
    
    for (f in 1:5){
      # counter for debugging  
      print(paste("f = ", f))
      
      testfold <- i %>%
        filter(as.integer(SampleID) %in% flds[[f]])
      
      trainfold <- i %>%
        filter(!(as.integer(SampleID) %in% flds[[f]]))
      
      # dataframe to store predictions in optimization alpha and lambda
      nrows_al <- length(alphas) * length(labdas) * 4  # 4 folds
      dfpred2 <- as.data.frame(matrix(data = NA, nrow = nrows_al, ncol = 7))  # nrows depend on nr of alpha and lambda combinations and number of folds
      names(dfpred2) <- c("f2", "a", "l", "ca", "cb", "cc", "cd")
      dfpred2_count <- 0
      
      # CV to find optimal values for alpha and lambda by using the Balanced Error Rate (BER)
      for (f2 in (1:5)[(1:5)!=f]){
        # counter for debugging  
        print(paste("f2 = ", f2))
        
        testfold2 <- i %>%
          filter(as.integer(SampleID) %in% flds[[f2]])
        
        trainfold2 <- i %>%
          filter(!(as.integer(SampleID) %in% flds[[f2]])) %>%
          filter(!(as.integer(SampleID) %in% flds[[f]]))
        
        for (a in alphas){
          for (l in labdas){
            #Counter for saving data
            dfpred2_count <- dfpred2_count + 1
            
            #Build elastic net model with alpha and labda, save variables nessecary for Balanced error rate determination (confusion matrix, ROC curve and accuracy curve. )
            elnet2 <- glmnet(x= as.matrix(trainfold2[, metrange]), y= as.matrix(trainfold2[, "test_groups"]), family = "binomial",alpha = a, lambda = l)
            
            p2 <- predict(elnet2, as.matrix(testfold2[, metrange]), type = "response")
            
            predclasses2 <- factor((p2>0.5)+1, levels = c(1,2))
            
            # confustion matrix
            cmatrix2 <- caret::confusionMatrix(predclasses2, as.factor(as.numeric(testfold2$test_groups)))
            
            #Save in dataframe: nr correctly predicted, nr wrongly predicted, fold size
            dfpred2[dfpred2_count,] <- c(f2, a, l, ca = cmatrix2$table[1,1], cb = cmatrix2$table[2,1], cc = cmatrix2$table[1,2], cd = cmatrix2$table[2,2])
          }
        }
      }
      
      ## Determine optimal alpha and lamda for trainingfold
      # Summarize the number of cases right and wrongly predicted for equal alpha's and labda's (over different folds)
      # Keep alpha's and labda's in the dataframe as reference
      dfpred2_sum <- dfpred2[complete.cases(dfpred2),] %>%
        group_by(a, l) %>%
        summarise(ca_sum = sum(ca), cb_sum = sum(cb), cc_sum = sum(cc), cd_sum = sum(cd)) %>%
        ungroup(a, l)
      
      # Calculate BER per alpha and lambda combination
      dfpred2_BER <- dfpred2_sum %>%
        mutate(BER = 0.5 * ( cb_sum / (ca_sum + cb_sum) + cc_sum / (cc_sum + cd_sum))) %>%
        select(a, l, BER)
      
      # Save dfpred2_BER per fold
      data_f_a_l_BER[[i_count]][[f]] <- dfpred2_BER
      
      # Select optimal lambda (biggest at minimal BER) and alpha (set to max but for no reason other than that we need to pick one value)
      BER_opt <- min(dfpred2_BER$BER)
      l_opt <- max(dfpred2_BER$l[dfpred2_BER$BER == BER_opt])
      a_opt <- max(dfpred2_BER$a[dfpred2_BER$BER == BER_opt & dfpred2_BER$l == l_opt ])
      
      # Save optima for fold into dataframe
      op_count <- (i_count-1)*5 + f
      optimal_parameters[op_count,] <- c(comparison = i_count, fold = f, alpha = a_opt, lambda = l_opt, BER = BER_opt)
      
      
      ## Build elastic net model to predict testfold
      #Build elastic net model with optimized alpha and labda, save variables nessecary for Balanced error rate determination
      elnet <- glmnet(x= as.matrix(trainfold[, metrange]), y= as.matrix(trainfold[, "test_groups"]), family = "binomial",alpha = a_opt, lambda = l_opt)
      # Predict testfold with chosen alpha and lambda
      p <- predict(elnet, as.matrix(testfold[, metrange]), type = "response")
      predclasses <- factor((p>0.5)+1, levels = c(1,2))
      # confustion matrix
      cmatrix <- caret::confusionMatrix(predclasses, as.factor(as.numeric(testfold$test_groups)))
      
      # Store elastic net model inlist
      elnet_model[[i_count]][[f]] <- elnet
      
      # Store coefficients in list
      elnet_coef[[i_count]][[f]] <- elnet$beta
      
      # Store predictions in dataframe: comparison, testfoldnr, nr correctly predicted (ca, cd), nr wrongly predicted (cb, cc)
      elnet_data[op_count,] <- c(comparison = i_count, fold = f, ca = cmatrix$table[1,1], cb = cmatrix$table[2,1], cc = cmatrix$table[1,2], cd = cmatrix$table[2,2])
      #Save real class and predictions in dataframe
      for ( prednr in 1:length(predclasses) ) {
        ii <- (length(elnet_pred[!is.na(elnet_pred$comparison),1])+1)
        elnet_pred[ii,] <- c(comparison = i_count, fold = f, y = testfold$test_groups[prednr], y_hat = predclasses[prednr])
      }
      
    }
  }
  
  # Save all dataframes and lists in repeat_data list
  repeat_data[[repeatnr]] <- list(data_f_a_l_BER, optimal_parameters, elnet_model, elnet_coef, elnet_data, elnet_pred)
  names(repeat_data[[repeatnr]]) <- c("data_f_a_l_BER", "optimal_parameters", "elnet_model", "elnet_coef", "elnet_data", "elnet_pred")
  
}


## Output --------------------------------------------------------------------------------------

# Create results folder
dir.create(paste(datafolder,"Results", sep=""))

# Save results elastic net regressions
save(repeat_data, file = paste(datafolder, "Results/", "03_Elastic_net_results.Rdata", sep=""))






