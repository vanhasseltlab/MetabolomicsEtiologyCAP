ElasticNetCV <- function(dat, metrange, outcome, alphas, lambdas, repeats, folds){ 

  # The function ElasticNetCV performes a nested, 5-fold, cross validation 
  #  of an elastic net model on metabolomics data.
  
  ## Source dependencies
  source("functions/MakeFolds.R")
  
  ## List to store dataframes and lists per repeat
  repeat_data <- list()
  
  # For r repeats
  for (r in 1:repeats) {    
    print(paste("repeat = ", r))
    
    ## List to store data of optimization of alpha and lambda inside nested loop
    data_f_a_l_BER <- list()
    
    ## Dataframe to store data of optimizated of alpha and lambda 
    optimal_parameters <- as.data.frame(matrix(data = NA, nrow = folds, ncol = 4))
    names(optimal_parameters) <- c("fold", "alpha", "lambda", "BER")
    
    ## List to store elastic net model and its coefficients
    elnet_model <- list()
    elnet_coef <- list()
    
    ## Dataframe to store prediction data from elastic net model
    elnet_data <- as.data.frame(matrix(data = NA, nrow = folds, ncol = 7))
    # In caret confusion matrix: TP = A, FP = B, FN = C, TN = D
    names(elnet_data) <- c("fold", "TP", "FP", "FN", "TN", "BER", "AUC")
    
    # Dataframe to store ROC data
    rocdata.list <- list()

    ## Dataframe to store real data & predictions of testfold elastic net model
    elnet_pred <- as.data.frame(matrix(data = NA, nrow = nrow(dat), ncol = 3))
    names(elnet_pred) <- c("fold", "y", "y_hat")
    
    ## n-fold nested cross validation for one comparison
    # Make 5 folds
    flds <- MakeFolds(r, folds, dat, outcome, "sample.id")
    
    for (f in 1:folds){
      # counter for debugging  
      print(paste("f = ", f))
      
      # Make test and trainfold
      testfold <- filter(dat, as.integer(sample.id) %in% flds[[f]])
      trainfold <- filter(dat, !(as.integer(sample.id) %in% flds[[f]]))
      
      # Make dataframe to store predictions in optimization alpha and lambda
      nrows_al <- length(alphas)*length(lambdas)*(folds-1)
      dfpred2 <- as.data.frame(matrix(data = NA, nrow = nrows_al, ncol = 7))  
      names(dfpred2) <- c("f2", "a", "l", "TP", "FP", "FN", "TN")
      dfpred2_count <- 0
      
      # CV to find optimal values for alpha and lambda by using the Balanced Error Rate (BER)
      for (f2 in (1:folds)[(1:folds) != f]) {
        # counter for debugging  
        print(paste("f2 = ", f2))
        
        testfold2 <- filter(dat, as.integer(sample.id) %in% flds[[f2]])
        
        trainfold2 <- filter(dat, (!(as.integer(sample.id) %in% flds[[f2]]) & 
                                   !(as.integer(sample.id) %in% flds[[f]])))
        
        for (a in alphas){
          for (l in lambdas){
            #Counter for saving data
            dfpred2_count <- dfpred2_count + 1
            
            # Build elastic net model with alpha and labda, save variables nessecary for 
            #  balanced error rate determination (confusion matrix)
            elnet2 <- glmnet(x = as.matrix(trainfold2[, metrange]), 
                             y = as.matrix(trainfold2[, outcome]), 
                             family = "binomial", alpha = a, lambda = l)
            
            p2 <- predict(elnet2, as.matrix(testfold2[, metrange]), type = "response")
            
            predclasses2 <- factor((p2 > 0.5) + 1, levels = c(1, 2))
            
            # confustion matrix
            cmatrix2 <- caret::confusionMatrix(predclasses2,
                                               as.factor(as.numeric(as.factor(testfold2[, outcome]))))
            
            #Save in dataframe: nr correctly predicted, nr wrongly predicted, fold size
            # Caret confusion matrix:
            #  TP = A = matrix[1,1], FP = B = matrix[1,2], 
            #  FN = C = matrix[2, 1], TN = D = matrix[2, 2]
            dfpred2[dfpred2_count, ] <- c(f2, a, l, 
                                          TP = cmatrix2$table[1, 1], FP = cmatrix2$table[1, 2], 
                                          FN = cmatrix2$table[2, 1], TN = cmatrix2$table[2, 2])
          }
        }
        
        ## Determine optimal alpha and lamda for trainingfold
        # Summarize the number of cases right and wrongly predicted for 
        #  equal alpha's and labda's (over different folds)
        # Keep alpha's and labda's in the dataframe as reference
        dfpred2_sum <- dfpred2[complete.cases(dfpred2), ] %>%
          group_by(a, l) %>%
          summarise(TP_sum = sum(TP), FP_sum = sum(FP), FN_sum = sum(FN), TN_sum = sum(TN)) %>%
          ungroup(a, l)
        
        # Calculate BER per alpha and lambda combination
        dfpred2_BER <- dfpred2_sum %>%
          mutate(BER = 0.5 * ((FP_sum / (TN_sum + FP_sum)) + (FN_sum / (FN_sum + TP_sum)))) %>%
          select(a, l, BER)
        
        # Save dfpred2_BER per fold
        data_f_a_l_BER[[f]] <- dfpred2_BER
        
        # Select optimal lambda (biggest at minimal BER) and alpha 
        # (set to max but for no reason other than that we need to pick one value)
        BER_opt <- min(dfpred2_BER$BER)
        l_opt <- max(dfpred2_BER$l[dfpred2_BER$BER == BER_opt])
        a_opt <- max(dfpred2_BER$a[dfpred2_BER$BER == BER_opt & dfpred2_BER$l == l_opt ])
        
        # Save optima for fold into dataframe
        optimal_parameters[f, ] <- c(fold = f, alpha = a_opt, 
                                     lambda = l_opt, BER = BER_opt)
        
        
        ## Build elastic net model to predict testfold
        # Build elastic net model with optimized alpha and labda, save variables nessecary for 
        #  balanced error rate determination
        elnet <- glmnet(x = as.matrix(trainfold[, metrange]), 
                        y = as.matrix(trainfold[, outcome]), 
                        family = "binomial", alpha = a_opt, lambda = l_opt)
        # Predict testfold with chosen alpha and lambda
        p <- predict(elnet, as.matrix(testfold[, metrange]), type = "response")
        predclasses <- factor((p > 0.5) + 1, levels = c(1, 2))
        # confustion matrix
        cmatrix <- caret::confusionMatrix(predclasses, 
                                          as.factor(as.numeric(as.factor(testfold[, outcome]))))
        
        # Calculate AUC
        pred <- prediction(p, testfold[, outcome])
        auc <- performance(pred, "auc")@y.values[[1]] #auc = area under curve
        
        
        # Calculate true and false positive rate
        perf <- performance(pred, "tpr", "fpr")
        # Extract roc data 
        rocdata <- data.frame("fpr" = perf@x.values[[1]], 
                              "tpr" = perf@y.values[[1]])
        
        # Store ROC data
        rocdata.list[[f]] <- rocdata
        
        # Store elastic net model inlist
        elnet_model[[f]] <- elnet
        
        # Store coefficients in list
        elnet_coef[[f]] <- elnet$beta
        
        # Store predictions in dataframe: comparison, testfoldnr, 
        # Caret confusion matrix:
        #  TP = A = matrix[1,1], FP = B = matrix[1,2], 
        #  FN = C = matrix[2, 1], TN = D = matrix[2, 2]
        elnet_data[f, ] <- c(fold = f, 
                              TP = cmatrix$table[1, 1], FP = cmatrix$table[1, 2], 
                              FN = cmatrix$table[2, 1], TN = cmatrix$table[2, 2], 
                              BER = 0.5*((cmatrix$table[1, 2] / (cmatrix$table[2, 2] + cmatrix$table[1, 2])) 
                                          + (cmatrix$table[2, 1] / (cmatrix$table[2, 1] + cmatrix$table[1, 1]))),
                              AUC = auc)
        
        #Save real class and predictions in dataframe
        for ( prednr in 1:nrow(dat)) {
          elnet_pred[prednr, ] <- c(fold = f, 
                                    y = as.factor(as.numeric(as.factor(testfold[, outcome][prednr]))), 
                                    y_hat = predclasses[prednr])
        }
      }
    }
    
    # Save all dataframes and lists in repeat_data list
    repeat_data[[r]] <- list(data_f_a_l_BER, optimal_parameters, elnet_model, 
                             elnet_coef, elnet_data, pred, elnet_pred, rocdata.list)
    names(repeat_data[[r]]) <- c("data_f_a_l_BER", "optimal_parameters", "elnet_model", 
                                 "elnet_coef", "elnet_data", "pred", "elnet_pred", "roc_data")
  }
  return(repeat_data)
}
