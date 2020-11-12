LogisticRegressionFolds <- function(dat, variable, outcome, repeats, folds){ 
  
  # This function performs a logistic regression on the data and outputs
  #  the model, auc value and roc curve

  # Source functions
  source("functions/MakeFolds.R")
  
  # Make dataframe to save AUCs
  logit.auc.data <- as.data.frame(matrix(nrow = repeats*folds, ncol = 1))
  colnames(logit.auc.data) <- "AUC"
  
  # Make list to save ROC plot values
  logit.roc.data <- list()
  
  # Make list saving performance data
  logit.data <- as.data.frame(matrix(data = NA, nrow = repeats*folds, ncol = 4))
  names(logit.data) <- c("Sensitivity", "Specificity", "BER", "AUC")
  
  k = 0
  # For 100 repeats
  for (s in 1:100){   
    
    # Make 5 folds
    flds <- MakeFolds(s, 5, dat, outcome, "sample.id")
    
    logit_dats <- as.data.frame(matrix(data = NA, nrow = folds, ncol = 7))
    # In caret confusion matrix: TP = A, FP = B, FN = C, TN = D
    names(logit_dats) <- c("fold", "TP", "FP", "FN", "TN", "BER", "AUC")
    
    for (f in 1:5){
      k = k + 1
      
      # Select train and testfold
      testfold <- filter(dat, as.integer(sample.id) %in% flds[[f]])
      trainfold <- filter(dat, !(as.integer(sample.id) %in% flds[[f]]))
      
      # Run model over data
      variables <- paste("as.factor(", outcome, ") ~", paste(variable, collapse = " + "))  
      logit <- glm(as.formula(variables), data = trainfold, family = "binomial")
      
      # Apply model to data
      prob <- predict(logit, testfold, type = "response")
      predclasses <- factor((prob > 0.5) + 1, levels = c(1, 2))
      #confusion matrix
      cmatrix <- caret::confusionMatrix(predclasses, 
                                        as.factor(as.numeric(as.factor(testfold[, outcome]))))
      pred <- prediction(prob, testfold[, outcome])
      
      # Calculate AUC (area under curve)
      auc.measure <- performance(pred, "auc")
      auc <- auc.measure@y.values[[1]]
      logit.auc.data[k, ] <- auc

      # Calculate sensitivity, specificity and BER
      # Caret confusion matrix:
      #  TP = A = matrix[1,1], FP = B = matrix[1,2], 
      #  FN = C = matrix[2, 1], TN = D = matrix[2, 2]
      logit_dats[k, ] <- c(iteration = k, 
                           TP = cmatrix$table[1, 1], FP = cmatrix$table[1, 2], 
                           FN = cmatrix$table[2, 1], TN = cmatrix$table[2, 2], 
                           BER = 0.5*((cmatrix$table[1, 2] / (cmatrix$table[2, 2] + cmatrix$table[1, 2])) 
                                      + (cmatrix$table[2, 1] / (cmatrix$table[2, 1] + cmatrix$table[1, 1]))),
                           AUC = auc)
      df <- logit_dats[k, ]
      df$Sensitivity <- df$TP / (df$TP + df$FN)
      df$Specificity <- df$TN / (df$FP + df$TN)
      logit.data[k, ] <- subset(df, select = c("Sensitivity", "Specificity", "BER", "AUC"))
      
            
      # Calculate true and false positive rate
      perf <- performance(pred, "tpr", "fpr")
      
      # Plot ROC curve
      rocdata <- data.frame("fpr" = perf@x.values[[1]],
                             "tpr" = perf@y.values[[1]])
      logit.roc.data[[k]] <- rocdata
    }
  }
  
  # Mean and SD of logit data
  logit.data.mean <- logit.data %>% 
    summarise_if(is.numeric, mean) 
  logit.data.sd <- logit.data %>%
    summarise_if(is.numeric, sd)

  performance.overview <- 
    data.frame(Comparison = outcome, 
               Model = "Logistic Regression", 
               Variables = paste(variable, collapse = " + "),
               AUC = paste(signif(logit.data.mean[[4]], digits = 2), 
                           " (", signif(logit.data.sd[[4]], digits = 2), ")", sep=""),
               Sensitivity = paste(signif(logit.data.mean[[1]], digits = 2), 
                               " (", signif(logit.data.sd[[1]], digits = 2), ")", sep=""),
               Specitivity = paste(signif(logit.data.mean[[2]], digits = 2),
                               " (",signif(logit.data.sd[[2]], digits = 2), ")", sep=""),
               BER = paste(signif(logit.data.mean[[3]], digits = 2), 
                       " (", signif(logit.data.sd[[3]], digits = 2), ")", sep=""))

  # Mean (SD) AUC
  logit.auc.mean <- mean(logit.auc.data[,1])
  
  # Mean ROC values
  total.logit.roc.data <- bind_rows(logit.roc.data, .id = "repetition")
  
  # Plot ROC curve
  mean.logit.roc.data <- tapply(total.logit.roc.data$tpr, total.logit.roc.data$fpr, mean)
  mean.logit.roc.data <- data.frame(fpr.mean = as.numeric(names(mean.logit.roc.data)), tpr.mean = mean.logit.roc.data)
  
  sd.logit.roc.data <- tapply(total.logit.roc.data$tpr, total.logit.roc.data$fpr, sd)
  sd.logit.roc.data <- data.frame(fpr.sd = as.numeric(names(sd.logit.roc.data)), tpr.sd = sd.logit.roc.data)
  
  logit.rocplotdata <- data.frame(mean.logit.roc.data, sd.logit.roc.data)
  logit.roc.plot <- ggplot(data = logit.rocplotdata)+
    geom_ribbon(aes(x = fpr.mean, ymin = tpr.mean - tpr.sd, ymax = tpr.mean + tpr.sd), 
                alpha = 0.4, color = "gray85", fill = "gray85")+
    geom_step(aes(x = fpr.mean, y = tpr.mean), color = "black", size = 0.5, show.legend = F)+
    geom_segment(aes(x = 0, y = 0, xend = 1, yend = 1), 
                 linetype = "dashed", color = "gray", size = 0.25, alpha = 0.4)+
    labs(x = "False Positive Rate", y = "True Positive Rate", 
         title = paste(paste(variable, collapse = " + "), "\n", "Logistic regression", sep = ""))+
    theme_bw()
  
  return(list(logit.auc.mean, logit.roc.plot, logit.rocplotdata, performance.overview))
}
