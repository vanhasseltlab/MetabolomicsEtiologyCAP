ElasticNetPerformance <- function(dat, metrange, outcome, repeats, folds, repeat_data){

  #Save all data from all folds and repeats into one dataframe
  innerCV_oppars_BER <- data.frame(matrix(data=NA, nrow = folds*repeats, ncol = 5))
  names(innerCV_oppars_BER) <- c("repeatnr", "fold", "alpha", "lambda", "BER")
  
  outerCV_sens_spec <- data.frame(matrix(data=NA, nrow = folds*repeats, ncol = 7))
  names(outerCV_sens_spec) <- c("repeatnr", "fold", "Sensitivity", "Specificity", "BER", "AUC", "nvars")
  
  for (r in 1:repeats){
    start <- 1 + (r-1)*folds
    stop <- r*folds
    
    #Store all optimal parameters (a, l, BER) from innerCV for all folds and repeats in one dataframe
    innerCV_oppars_BER[start:stop,] <- cbind("repeatnr" = rep(r, folds), repeat_data[[r]]$optimal_parameters)
    
    #Calculate sensitivity, specificity and BER for all folds and repeats in outer cross validation loop
    df <- repeat_data[[r]]$elnet_data
    df$Sensitivity <- df$TP / (df$TP + df$FN)
    df$Specificity <- df$TN / (df$FP + df$TN)
    df$BER <- 0.5 * (df$FP / (df$TN + df$FP) + (df$FN / (df$FN + df$TP)))
    df2 <- subset(df, select = c("fold", "Sensitivity", "Specificity", "BER", "AUC"))
    
    #Extract number of variables selected in outer CV elastic net model
    vars.data <- data.frame("nvars" = matrix(ncol = 1, nrow = 5))
    for (f in 1:folds){
      vars.data[f,] <- sum(repeat_data[[r]]$elnet_model[[f]]$beta[,1] != 0)
    }
    
    outerCV_sens_spec[start:stop, ] <- cbind("repeatnr" = rep(r, folds), df2, "nvars"= vars.data[,1])
  }
  
  # Sensitivity: mean and sd
  sens.mean <- mean(outerCV_sens_spec$Sensitivity)
  sens.sd <- sd(outerCV_sens_spec$Sensitivity)
  
  # Specificity: mean and sd
  spec.mean <- mean(outerCV_sens_spec$Specificity)
  spec.sd <- sd(outerCV_sens_spec$Specificity)
  
  # BER: mean and sd
  ber.mean <- mean(outerCV_sens_spec$BER)
  ber.sd <- sd(outerCV_sens_spec$BER)
  
  # AUC: mean and sd
  auc.mean <- mean(outerCV_sens_spec$AUC)
  auc.sd <- sd(outerCV_sens_spec$AUC)
  
  # Number of variables: mean and sd
  nvars.mean <- mean(outerCV_sens_spec$nvars)
  nvars.sd <- sd(outerCV_sens_spec$nvars)

  
  performance.overview <- 
    data.frame(Comparison = outcome,
               Model = "Elastic net regression",
               Variables = paste(signif(nvars.mean, digits = 2), 
                                   " (", signif(nvars.sd, digits = 2), ")", sep = ""),
               AUC = paste(signif(auc.mean, digits = 2), " (",
                           signif(auc.sd, digits = 2), ")", sep=""),
               Sensitivity = paste(signif(sens.mean, digits = 2),
                                   " (", signif(sens.sd, digits = 2), ")", sep=""),
               Specitivity = paste(signif(spec.mean, digits = 2), 
                               " (",signif(spec.sd, digits = 2), ")", sep=""),
               BER = paste(signif(ber.mean, digits = 2), " (",
                           signif(ber.sd, digits = 2), ")", sep=""))
                               
  return(list(performance.overview, auc.mean))  
}
