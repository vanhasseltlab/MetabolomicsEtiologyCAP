ElnetROC <- function(dat, elnet_results, repeats, nfolds){
  ## Script to generate a ROC curve of elastic net models
  
  
  # Make list to combine ROC data in for plotting
  elnet.roc.data <- list()
  
  k = 0
  for (r in 1:repeats){
    for (f in 1:nfolds){
      k = k + 1
      
      # Extract ROC data from elnet results
      elnet.roc.data[[k]] <- elnet_results[[r]]$roc_data[[f]]
  }}
  
  
  # Mean ROC values
  total.elnet.roc.data <- bind_rows(elnet.roc.data, .id = "repetition")
  
  # Plot ROC curve
  mean.elnet.roc.data <- tapply(total.elnet.roc.data$tpr, total.elnet.roc.data$fpr, mean)
  mean.elnet.roc.data <- data.frame(fpr.mean = as.numeric(names(mean.elnet.roc.data)), tpr.mean = mean.elnet.roc.data)
  
  sd.elnet.roc.data <- tapply(total.elnet.roc.data$tpr, total.elnet.roc.data$fpr, sd)
  sd.elnet.roc.data <- data.frame(fpr.sd = as.numeric(names(sd.elnet.roc.data)), tpr.sd = sd.elnet.roc.data)
  
  elnet.rocplotdata <- data.frame(mean.elnet.roc.data, sd.elnet.roc.data)
  elnet.roc.plot <- ggplot(data = elnet.rocplotdata)+
    geom_ribbon(aes(x = fpr.mean, ymin = tpr.mean - tpr.sd, ymax = tpr.mean + tpr.sd), 
                alpha = 0.4, color = "gray85", fill = "gray85")+
    geom_step(aes(x = fpr.mean, y = tpr.mean), color = "black", size = 0.5, show.legend = F)+
    geom_segment(aes(x = 0, y = 0, xend = 1, yend = 1), 
                 linetype = "dashed", color = "gray", size = 0.25, alpha = 0.4)+
    labs(x = "False Positive Rate", y = "True Positive Rate", 
         title = paste("All metabolites", "\n", "elastic net", sep = ""))+
    theme_bw()
  
  return(list(elnet.roc.plot, elnet.rocplotdata))
}



