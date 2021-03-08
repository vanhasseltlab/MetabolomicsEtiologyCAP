PCAplotting <- function(dat, metrange, outcome, colorpalette = "three"){
  
  # PCAplot is a function that makes a PCA plot including amount of variance 
  # explained for the first two principal components. 
  
  # perform PCA
  pca.data <- prcomp(dat[, metrange], scale = FALSE)
  
  # Extract proportion of variance explained for first six components
  variance <- round(summary(pca.data)$importance[2, 1:6]*100, 1)
  
  # Make dataset for plotting
  plot.data <- cbind(as.data.frame(pca.data$x[,1:6]), dat[, outcome])
  names(plot.data) <- c("comp_1", "comp_2", "comp_3","comp_4", "comp_5", "comp_6", outcome)
  
  # Make PCA plot
  pca.plot.c12 <- ggplot(plot.data, aes(y = comp_2, x = comp_1, colour = dat[, outcome]))+
    geom_point(aes(), size = 2, alpha = 0.8)+
    scale_color_lei(palette = colorpalette)+
    labs(x = paste("PC 1 (", variance[1], "% )"),
         y = paste("PC 2 (", variance[2], "% )"),
         colour = outcome, shape = outcome) +
    theme_bw()
  
  return(pca.plot.c12)
}