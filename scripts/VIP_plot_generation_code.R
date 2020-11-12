## Supplement code Figure 2

# This is code to replicate Figure 2 the paper:
# "Metabolomic profiling of microbial disease etiology in community-acquired pneumonia"  
# Code developed by Ilona den Hartog, Laura Zwep, and Coen van Hasselt. 

# # Load libraries
# library(tidyr)
# library(cowplot)
# library(table1)



## Input ---------------------------------------------------------------------------------------

# # Load elastic net results
# load("results/elnet.atyp.other.results.R")

# Load metabolite names and metabolite class information of remaining metabolites
metabolite_info <- read_excel("data/Supplementary material II, Table S3, metabolite names and identifiers.xlsx")

## Data analysis: generate information for figure from elastic net results----------------------

## Put all beta's for the 5 folds in one dataframe
# Extract all beta values for the detected metabolites and add them to the 'betas' dataframe
betas <- data.frame(matrix(nrow = length(metrange), ncol = 501))
betas[,1] <- metrange
names(betas)[1] <- c("Metabolite")

c = 2
for (r in 1:100){
    betas[,(c):(c+4)] <- data.frame(
      data.frame(matrix(unlist(elnet.atyp.other.results[[r]]$elnet_model[[1]]$beta))),
      data.frame(matrix(unlist(elnet.atyp.other.results[[r]]$elnet_model[[2]]$beta))),
      data.frame(matrix(unlist(elnet.atyp.other.results[[r]]$elnet_model[[3]]$beta))),
      data.frame(matrix(unlist(elnet.atyp.other.results[[r]]$elnet_model[[4]]$beta))),
      data.frame(matrix(unlist(elnet.atyp.other.results[[r]]$elnet_model[[5]]$beta))))
    names(betas)[c:(c+4)] <- c(paste("Beta_R", r, "_F1", sep = ""),
                               paste("Beta_R", r, "_F2", sep = ""),
                               paste("Beta_R", r, "_F3", sep = ""),
                               paste("Beta_R", r, "_F4", sep = ""),
                               paste("Beta_R", r, "_F5", sep = ""))
    c = c + 5
  }


# Select only metabolites for plotting that have at least one non-zero value in one fold
betas_filtered <- betas[rowSums(betas[, 2:501] > 0) != 0, ]

# To make a bar plot of metabolites percentage influence on prediction (variable importance in prediction, VIP):
# Make a vector with the sum of beta's per fold/repeat
sumB <- data.frame(matrix(nrow = 500, ncol = 1))
names(sumB) <- "sumB"
for (i in 3:ncol(betas_filtered)){
  sumB[(i-2),] <- summarize(betas_filtered, sumBetas = sum(abs(betas_filtered[,i])))
}

# Dataframe with percentage of influence per metabolite per fold/repeat
betas_perc <- betas_filtered
for (j in 3:ncol(betas_filtered)){
  for (i in 1:nrow(betas_filtered)){
    betas_perc[i,j] <- betas_filtered[i,j]/sumB[(j-2),]*100
  }
}

# Arrange the metabolites by importance (higher mean variable importance in prediction)
betas_perc_arr <- betas_perc %>%
  mutate(Beta_perc_mean = rowMeans(betas_perc[,3:ncol(betas_perc)])) %>%
  arrange(desc(Beta_perc_mean))
# Re-arrange factors on mean
betas_perc_arr$Metabolite <- factor(betas_perc_arr$Metabolite, levels = betas_perc_arr$Metabolite[order(betas_perc_arr$Beta_perc_mean, decreasing = FALSE)])

# Select metabolites with percentage mean > 1% 
betas_plotdat2 <- betas_perc_arr[abs(betas_perc_arr$Beta_perc_mean) >= 1,]
betas_plotdat2$Metabolite <- as.character(factor(betas_plotdat2$Metabolite))

# Select columns from metabolite information used in the VIP plot: metabolite class en official metabolite name. 
metabolite_info_2 <- select(metabolite_info, Measurement_platform, Metabolite_class, Targeted_metabolite, Detected_metabolite_name_in_R)
# Subset metabolite info maka a dataset of that can be joined to the betas data of metabolites with VIP > 1%
metabolite_plotdat <- filter(metabolite_info_2, metabolite_info_2$Detected_metabolite_name_in_R %in% betas_plotdat2$Metabolite)
metabolite_plotdat$Detected_metabolite_name_in_R <- as.character(factor(metabolite_plotdat$Detected_metabolite_name_in_R))
metabolite_plotdat$Targeted_metabolite <- as.character(factor(metabolite_plotdat$Targeted_metabolite))
# Join beta's data and metabolite data for plotting
betas_plotdat3 <- left_join(betas_plotdat2, metabolite_plotdat, by = c("Metabolite" = "Detected_metabolite_name_in_R"))
# Drop unused metabolite class levels
betas_plotdat3$Metabolite_class <- factor(betas_plotdat3$Metabolite_class)

# Re-arrange factors on mean again
betas_plotdat3$Targeted_metabolite <- factor(betas_plotdat3$Targeted_metabolite, levels = betas_plotdat3$Targeted_metabolite[order(betas_plotdat3$Beta_perc_mean, decreasing = FALSE)])

# Re-arrange levels of metabolite class for nice plotting
betas_plotdat3$Metabolite_class <- factor(betas_plotdat3$Metabolite_class, levels = c("Amino acids", 
                                                                                      "Biogenic amines", 
                                                                                      "Acylcarnitines", 
                                                                                      "Betaines", 
                                                                                      "Organic acids", 
                                                                                      "Fatty acids", 
                                                                                      "Diacyl-phosphatidylcholine", 
                                                                                      "Diacyl-phosphatidylethanolamine", 
                                                                                      "Sphingomyelin", 
                                                                                      "Ceramides", 
                                                                                      "Lysophospholipids", 
                                                                                      "Bile acids and other steroids", 
                                                                                      "Endocanbinoids"))  

# Convert the data to a long format that can be used to plot beta's with error bar for folds
betas_plotdat3_long <- gather(betas_plotdat3, key = "Repeat/Fold", value = "Beta",  -Metabolite, -Beta_perc_mean, -Measurement_platform, -Metabolite_class, -Targeted_metabolite)


## Figure generation: Boxplot of elastic net VIP values per metabolite -------------------------
# - For comparison atypical bacteria versus S.pneumoniae + virusses
# - Colored by metabolite class

# Boxplot VIP values of selected metabolites
betas_plot_VIP <- ggplot(data = betas_plotdat3_long)+
  geom_boxplot(aes(x = Targeted_metabolite, y = Beta, color = Metabolite_class), outlier.size = 0.25)+
  scale_color_lei(palette = "mixed")+
  labs(x = "Metabolite", y = "Variable importance in prediction (%)", color = "Metabolite class")+
  coord_flip()+
  theme_bw()

## Output --------------------------------------------------------------------------------------
# Save figure
ggexport(betas_plot_VIP, filename = "figures/VIP_plot_atyp_other.png",
         width = 2200, height = 1100, res = 250)
