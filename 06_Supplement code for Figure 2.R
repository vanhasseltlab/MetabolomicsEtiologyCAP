## Supplement code Figure 2

# This is code to replicate the Figure 2 the paper:
# "Metabolomic profiling of microbial disease etiology in community-acquired pneumonia"  
# Code developed by Ilona den Hartog, Laura Zwep, and Coen van Hasselt. 

# Code contents:
# - Creating color palette for figures
# - Input data
# - Data analysis: generate information for figure from elastic net results
# - Figure generation
# - Output figure

## Setup R --------------------------------------------------------------------------------------
# R version: 3.6.3
# Rstudio version: 1.1.463

# Clear global environment
rm(list=ls())

# Set working directory and datafolder
setwd("~/ZonMW/Supplement")
datafolder <- "~/ZonMW/Supplement/"

# Load libraries
library(readxl)
library(tidyr)
library(ggplot2)
library(tidyverse)
library(RColorBrewer)
library(cowplot)
library(ggpubr)
library(table1)

## Create color palette for figures ------------------------------------------------------------
#   Using the colors of Leiden University and Leiden Academic Centre for Drug Research
#   source: https://drsimonj.svbtle.com/creating-corporate-colour-palettes-for-ggplot2

# Create vector with 8 Leiden university / LACDR colors

lei_colors <- c(
  `darkblue`   = "#001158",
  `orange`     = "#FF9933", 
  `red`        = "#be1908",
  `yellow`     = "#F2F45F", # or FFEE33
  `lightorange`= "#FFC233",
  `lightgreen` = "#aaad00",
  `darkgreen`  = "#2c712d",
  `turquoise`  = "#34a3a9",
  `blue`       = "#5cb1eb", 
  `lightblue`  = "#AFD6FC",
  `pink`       = "#FE4ED7",
  `violet`     = "#b02079", 
  `purple`     = "#780B98")

## Function to extract the hex codes from this vector by name
#' Function to extract lei colors as hex codes
#' @param ... Character names of lei_colors 
lei_cols <- function (...){
  cols <- c(...)
  if (is.null(cols))
    return (lei_colors)
  lei_colors[cols]
}

lei_palettes <- list(
  `main`  = lei_cols("darkblue", "orange"),
  `three` = lei_cols("darkblue", "orange", "darkgreen"),
  `cool`  = lei_cols("darkblue", "lightblue", "turquoise", "lightgreen", "darkgreen"),
  `hot`   = lei_cols("violet", "red", "orange"),
  `mixed` = lei_cols("darkblue", "lightblue", "turquoise", "darkgreen", "lightgreen", "orange", "red", "violet"), 
  `extended` = lei_cols("darkblue", "blue", "lightblue", "turquoise", "darkgreen", "lightgreen", "yellow", "lightorange", "orange", "red", "pink", "violet", "purple"))

#' Return function to interpolate a lei color palette
#' @param palette Character name of palette in lei_palettes
#' @param reverse Boolean indicating whether the palette should be reversed
#' @param ... Additional arguments to pass to colorRampPalette(), such as an alpha value
lei_pal <- function(palette = "main", reverse = FALSE, ...) {
  pal <- lei_palettes[[palette]]
  
  if (reverse) pal <- rev(pal)
  
  colorRampPalette(pal, ...)
}

# Now return a function for any palette, for example `cool`:
lei_pal("cool")
# The returned function will interpolate the palette colors for a certain number of levels, making it possible to create shades between our original colors. To demonstrate, we can interpolate the "cool" palette to a length of 10:
lei_pal("cool")(10)
# This is what we need to create custom ggplot2 scales

## Create custom color and fill scales for ggplot2 by creating one function for color and one for fill. 
#' Color scale constructor for lei colors
#' @param palette Character name of palette in lei_palettes
#' @param discrete Boolean indicating whether color aesthetic is discrete or not
#' @param reverse Boolean indicating whether the palette should be reversed
#' @param ... Additional arguments passed to discrete_scale() or
#'            scale_color_gradientn(), used respectively when discrete is TRUE or FALSE
#'
scale_color_lei <- function(palette = "main", discrete = TRUE, reverse = FALSE, ...) {
  pal <- lei_pal(palette = palette, reverse = reverse)
  if (discrete) {
    discrete_scale("colour", paste0("lei_", palette), palette = pal, ...)
  } else {
    scale_color_gradientn(colours = pal(256), ...)
  }
}

#' Fill scale constructor for lei colors
#'
#' @param palette Character name of palette in lei_palettes
#' @param discrete Boolean indicating whether color aesthetic is discrete or not
#' @param reverse Boolean indicating whether the palette should be reversed
#' @param ... Additional arguments passed to discrete_scale() or
#'            scale_fill_gradientn(), used respectively when discrete is TRUE or FALSE
#'
scale_fill_lei <- function(palette = "main", discrete = TRUE, reverse = FALSE, ...) {
  pal <- lei_pal(palette = palette, reverse = reverse)
  
  if (discrete) {
    discrete_scale("fill", paste0("lei_", palette), palette = pal, ...)
  } else {
    scale_fill_gradientn(colours = pal(256), ...)
  }
}


## Input ---------------------------------------------------------------------------------------

# Load elastic net data
load("~/ZonMW/Supplement/Results/03_Elastic_net_results.Rdata")

# Load metabolite names and metabolite class information of remaining metabolites
metabolite_info <- read_excel(paste0(datafolder, "Data/", "Supplementary material II, Table S3, metabolite names and identifiers.xlsx"))

## Data analysis: generate information for figure from elastic net results----------------------
# Extract metabolite names (same for all folds and repeats)
metnames_elnet_model <- data.frame(matrix(unlist(dimnames(repeat_data[[1]]$elnet_model[[1]][[1]]$beta)[[1]])))

## Put all beta's for the 5 folds in one dataframe
# Extract all beta values for the detected metabolites and add them to the 'betas' dataframe
betas <- data.frame(matrix(nrow = 347, ncol = 102))
betas[,1] <- metnames_elnet_model
names(betas)[1] <- c("Metabolite")

c = 2
for (r in 1:100){
  for (i in 2){
    betas[,(c):(c+4)] <- data.frame(
      data.frame(matrix(unlist(repeat_data[[r]]$elnet_model[[i]][[1]]$beta))),
      data.frame(matrix(unlist(repeat_data[[r]]$elnet_model[[i]][[2]]$beta))),
      data.frame(matrix(unlist(repeat_data[[r]]$elnet_model[[i]][[3]]$beta))),
      data.frame(matrix(unlist(repeat_data[[r]]$elnet_model[[i]][[4]]$beta))),
      data.frame(matrix(unlist(repeat_data[[r]]$elnet_model[[i]][[5]]$beta))))
    names(betas)[c:(c+4)] <- c(paste("Beta_R", r, "_C", i, "_F1", sep = ""),
                               paste("Beta_R", r, "_C", i, "_F2", sep = ""),
                               paste("Beta_R", r, "_C", i, "_F3", sep = ""),
                               paste("Beta_R", r, "_C", i, "_F4", sep = ""),
                               paste("Beta_R", r, "_C", i, "_F5", sep = ""))
    c = c + 5
  }
  
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
  geom_boxplot(aes(x = Targeted_metabolite, y = Beta, color = Metabolite_class))+
  scale_color_lei(palette = "extended")+
  labs(x = "Metabolite", y = "Variable importance in prediction (%)", color = "Metabolite class")+
  coord_flip()+
  theme_bw()




## Output --------------------------------------------------------------------------------------

# Create figures folder
dir.create(paste(datafolder,"Figures", sep=""))

# Save figure
svg(filename = "Figures/Figure2.svg", width = 9, height = 4)
betas_plot_VIP
dev.off()
