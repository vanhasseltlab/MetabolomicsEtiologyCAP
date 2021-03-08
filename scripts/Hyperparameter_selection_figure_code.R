## Figure hyperparameter selection (supplementary information)

# Optimization of alpha and lambda to reach minimal balanced error rate (BER)
# 
# * Plot of all alpha and labdas tested in inner CV against mean BER
# * Plot of the optimal alpha and lambda combinations chosen in the inner CV against their BER in the outer CV
# * Amout of variables selected in elastic net model in outer CV

## Generation of data -----
## InnerCV
# Dataframe containing a, l and BER for all folds and repeats in InnerCV
innerCV_BER_a_l <- data.frame(matrix(data = NA, nrow = 441*5*100, ncol = 5))
names(innerCV_BER_a_l) <- c("repeatnr", "fold", "alpha", "lambda", "BER")
for (r in 1:100){
    for (f in 1:5){
      start <- 1 + (f-1)*441 + (r-1)*2205
      stop <- f*441 + (r-1)*2205
      
      df <- data.frame(matrix(data= NA, nrow = 441, ncol = 5))
      names(df) <- c("repeatnr", "fold", "alpha", "lambda", "BER")
      df$repeatnr <- r
      df$fold <- f
      df$alpha <- elnet.atyp.other.results[[r]]$data_f_a_l_BER[[f]]$a
      df$lambda <- elnet.atyp.other.results[[r]]$data_f_a_l_BER[[f]]$l
      df$BER <- elnet.atyp.other.results[[r]]$data_f_a_l_BER[[f]]$BER
      
      innerCV_BER_a_l[start:stop,] <- df
    }
}

#Put all optimal parameters (a, l and BER) from innerCV for all folds and repeats in one dataframe
innerCV_oppars_BER <- data.frame(matrix(data=NA, nrow = 5*100, ncol = 5))
names(innerCV_oppars_BER) <- c("repeatnr", "fold", "alpha", "lambda", "BER")
for (r in 1:100){
  start <- 1 + (r-1)*5
  stop <- r*5
  innerCV_oppars_BER[start:stop,] <- cbind("repeatnr" = rep(r, 5), elnet.atyp.other.results[[r]]$optimal_parameters)
}


## Outer CV
# Calculate BER for all folds and repeats (per alpha and lambda combination resulting from innerCV)
outerCV_BER_a_l <- data.frame(matrix(data=NA, nrow = 5*100, ncol = 3))
names(outerCV_BER_a_l) <- c("repeatnr", "fold", "BER")
for (r in 1:100){
  start <- 1 + (r-1)*5
  stop <- r*5
  
  df <- elnet.atyp.other.results[[r]]$elnet_data
  df$BER <- 0.5 * (df$FP / (df$TN + df$FP) + (df$FN / (df$FN + df$TP)))
  df2 <- subset(df, select = c("fold", "BER"))
  
  outerCV_BER_a_l[start:stop,] <- cbind("repeatnr" = rep(r, 5), df2)
}
outerCV_BER_a_l <- outerCV_BER_a_l %>% 
  mutate(alpha = innerCV_oppars_BER$alpha) %>% 
  mutate(lambda = innerCV_oppars_BER$lambda)


# Make dataframe that includes the amount of variables selected:
outerCV_nvars <- data.frame(matrix(data = NA, nrow = 5*100, ncol = 3))
names(outerCV_nvars) <- c("repeatnr", "fold", "nvars")
for (r in 1:100){
    for (f in 1:5){
      start <- 1 + (f-1) + (r-1)*5
      stop <- f + (r-1)*5
      
      nvars <- as.numeric(as.character(sum(elnet.atyp.other.results[[r]]$elnet_model[[f]]$beta[,1] != 0)))
      outerCV_nvars[start:stop,] <- cbind("repeatnr" = r, "fold" = f, "nvars" = nvars)
    }
}

# Make dataframe that includes a, l, BER and nvars outerCV
outerCV_BER_nvars <- data.frame(cbind(outerCV_BER_a_l, nvars = outerCV_nvars$nvars))

## Generation of figure hyperparameter selection ----
# Plot all tested alpha's against mean BER in inner CV for comparison atyp - s.pneu, color by lambda
repeat_a_l_BER_alll <- innerCV_BER_a_l %>%
  filter(lambda != 0) %>%
  group_by(alpha, lambda) %>%
  mutate(BERmean = mean(BER)) %>%
  mutate(BERsd = sd(BER)) %>%
  ungroup() %>%
  distinct(alpha ,lambda, .keep_all = T)

alpha_plot_alll <- ggplot(repeat_a_l_BER_alll) +
  geom_point(aes(x = alpha, y = BERmean, color = log(lambda), group = log(lambda)))+
  geom_line(aes(x = alpha, y = BERmean, color = log(lambda), group = log(lambda))) +
  scale_color_lei(palette = "cool", discrete = FALSE)+
  labs(x = expression(alpha), y = "Balanced Error Rate", color = expression(paste("log (", lambda, ")")), 
       title = expression(paste("Mean BER of all ", alpha,  " and ",  lambda,  " combinations")))+
  ylim(0, 0.8)+
  theme_bw()+
  theme(plot.title = element_text(size=10), axis.title.x = element_text(size=10), axis.title.y = element_text(size=10))


# plot optimal alpha's and lambda's
outerCV_nvars_plotdata <- outerCV_BER_nvars %>%
  filter(lambda != 0) 

alpha_plot_opt <- ggplot(outerCV_nvars_plotdata) +
  geom_point(aes(x = alpha, y = BER, color = log(lambda), group = (lambda)))+
  scale_color_lei(palette = "cool", discrete = FALSE)+
  labs(x = expression(alpha), y = "Balanced Error Rate", color = expression(paste("log (", lambda, ")")), 
       title = expression(paste("BER of the optimal ", alpha, " and ", lambda, " combinations")))+
  ylim(0, 0.8)+
  theme_bw()+
  theme(plot.title = element_text(size=10), axis.title.x = element_text(size=10), axis.title.y = element_text(size=10))

#plot optimal alpha's and lambda's against nr variables selected
outerCV_nvars_plotdata <- outerCV_BER_nvars %>%
  filter(lambda != 0)

a_vs_nvars <- ggplot(outerCV_nvars_plotdata) +
  geom_point(aes(x = alpha, y = nvars, color = log(lambda), group = log(lambda)))+
  scale_color_lei(palette = "cool", discrete = FALSE)+
  labs(x = expression(alpha), y = "Number of variables selected", color = expression(paste("log (", lambda, ")")),
       title = expression(paste(" # variables selected per optimal ", alpha, " and ", lambda, " combination")))+
  theme_bw()+
  theme(plot.title = element_text(size=10), axis.title.x = element_text(size=10), axis.title.y = element_text(size=10))

#Combine the three plots above: (saved as JPEG, w1448 h315, w1122 h472)
ggarrange(alpha_plot_alll, alpha_plot_opt, a_vs_nvars,
          labels = c("A", "B", "C"),
          ncol = 3, nrow = 1,
          #font.label = list(size=8),
          #align = c("hv"),
          common.legend = TRUE) %>%
            ggexport(filename = "figures/hyperparameter_selection_plots.png", 
                     width = 3000, height = 1200, res = 250)



# Figure: Number of metabolites selected and performance of prediction model----
# * Boxplot visualising the performance (BER)) of the predictive models in 
#    relation to the number of variables selected. 
# * Histogram of the number of variables selected by the elastic net method.

# Plot amount of variables selected vs BER of the first comparison (s.pneu - atyp)
nvars_BER_boxplot <- ggplot(outerCV_BER_nvars)+
  geom_boxplot(aes(x=nvars, y = BER, group = nvars))+
  scale_x_continuous(breaks = c(0, 50, 100, 150, 200, 250, 300, 350))+
  labs(x = "Number of variables selected", y = "Balanced error rate (BER)")+
  theme_bw()

#Histogram nr variables selected
nvars_hist <- ggplot(outerCV_BER_nvars)+
  geom_bar(aes(x = nvars))+
  scale_x_continuous(breaks = c(0, 50, 100, 150, 200, 250, 300, 350))+
  labs(x = "Number of variables selected")+
  theme_bw()

# Histogram and boxplot combined (saved as JPEG, width 636, height 527 / w800 h525)
ggarrange(nvars_BER_boxplot, nvars_hist,
          labels = c("A", "B"),
          ncol = 1, nrow = 2) %>%
  ggexport(filename = "figures/boxplot_BER_nvars.png", 
           width = 1350, height = 1350, res = 250)




