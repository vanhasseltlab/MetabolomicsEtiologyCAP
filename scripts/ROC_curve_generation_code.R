## ROC curves of cross validated elastic net models ----------------
## Elastic net without covariates
# Atypical versus other
elnet.atyp.other.roc <- ElnetROC(dat, elnet.atyp.other.results, 
                                 repeats = 100, nfolds = 5)
# Spneumoniae versus other
elnet.spneu.other.roc <- ElnetROC(dat, elnet.spneu.other.results, 
                                  repeats = 100, nfolds = 5)
# Viral versus bacterial
elnet.bac.vir.roc <- ElnetROC(dat, elnet.bac.vir.results, 
                              repeats = 100, nfolds = 5)

## Elastic net with covariates age and sex
# Atypical versus other
cov.elnet.atyp.other.roc <- ElnetROC(dat, cov.elnet.atyp.other.results, 
                                     repeats = 100, nfolds = 5)
# Spneumoniae versus other
cov.elnet.spneu.other.roc <- ElnetROC(dat, cov.elnet.spneu.other.results, 
                                      repeats = 100, nfolds = 5)
# Viral versus bacterial
cov.elnet.bac.vir.roc <- ElnetROC(dat, cov.elnet.bac.vir.results, 
                                  repeats = 100, nfolds = 5)

## Elastic net with all covariates
# Atypical versus other
cov.all.elnet.atyp.other.roc <- ElnetROC(dat, cov.all.elnet.atyp.other.results, 
                                     repeats = 100, nfolds = 5)
# Spneumoniae versus other
cov.all.elnet.spneu.other.roc <- ElnetROC(dat, cov.all.elnet.spneu.other.results, 
                                      repeats = 100, nfolds = 5)
# Viral versus bacterial
cov.all.elnet.bac.vir.roc <- ElnetROC(dat, cov.all.elnet.bac.vir.results, 
                                  repeats = 100, nfolds = 5)


## Visualization of the results ----------------------------------------------

## Plot ROC curves of all models together ----
# for the comparison atyp-other 
roc.plotdat.atyp.other <- bind_rows(data.frame(logit.glycylglycine[[3]], 
                                               model = "LR: Glycine"), 
                                    data.frame(logit.sdma[[3]], 
                                               model = "LR: SDMA"),
                                    data.frame(logit.lpi18.1[[3]], 
                                               model = "LR: LPI(18:1)"),
                                    data.frame(logit.sigatyp[[3]], 
                                               model = "LR: Glycine, SDMA & LPI(18:1)"),
                                    data.frame(logit.sigatyp.age.sex[[3]], 
                                               model = "LR: Glycine, SDMA, LPI(18:1), age & sex"),
                                    data.frame(logit.sigatyp.cov[[3]], 
                                               model = "LR: Glycine, SDMA, LPI(18:1) & all covariates"),
                                    data.frame(elnet.atyp.other.roc[[2]],
                                               model = "EN: All metabolites"),
                                    data.frame(cov.elnet.atyp.other.roc[[2]],
                                               model = "EN: All metabolites, age & sex"), 
                                    data.frame(cov.all.elnet.atyp.other.roc[[2]],
                                               model = "EN: All metabolites & all covariates"))

roc.plotdat.atyp.other$model <- factor(roc.plotdat.atyp.other$model, 
                                       levels = c("LR: Glycine", "LR: SDMA", "LR: LPI(18:1)",
                                                  "LR: Glycine, SDMA & LPI(18:1)", 
                                                  "LR: Glycine, SDMA, LPI(18:1), age & sex",
                                                  "LR: Glycine, SDMA, LPI(18:1) & all covariates",
                                                  "EN: All metabolites",
                                                  "EN: All metabolites, age & sex", 
                                                  "EN: All metabolites & all covariates"))

roc.plot.atyp.other <- ggplot()+
  geom_step(data = roc.plotdat.atyp.other, aes(x = fpr.mean, y = tpr.mean, 
                                               color = model, linetype = model))+
  geom_segment(aes(x = 0, y = 0, xend = 1, yend = 1), 
               linetype = "dotted", color = "gray")+
  scale_color_lei(palette = "nine")+
  scale_linetype_manual(values = c(rep(2, 3), rep(5, 3), rep(1, 3)))+
  labs(x = "False Positive Rate", y = "True Positive Rate", 
       colour = "Model", linetype = "Model",
       title = "ROC curves: atypical pathogen models")+ 
  theme(legend.title = element_text(size = 18),
      legend.text = element_text(size = 18))+
  theme_bw() 

ggexport(roc.plot.atyp.other, filename = "figures/ROC_curves_atyp_other_combined.png",
         width = 1500, height = 800, res = 250)

## for comparison s.pneu - other
roc.plotdat.spneu.other <- bind_rows(data.frame(elnet.spneu.other.roc[[2]],
                                                model = "EN: All metabolites"),
                                     data.frame(cov.elnet.spneu.other.roc[[2]],
                                                model = "EN: All metabolites, age & sex"), 
                                     data.frame(cov.all.elnet.spneu.other.roc[[2]], 
                                                model = "EN: All metabolites & all covariates"))
roc.plotdat.spneu.other$model <- factor(roc.plotdat.spneu.other$model, 
                                        levels = c("EN: All metabolites",
                                                   "EN: All metabolites, age & sex", 
                                                   "EN: All metabolites & all covariates"))
roc.plot.spneu.other <- ggplot()+
  geom_step(data = roc.plotdat.spneu.other, aes(x = fpr.mean, y = tpr.mean, 
                                                color = model))+
  geom_segment(aes(x = 0, y = 0, xend = 1, yend = 1), 
               linetype = "dotted", color = "gray")+
  scale_color_lei(palette = "three")+
  scale_linetype_manual(values = c(rep(2, 5), 1, 1))+
  labs(x = "False Positive Rate", y = "True Positive Rate", 
       colour = "Model", linetype = "Model",
       title = "ROC curves: S.pneumoniae models")+
  theme_bw() 

ggexport(roc.plot.spneu.other, filename = "figures/ROC_curves_spneu_other_combined.png",
         width = 1500, height = 800, res = 250)

## for comparison viral - other
roc.plotdat.bac.vir <- bind_rows(data.frame(elnet.bac.vir.roc[[2]],
                                            model = "EN: All metabolites"),
                                 data.frame(cov.elnet.bac.vir.roc[[2]],
                                            model = "EN: All metabolites, age & sex"), 
                                 data.frame(cov.all.elnet.bac.vir.roc[[2]], 
                                            model = "EN: All metabolites & all covariates"))
roc.plotdat.bac.vir$model <- factor(roc.plotdat.bac.vir$model, 
                                    levels = c("EN: All metabolites",
                                               "EN: All metabolites, age & sex", 
                                               "EN: All metabolites & all covariates"))
roc.plot.bac.vir <- ggplot()+
  geom_step(data = roc.plotdat.bac.vir, aes(x = fpr.mean, y = tpr.mean, 
                                            color = model))+
  geom_segment(aes(x = 0, y = 0, xend = 1, yend = 1), 
               linetype = "dotted", color = "gray")+
  scale_color_lei(palette = "three")+
  scale_linetype_manual(values = c(rep(2, 5), 1, 1))+
  labs(x = "False Positive Rate", y = "True Positive Rate", 
       colour = "Model", linetype = "Model",
       title = "ROC curves: viral pathogen models")+
  theme_bw() 

ggexport(roc.plot.bac.vir, filename = "figures/ROC_curves_vir_other_combined.png",
         width = 1500, height = 800, res = 250)

# Plot all three comparions in same window
ggarrange(plotlist = list(roc.plot.atyp.other, 
                          roc.plot.spneu.other, 
                          roc.plot.bac.vir), 
          ncol = 3, 
          align = "h", 
          font.label = list(size = 17),
          common.legend = TRUE,
          legend = "bottom") %>%
  ggexport(filename = "figures/ROC_curves_all_combined.png", 
           width = 2800, height = 1100, res = 250)

# Plot seperate ROC curves of all models (supplementary information)
ggarrange(logit.glycylglycine[[2]], logit.sdma[[2]],
          logit.lpi18.1[[2]], logit.sigatyp[[2]],
          logit.sigatyp.age.sex[[2]], logit.sigatyp.cov[[2]],
          elnet.atyp.other.roc[[1]], cov.elnet.atyp.other.roc[[1]], 
          cov.all.elnet.atyp.other.roc[[1]]) %>%
  ggexport(filename = "figures/ROC_curves_atyp_other_separate.png",
           width = 2500, height = 2250, res = 250)
