Ttesting <- function(metrange, outcome_column, outcome_1, outcome_2, dat){
  
  # This function performs an t-test on multiple metabolites
  
  library(data.table)
  
  # Make dataframe for t-test p-values
  pval <- as.data.frame(matrix(nrow = 2, ncol = 347), 
                        row.names = c("p.value", "q.value"))
  names(pval) <- c(metrange)
  
  # Dataframe for two outcomes
  dat.1 <- dat[outcome_column == outcome_1, ]
  dat.2 <- dat[outcome_column == outcome_2, ]
  
  for (i in 1:length(metrange)){ 
    # define metabolite to test
    met <- metrange[[i]]
    # Make vectors to compare
    vector.1 <- dat.1[, met]
    vector.2 <- dat.2[, met]
    # T.test of vector 1 and 2
    pval[1, i] <- t.test(vector.1, vector.2)$p.value
  }

  # Perform FDR multiple testing correction
  pval[2, ] <- p.adjust(pval[1, ], method = "fdr")
  
  # Transpose dataframe
  pval_transposed <- transpose(pval)
  # get row and colnames in order
  colnames(pval_transposed) <- rownames(pval)
  rownames(pval_transposed) <- colnames(pval)
  
  return(pval_transposed)
}
