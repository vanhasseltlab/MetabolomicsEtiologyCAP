MakeFolds <- function(seed, nfolds, dat, outcome, sampleID){
  
  # Set seed
  set.seed(seed)
  
  # Determine fold size
  fld_size_G1 <- round(sum(dat[, outcome] == levels(dat[, outcome])[1])*(1/nfolds))
  fld_size_G2 <- round(sum(dat[, outcome] == levels(dat[, outcome])[2])*(1/nfolds))
  
  # Make 5 stratified folds and save in folds list
  flds <- list()
  
  # Fold 1 - 4
  rest <- dat
  for (n in 1:(nfolds-1)){
    # Sample for testset by MNC name per group
    testfold_G1 <- sample(rest[, sampleID][rest[, outcome] == levels(dat[, outcome])[1]], 
                         size = fld_size_G1)
    
    testfold_G2 <- sample(rest[, sampleID][rest[, outcome] == levels(dat[, outcome])[2]], 
                          size = fld_size_G2)
    
    # Save in folds list
    flds[[n]] <- c(testfold_G1, testfold_G2)
    # gather remaining samples
    rest <- filter(rest, !(rest[, sampleID] %in% testfold_G1 | 
                           rest[, sampleID] %in% testfold_G2))
    
    # Name fold
    names(flds)[[n]] <- paste("Fold", n, sep = "")
  }
  
  # last fold
  flds[[nfolds]] <- as.integer(rest[, sampleID])
  # Name last fold
  names(flds)[[nfolds]] <- paste("Fold", nfolds, sep = "")
 
  return(flds) 
}
