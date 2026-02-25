ComplemIons_SQL <- function(prod_ion.num, intens_ion.num,
                            neutprod.num, err) {
  
  comp.list <- list()
  k <- 1
  
  for (i in seq_along(neutprod.num)) {
    
    complem_ion <- neutprod.num[i] - 1.008
    low  <- complem_ion - err
    high <- complem_ion + err
    
    for (j in seq_along(prod_ion.num)) {
      
      if (prod_ion.num[j] >= low &&
          prod_ion.num[j] <= high) {
        
        comp.list[[k]] <- c(
          prod_ion.num[i],
          intens_ion.num[i],
          prod_ion.num[j],
          intens_ion.num[j]
        )
        k <- k + 1
      }
    }
  }
  
  # Si aucun complémentaire trouvé
  if (length(comp.list) == 0) {
    cat("No complementary ions found\n")
    return(invisible(NULL))
  }
  
  comp.ion <- do.call(rbind, comp.list)
  
  # Retirer les doublons miroir (comme dans l'original)
  retain_rows <- floor(nrow(comp.ion) / 2)
  comp.ion <- comp.ion[1:retain_rows, , drop = FALSE]
  
  colnames(comp.ion) <- c("Ion 1",
                          "Intensity 1",
                          "Ion 2",
                          "Intensity 2")
  
  print(comp.ion)
}
