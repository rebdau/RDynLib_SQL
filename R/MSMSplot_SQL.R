MSMSplot_SQL <- function(sql_path,
                         dbkey,
                         err = NULL,
                         minum = NULL,
                         oldpar = NULL,
                         nl = NULL) {
  
  library(DBI)
  library(RSQLite)
  
  if (is.null(err)) err <- 0.015
  if (is.null(minum)) minum <- 2
  
  con <- dbConnect(SQLite(), sql_path)
  on.exit(dbDisconnect(con), add = TRUE)
  
  # Get precursor m/z 
  compound <- dbGetQuery(con, sprintf(
    "SELECT compound_id, precursor_mz
     FROM ms_compound
     WHERE compound_id = %d",
    dbkey
  ))
  
  if (nrow(compound) == 0)
    stop("Compound not found")
  
  precursor_mz <- compound$precursor_mz
  
  #Get product ions and intensities 
  peaks <- dbGetQuery(con, sprintf(
    "SELECT p.Mz, p.Intensity
     FROM msms_spectrum s
     JOIN msms_spectrum_peak p ON s.spectrum_id = p.spectrum_id
     WHERE s.compound_id = %d
     ORDER BY p.Mz",
    dbkey
  ))
  
  if (nrow(peaks) == 0)
    stop("No MS/MS peaks found for this compound")
  
  prod_ion.num <- round(peaks$Mz, 2)
  intens_ion.num <- peaks$Intensity
  
  # Candidate product ions 
  cat("\nCandidate product ions:\n")
  
  # You must adapt these to SQLite if needed
  # For now assumed available
  ProdIonMatch(prod_ion.num, err = err)
  
  minum_sel <- which(intens_ion.num >= minum)
  
  # Neutral losses 
  neutprod.num <- round(precursor_mz - prod_ion.num, 2)
  
  cat("\nCandidate neutral losses:\n")
  NeutLossMatch(neutprod.num, err = err)
  
  # Complementary ions
  cat("\nComplementary product ions:\n")
  ComplemIons(prod_ion.num,
              intens_ion.num,
              neutprod.num,
              err)
  
  #  Plot 
  if (is.null(oldpar)) oldpar <- par(no.readonly = TRUE)
  
  par(mfrow = c(1, 1))
  par(cex = 0.6)
  
  adj_intens <- max(intens_ion.num) * 1.1
  adj_prod_min <- min(prod_ion.num) * 0.9
  adj_prod_max <- max(prod_ion.num) * 1.1
  
  plot(prod_ion.num,
       intens_ion.num,
       type = "h",
       xlab = "m/z",
       ylab = "ion intensity",
       main = round(precursor_mz, 2),
       xlim = c(adj_prod_min, adj_prod_max),
       ylim = c(0, adj_intens))
  
  text(prod_ion.num[minum_sel],
       intens_ion.num[minum_sel],
       prod_ion.num[minum_sel],
       pos = 3)
  
  #Interactive neutral loss mode 
  if (!is.null(nl)) {
    
    text(prod_ion.num[minum_sel],
         intens_ion.num[minum_sel],
         neutprod.num[minum_sel],
         pos = 3,
         offset = 1.5,
         col = 2)
    
    cat("\nEnter 0 to exit.\n")
    i <- 2.5
    
    repeat {
      chosen_ion <- as.numeric(
        readline("From which product ion do you want to see the neutral losses? ")
      )
      
      if (chosen_ion == 0) break
      
      sl <- which(prod_ion.num %in% chosen_ion)
      if (length(sl) == 0) next
      
      neutprod.num <- round(prod_ion.num[sl] - prod_ion.num, 2)
      
      NeutLossMatch(neutprod.num, err = err)
      
      text(prod_ion.num[minum_sel],
           intens_ion.num[minum_sel],
           neutprod.num[minum_sel],
           pos = 3,
           offset = i,
           col = 3)
      
      i <- i + 1
    }
  }
  
  if (!is.null(oldpar)) par(oldpar)
}
