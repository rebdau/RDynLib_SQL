allCSPP_SQL <- function(sql_path,
                        exp.id,
                        nr_col = 35,
                        thr1 = 3,
                        thr2 = 0.1,
                        thr3 = 0.1) {
  
  con <- dbConnect(RSQLite::SQLite(), sql_path)
  on.exit(dbDisconnect(con), add = TRUE)
  
  ##Get compound IDs for experiment 
  compounds <- dbGetQuery(
    con,
    sprintf("SELECT compound_id FROM ms_compound WHERE expid = %d", exp.id)
  )
  
  if (nrow(compounds) == 0)
    return(data.frame())
  
  comp_ids <- compounds$compound_id
  
  ## Load compound_add table 
  compound_add <- dbGetQuery(con, "SELECT * FROM compound_add")
  
  ## Keep only rows for experiment compounds
  compound_add <- compound_add[compound_add[,1] %in% comp_ids, ]
  
  cspp.res <- data.frame(
    compid.sub = integer(),
    compid.prod = integer(),
    conv.type = character(),
    ions.prod = integer(),
    ave.common = numeric(),
    ave.dot = numeric(),
    stringsAsFactors = FALSE
  )
  
  k <- 1
  
  for (i in seq_len(nrow(compound_add))) {
    
    compid.sub <- as.integer(compound_add[i,1])
    
    for (j in 3:nr_col) {
      
      cell_value <- compound_add[i,j]
      if (is.na(cell_value) || cell_value == "") next
      
      parts <- strsplit(cell_value, "!!")[[1]]
      if (length(parts) < 3) next
      
      compid.prod <- as.integer(parts[3])
      if (is.na(compid.prod)) next
      
      subparts <- strsplit(parts[2], "!")[[1]]
      if (length(subparts) < 3) next
      
      cspp.res[k,] <- list(
        compid.sub,
        compid.prod,
        colnames(compound_add)[j],
        as.integer(subparts[1]),
        as.numeric(subparts[2]),
        as.numeric(subparts[3])
      )
      
      k <- k + 1
    }
  }
  
  ##Apply thresholds 
  cspp.res <- cspp.res[
    cspp.res$ions.prod > thr1 &
      cspp.res$ave.common > thr2 &
      cspp.res$ave.dot > thr3, ]
  
  cspp.res <- cspp.res[order(cspp.res$compid.sub), ]
  
  return(cspp.res)
}
