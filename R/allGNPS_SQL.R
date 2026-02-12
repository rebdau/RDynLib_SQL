allGNPS_SQL <- function(sql_path,
                        exp.id,
                        nr_col = 6,
                        thr1 = 3,
                        thr2 = 0.1,
                        thr3 = 0.1) {
  
  con <- dbConnect(RSQLite::SQLite(), sql_path)
  on.exit(dbDisconnect(con), add = TRUE)
  
  compounds <- dbGetQuery(
    con,
    sprintf("SELECT compound_id FROM ms_compound WHERE expid = %d", exp.id)
  )
  
  if (nrow(compounds) == 0)
    return(data.frame())
  
  comp_ids <- compounds$compound_id
  
  gnps_add <- dbGetQuery(con, "SELECT * FROM gnps_add")
  gnps_add <- gnps_add[gnps_add[,1] %in% comp_ids, ]
  
  gnps.res <- data.frame(
    compid.sub = integer(),
    compid.prod = integer(),
    mass.diff = numeric(),
    ions.prod = integer(),
    ave.common = numeric(),
    ave.dot = numeric(),
    stringsAsFactors = FALSE
  )
  
  k <- 1
  
  for (i in seq_len(nrow(gnps_add))) {
    
    compid.sub <- as.integer(gnps_add[i,1])
    
    for (j in 2:nr_col) {
      
      cell_value <- gnps_add[i,j]
      if (is.na(cell_value) || cell_value == "") next
      
      parts <- strsplit(cell_value, "!!")[[1]]
      if (length(parts) < 3) next
      
      subparts <- strsplit(parts[2], "!")[[1]]
      if (length(subparts) < 4) next
      
      compid.prod <- as.integer(subparts[4])
      if (is.na(compid.prod)) next
      
      gnps.res[k,] <- list(
        compid.sub,
        compid.prod,
        as.numeric(parts[3]),
        as.integer(subparts[1]),
        as.numeric(subparts[2]),
        as.numeric(subparts[3])
      )
      
      k <- k + 1
    }
  }
  
  ##Apply thresholds 
  gnps.res <- gnps.res[
    gnps.res$ions.prod > thr1 &
      gnps.res$ave.common > thr2 &
      gnps.res$ave.dot > thr3, ]
  
  gnps.res <- gnps.res[order(gnps.res$compid.sub), ]
  
  return(gnps.res)
}
