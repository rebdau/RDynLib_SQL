#' @title generating compound pairs with their associated similarity metrics.
#' 
#' @description
#' allCSPP_SQL() connects to a SQLite database, retrieves CSPP-based 
#' relationships between compounds for a given experiment (exp.id), parses 
#' pairwise similarity information from the gnps_add table, filters the results 
#' using defined thresholds (thr1, thr2, thr3).
#' 
#' @param sql_path 'character(1)' path to the sqlite database.
#' 
#' @param exp.id 'number(1)' experiment to use for network generation.
#' 
#' @param nr_col number of columns in the compound_add table.
#' 
#' @param thr1 'number(1)' minimum number of product ions (varying between 0 
#' and 1) of one CID spectrum that can be traced in another CID spectrum,
#' set on 1 by default.
#' 
#' @param thr2 'number(1)' dot product threshold for the common ions 
#' between two CID spectra, set on 0.9 by default.
#' 
#' @param thr3 'number(1)' average of the dot products obtained for the common 
#' product ions and the common neutral losses, by default set on 0.4.
#' 
#' @returns a data frame of compound pairs with their associated CSPP similarity
#'          metrics.
#'          
#' @import DBI
#' @import RSQLite
#' 
#' @author Ahlam Mentag
#' 
#' @export

allCSPP_SQL <- function(sql_path,
                        exp.id,
                        nr_col = NULL,
                        thr1 = 1,
                        thr2 = 0.9,
                        thr3 = 0.4) {
  
  con <- dbConnect(RSQLite::SQLite(), sql_path)
  on.exit(dbDisconnect(con), add = TRUE)
  
  if (is.null(nr_col)) {
    table_info <- dbGetQuery(con, "PRAGMA table_info(compound_add)")
    nr_col <- nrow(table_info)
  }
  
  compounds <- dbGetQuery(
    con,
    sprintf("SELECT compound_id FROM ms_compound WHERE expid = %d", exp.id)
  )
  
  if (nrow(compounds) == 0) 
    return(data.frame())
  
  comp_ids <- compounds$compound_id
  
  compound_add <- dbGetQuery(con, "SELECT * FROM compound_add")
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
  max_col <- min(nr_col, ncol(compound_add))
  
  for (i in seq_len(nrow(compound_add))) {
    
    compid.sub <- as.integer(compound_add[i,1])
    
    for (j in 3:max_col) {
      
      cell_value <- compound_add[i,j]
      if (is.na(cell_value) || cell_value == "") next
      
      parts <- strsplit(cell_value, "!!")[[1]]
      if (length(parts) < 3) next
      
      compid.prod <- as.integer(trimws(parts[3]))
      if (is.na(compid.prod)) next
      
      subparts <- strsplit(parts[2], "!")[[1]]
      if (length(subparts) < 3) next
      subparts <- trimws(subparts)
      
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
  
  cspp.res <- cspp.res[
    cspp.res$ions.prod >= thr1 &
      cspp.res$ave.common >= thr2 &
      cspp.res$ave.dot >= thr3, ]
  
  cspp.res <- cspp.res[order(cspp.res$compid.sub), ]
  
  return(cspp.res)
}