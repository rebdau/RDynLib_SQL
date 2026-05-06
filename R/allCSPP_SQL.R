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
                        thr1 = 0,
                        thr2 = 0,
                        thr3 = 0) {
  
  library(data.table)
  library(DBI)
  library(RSQLite)
  
  con <- dbConnect(RSQLite::SQLite(), sql_path)
  on.exit(dbDisconnect(con), add = TRUE)
  
  # MS compounds filter
  compounds <- dbGetQuery(
    con,
    sprintf("SELECT compound_id FROM ms_compound WHERE expid = %d", exp.id)
  )
  
  if (nrow(compounds) == 0) return(data.frame())
  
  comp_ids <- as.character(compounds$compound_id)
  
  compound_add <- as.data.table(dbGetQuery(con, "SELECT * FROM compound_add"))
  compound_add <- compound_add[compound_id %in% comp_ids]
  
  cols <- setdiff(names(compound_add), "compound_id")
  
  res_list <- list()
  k <- 1
  

  for (i in seq_len(nrow(compound_add))) {
    
    sub_id <- compound_add$compound_id[i]
    
    for (col in cols) {
      
      cell <- compound_add[[col]][i]
      if (is.na(cell) || cell == "") next
      
      # IMPORTANT: multiple entries are separated by |
      entries <- unlist(strsplit(cell, "\\|"))
      
      for (entry in entries) {
        
        entry <- trimws(entry)
        if (entry == "") next
        
        # expected format:
        # !!COMMON!FORW_REV!DOT!!PRODID
        parts <- unlist(strsplit(entry, "!!", fixed = TRUE))
        if (length(parts) < 3) next
        
        metrics <- unlist(strsplit(parts[2], "!", fixed = TRUE))
        if (length(metrics) < 3) next
        
        prod_id <- parts[3]
        
        res_list[[k]] <- data.table(
          compid.sub  = as.integer(sub_id),
          compid.prod = as.integer(prod_id),
          conv.type   = col,
          ions.prod   = as.numeric(metrics[1]),
          ave.common  = as.numeric(metrics[2]),
          ave.dot     = as.numeric(metrics[3])
        )
        
        k <- k + 1
      }
    }
  }
  
  if (length(res_list) == 0) return(data.frame())
  
  res <- rbindlist(res_list, fill = TRUE)
  

  res <- res[
    ions.prod >= thr1 &
      ave.common >= thr2 &
      ave.dot >= thr3
  ]
  
  setorder(res, compid.sub)
  
  return(res)
}