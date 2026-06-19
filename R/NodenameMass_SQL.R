#' @title extract compounds based on their nodename, compound name, and 
#'       precursorMz in a given experiment.
#' 
#' @description The function NodenameMass_SQL() look for compounds in a given 
#'  sqlite database. It lets the user query compounds in three possible ways: 
#'  
#'    1. by nodename.
#'    2. by compound name (partial match)
#'    3. by m/z value within a tolerance
#'      
#' @param sql_path 'character(1)' path to the sqlite database.
#'
#' @param expid 'numeric(1)' experiment number.
#' 
#' @param nodename 'character(1)' corresponds to the feature_id or nodename in 
#'  the format of "MxTy".
#'  
#' @param name 'character(1)'  name or a part of the name of the compounds.
#'  
#' @param mass_measured 'numeric(1)' the measured mass of the compounds.
#' 
#' @param mass_range the window error of the compound mass.
#' 
#' @param rt_range the window error of the compound retention time.
#'  
#' @return It returns the rows of the compound table that match the query.
#' 
#' @author Ahlam Mentag
#' 
#' @export
NodenameMass_SQL <- function(sql_path, expid = NULL,
                             nodename = NULL,
                             name = NULL,
                             mass_measured = NULL,
                             rt = NULL,
                             mass_range = 0.02,
                             rt_range = 1) {
  
  con <- DBI::dbConnect(RSQLite::SQLite(), sql_path)
  on.exit(DBI::dbDisconnect(con), add = TRUE)
  
  base_query <- "SELECT * FROM ms_compound"
  
  conditions <- c()
  params <- list()
  
  # Optional expid filter
  if (!is.null(expid)) {
    conditions <- c(conditions, "expid = ?")
    params <- c(params, expid)
  }
  
  # Priority logic (same as your original structure)
  if (!is.null(nodename)) {
    conditions <- c(conditions, "nodename = ?")
    params <- c(params, nodename)
    
  } else if (!is.null(name)) {
    conditions <- c(conditions, "name LIKE ?")
    params <- c(params, paste0("%", name, "%"))
    
  } else if (!is.null(mass_measured)) {
    conditions <- c(conditions, "mass_measured BETWEEN ? AND ?")
    params <- c(params, mass_measured - mass_range,
                mass_measured + mass_range)
    
  } else if (!is.null(rt)) {
    cat("the retention time should be in minutes\n")
    
    conditions <- c(conditions, "retention_time BETWEEN ? AND ?")
    params <- c(params, rt - rt_range, rt + rt_range)
  }
  
  # Build final query
  if (length(conditions) > 0) {
    query <- paste(base_query, "WHERE", paste(conditions, collapse = " AND "))
  } else {
    query <- base_query
  }
  
  DBI::dbGetQuery(con, query, params = params)
}