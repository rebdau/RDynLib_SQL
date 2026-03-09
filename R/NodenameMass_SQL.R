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
NodenameMass_SQL <- function(sql_path, expid,
                             nodename = NULL,
                             name = NULL,
                             mass_measured = NULL,
                             rt = NULL,
                             mass_range = 0.02,
                             rt_range = 1) {
  
  # Create connection
  con <- DBI::dbConnect(RSQLite::SQLite(), sql_path)
  
  on.exit(DBI::dbDisconnect(con), add = TRUE)
  
  if (!is.null(nodename)) {
    query <- "
    SELECT *
    FROM ms_compound
    WHERE expid = ?
    AND nodename = ?
    "
    
    return(DBI::dbGetQuery(con, query, params = list(expid, nodename)))
  }
  
  if (!is.null(name)) {
    query <- "
    SELECT *
    FROM ms_compound
    WHERE expid = ?
    AND name LIKE ?
    "
    
    return(DBI::dbGetQuery(con, query,
                           params = list(expid, paste0('%', name, '%'))))
  }
  
  if (!is.null(mass_measured)) {
    query <- "
    SELECT *
    FROM ms_compound
    WHERE expid = ?
    AND mass_measured BETWEEN ? AND ?
    "
    
    return(DBI::dbGetQuery(con, query,
                           params = list(expid, mass_measured - mass_range, 
                                         mass_measured + mass_range)))
  }
  
  if (!is.null(rt)) {
    
    cat("the retention time should be in minutes\n")
    
    query <- "
    SELECT *
    FROM ms_compound
    WHERE expid = ?
    AND retention_time BETWEEN ? AND ?
    "
    
    return(DBI::dbGetQuery(con, query,
                           params = list(expid, rt - rt_range, rt + rt_range)))
  }
  
}