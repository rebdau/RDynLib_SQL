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
#' @param precursor_Mz 'numeric(1)' the precursor mz of the compounds.
#' 
#' @param err the window error of the precursorMz.
#'  
#' @return It returns the rows of the compound table that match the query.
#' 
#' @author Ahlam Mentag
#' 
#' @export
NodenameMass_SQL <- function(sql_path, expid,
                             nodename = NULL,
                             name = NULL,
                             precursor_Mz = NULL,
                             err = 0.02) {
  
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
  
  if (!is.null(precursor_Mz)) {
    query <- "
    SELECT *
    FROM ms_compound
    WHERE expid = ?
    AND mass_measured BETWEEN ? AND ?
    "
    
    return(DBI::dbGetQuery(con, query,
                           params = list(expid, precursor_Mz - err, 
                                         precursor_Mz + err)))
  }
  
}