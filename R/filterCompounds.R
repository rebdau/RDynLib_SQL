#' @title Filter compounds in a DynLib database 
#' 
#' @description 
#' Retrieve compound information from a DynLib database, by nodename, compound name (partial match), 
#' measured m/z, and retention time.
#' 
#'
#' Multiple filtering criteria can be combined. Only compounds satisfying
#' all supplied criteria are returned. Arguments set to NULL are ignored.
#' If no filtering criteria are supplied, the entire compound table is
#' returned.
#'      
#' @param sql_path 'character(1)' path to the sqlite database.
#'
#' @param expid 'numeric(1)' experiment number.
#' 
#' @param nodename `character(1)` Feature identifier (feature_id).
#' In older experiments, nodenames followed the "MxTy" format.
#'  
#' @param name 'character(1)'  Compound name or a part of the compound name.
#'  
#' @param mass_measured `numeric(1)` Experimentally measured m/z.
#'
#' @param retention_time `numeric(1)` Retention time, in minutes.
#' 
#' @param mass_tolerance `numeric(1)` Absolute m/z tolerance applied
#' on either side of `mass_measured`. Defaults to 0.02.
#'
#' @param rt_tolerance `numeric(1)` Absolute retention time tolerance,
#' in minutes, applied on either side of `retention_time`.
#' Defaults to 1.
#'
#' @return A `data.frame` containing all columns of the `ms_compound`
#' table for compounds satisfying the supplied filtering criteria.
#'#' 
#' @author Ahlam Mentag, Rebecca Dauwe
#' 
#' @export
filterCompounds <- function(sql_path,
                            expid = NULL,
                            nodename = NULL,
                            name = NULL,
                            mass_measured = NULL,
                            retention_time = NULL,
                            mass_tolerance = 0.02,
                            rt_tolerance = 1) {
  
  con <- DBI::dbConnect(RSQLite::SQLite(), sql_path)
  on.exit(DBI::dbDisconnect(con), add = TRUE)
  
  base_query <- "SELECT * FROM ms_compound"
  
  conditions <- c()
  params <- list()
  
  # Filter by experiment identifier
  if (!is.null(expid)) {
    conditions <- c(conditions, "expid = ?")
    params <- c(params, expid)
  }
  
  # Combine all supplied filtering criteria
  if (!is.null(nodename)) {
    conditions <- c(conditions, "nodename = ?")
    params <- c(params, nodename)
  }
  
  if (!is.null(name)) {
    conditions <- c(conditions, "name LIKE ?")
    params <- c(params, paste0("%", name, "%"))
  }
  
  if (!is.null(mass_measured)) {
    conditions <- c(conditions, "mass_measured BETWEEN ? AND ?")
    params <- c(params, mass_measured - mass_tolerance,
                mass_measured + mass_tolerance) 
  }
  
  if (!is.null(retention_time)) {  
    conditions <- c(conditions, "retention_time BETWEEN ? AND ?")
    params <- c(params, retention_time - rt_tolerance,
                retention_time + rt_tolerance)
  }
  
  # Build final query
  if (length(conditions) > 0) {
    query <- paste(base_query, "WHERE", paste(conditions, collapse = " AND "))
    DBI::dbGetQuery(con, query, params = params)
  } else {
    DBI::dbGetQuery(con, base_query)  
  }
}