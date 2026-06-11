#' @title Adding feature matrix table to the sqlite database
#'
#' @details
#' Adding feature matrix table which contains the feature intensities in each 
#' sample to the sqlite database
#' 
#' @param dbfile `character(1)` the path to the SQL database.
#'
#' @param xcms_exp `character(1)` the path to the XCMSExperiment object.
#'
#' @return Updated SQL database with feature matrix table
#'
#' @importFrom DBI dbGetQuery
#'
#' @importFrom DBI dbConnect
#'
#' @importFrom DBI dbDisconnect
#'
#' @importFrom RSQLite SQLite
#'
#' @import xcms
#'
#' @import tidyr
#'
#' @import tibble
#' 
#' @import dplyr
#' 
#' @author Ahlam Mentag
#'
#' @export
add_FeatureMatrix_SQLite <- function(dbfile, xcms_exp) {
  
  con <- dbConnect(SQLite(), dbfile)
  on.exit(dbDisconnect(con))
  

  ##  Ensure table exists
  dbExecute(con, "
    CREATE TABLE IF NOT EXISTS feature_matrix (
      nodename TEXT,
      dataOrigin TEXT,
      Value REAL,
      expid INTEGER,
      FOREIGN KEY (expid) REFERENCES experiment(expid)
    );
  ")
  
  ##Extract feature values
  feature_values <- xcms::featureValues(xcms_exp)
  

  ##Convert to long format
  long_format <- feature_values %>%
    as.data.frame() %>%
    tibble::rownames_to_column("feature_id") %>%
    tidyr::pivot_longer(
      cols = -feature_id,
      names_to = "dataOrigin",
      values_to = "Value"
    )
  
  ##  Get expid safely
  exp_id <- dbGetQuery(con, "
    SELECT expid
    FROM experiment
    ORDER BY expid DESC
    LIMIT 1
  ")$expid
  
  if (length(exp_id) == 0 || is.na(exp_id)) {
    stop("No expid found in experiment table")
  }
  
  long_format$expid <- exp_id
  
  ## Rename columns
  long_format <- long_format %>%
    dplyr::rename(nodename = feature_id)
  
  ## Insert into SQLite
  dbWriteTable(
    con,
    "feature_matrix",
    long_format,
    append = TRUE,
    row.names = FALSE
  )
  
  invisible(exp_id)
}