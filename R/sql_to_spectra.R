#' @title convert an sql database to a spectra object 
#'
#' @param sql_path the path to the sql database
#' 
#' @return a spectra object of the sql database
#' 
#' @import DBI
#' @import RSQLite
#' @import CompoundDb
#' 
#' @author Ahlam Mentag
#' 
#' @export
sql_to_spectra <- function(sql_path){
  
  db <- dbConnect(SQLite(), sql_path)
  dbListTables(db)
  #dbExecute(db, "DROP TABLE IF EXISTS feature_matrix ")
  experiment_df <- dbGetQuery(db, "SELECT * FROM experiment")
  
  #here we remove "x", "x.1" columns from experiment table 
  experiment_df <- experiment_df[, !names(experiment_df) %in% c("x", "x.1")]
  dbExecute(db, "DROP TABLE experiment")
  dbWriteTable(db, "experiment", experiment_df)
  
  #the columns have been successfully deleted
  dbGetQuery(db, "PRAGMA table_info(experiment)")
  dbDisconnect(db)
  #we convert the sql database to a Compdb object
  cdb <- ExpCompDb(sql_path)
  #we convert the Compdb object to a spectra object
  dy_sps <- Spectra(cdb)
}