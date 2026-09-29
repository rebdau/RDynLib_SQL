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
sql_to_spectra <- function(sql_path) {
  
  db <- DBI::dbConnect(
    RSQLite::SQLite(),
    sql_path
  )
  
  ## Get available tables
  tables <- DBI::dbListTables(db)
  
  
  ## Clean experiment table if necessary
  
  if ("experiment" %in% tables) {
    
    experiment_df <- DBI::dbGetQuery(
      db,
      "SELECT * FROM experiment"
    )
    
    cols_to_remove <- intersect(
      c("x", "x.1"),
      names(experiment_df)
    )
    
    if (length(cols_to_remove) > 0) {
      
      experiment_df <- experiment_df[
        ,
        !names(experiment_df) %in% cols_to_remove,
        drop = FALSE
      ]
      
      DBI::dbExecute(
        db,
        "DROP TABLE experiment"
      )
      
      DBI::dbWriteTable(
        db,
        "experiment",
        experiment_df
      )
    }
  }
  
  
  ## Close connection
  DBI::dbDisconnect(db)
  
  
  ## Convert SQLite -> ExpCompDb
  
  cdb <- ExpCompDb(sql_path)
  
  
  ## Add DynLib-specific joins
  if ("compound_add" %in% tables) {
    
    cdb <- CompoundDb::addJoinDefinition(
      cdb,
      table_a = "ms_compound",
      table_b = "compound_add",
      column_a = "compound_id",
      column_b = "compound_id"
    )
  }
  
  
  if ("gnps_add" %in% tables) {
    
    cdb <- CompoundDb::addJoinDefinition(
      cdb,
      table_a = "ms_compound",
      table_b = "gnps_add",
      column_a = "compound_id",
      column_b = "compound_id"
    )
  }
  
  
  ## ExpCompDb -> Spectra
  dy_sps <- Spectra::Spectra(cdb)
  
  return(dy_sps)
}