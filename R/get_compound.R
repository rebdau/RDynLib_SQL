#' @title extracting the ms_compound table of given compound_id(s) 
#'
#' @param sql_path the file path to the sql database
#' @param compound_id a given or a list of compound_ids
#' 
#' @return the ms_compound table of a given compound_id 
#' 
#' @import DBI
#' @import RSQLite
#' @import CompoundDb
#' 
#' @author Ahlam Mentag, Rebecca Dauwe
#' 
#' @export
get_compound <- function(sql_path, compound_id) {
  
  con <- DBI::dbConnect(RSQLite::SQLite(), sql_path)
  on.exit(DBI::dbDisconnect(con), add = TRUE)
  
  # Convert to numeric vector
  compound_id <- as.numeric(unlist(compound_id))
  
  # Create ID list for SQL
  ids <- paste(compound_id, collapse = ", ")
  
  query <- sprintf(
    "SELECT *
         FROM ms_compound
         WHERE compound_id IN (%s)",
    ids
  )
  
  compounds <- DBI::dbGetQuery(con, query)
  
  return(compounds)
}