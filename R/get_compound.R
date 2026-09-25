#' @title extracting the ms_compound table of given compound_id(s) 
#'
#' @param qt_path the file path to the sql database
#' @param compound_id a given or a list of compound_ids
#' 
#' @return the ms_compound table of a given compound_id 
#' 
#' @import DBI
#' @import RSQLite
#' @import CompoundDb
#' 
#' @author Ahlam Mentag
#' 
#' @export
get_compound <- function(qt_path, compound_id) {
  
  con_qt <- DBI::dbConnect(RSQLite::SQLite(), qt_path)
  on.exit(DBI::dbDisconnect(con_qt), add = TRUE)
  
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
  
  compounds <- DBI::dbGetQuery(con_qt, query)
  
  return(compounds)
}