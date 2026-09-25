#' @title extracting the msms_spectrum table of a given compound_id(s)  
#'
#' @param qt_path the file path to the sql database
#' @param compound_id a given or a list of compound_ids
#' 
#' @return the msms_spectrum table of the given compound_id(s) 
#' 
#' @import DBI
#' @import RSQLite
#' @import CompoundDb
#' 
#' @author Ahlam Mentag
#' 
#' @export
get_spectrum_id <- function(qt_path, compound_id){
  con_qt <- dbConnect(SQLite(), qt_path)
  
  compound_id <- as.numeric(unlist(compound_id))
  
  # Create ID list for SQL
  ids <- paste(compound_id, collapse = ", ")
  
  query <- dbGetQuery(con_qt, sprintf("Select * from msms_spectrum 
                                            where compound_id IN (%s)", 
                                            ids))
  return(query)
}




