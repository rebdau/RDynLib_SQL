#' @title extracting the msms_spectrum table of a given compound_id(s)  
#'
#' @param sql_path the file path to the sql database
#' @param compound_id a given or a list of compound_ids
#' 
#' @return the msms_spectrum table of the given compound_id(s) 
#' 
#' @import DBI
#' @import RSQLite
#' @import CompoundDb
#' 
#' @author Ahlam Mentag, Rebecca Dauwe
#' 
#' @export
get_spectrum_id <- function(sql_path, compound_id){
  con <- dbConnect(SQLite(), sql_path)
  
  compound_id <- as.numeric(unlist(compound_id))
  
  # Create ID list for SQL
  ids <- paste(compound_id, collapse = ", ")
  
  query <- dbGetQuery(con, sprintf("Select * from msms_spectrum 
                                            where compound_id IN (%s)", 
                                            ids))
  return(query)
}




