#' @title extracting the available_spectrum_types table in an sql database.
#'
#' @param con connection to the database
#' 
#' @return the available spectrum_types 
#' 
#' @import DBI
#' @import RSQLite
#' @import CompoundDb
#' 
#' @author Ahlam Mentag
#' 
#' @export
available_spectrum_types <- function(con) {
  
  dbGetQuery(
    con,
    "SELECT DISTINCT spectrum_type
     FROM msms_spectrum
     ORDER BY spectrum_type"
  )
}