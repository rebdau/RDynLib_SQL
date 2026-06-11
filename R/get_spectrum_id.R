#' @title extracting the msms_spectrum table of a given compound_id  
#'
#' 
#' @return the msms_spectrum table of a given compound_id 
#' 
#' @import DBI
#' @import RSQLite
#' @import CompoundDb
#' 
#' @author Ahlam Mentag
#' 
#' @export
get_spectrum_id <- function(compound_id){
  compound_id <- dbGetQuery(con_qt, sprintf("Select * from msms_spectrum 
                                            where compound_id = %s", 
                                            compound_id))
}