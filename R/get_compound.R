#' @title extracting the ms_compound table of a given compound_id  
#'
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
get_compound <- function(qt_path, compound_id){
  con_qt <- dbConnect(SQLite(), qt_path)
  compound_id <- dbGetQuery(con_qt, sprintf("Select * from ms_compound 
                                            where compound_id = %s", 
                                            compound_id))
}