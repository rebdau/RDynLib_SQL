#' @title extracting the existing experiment id(s) in an sql database 
#'
#' @param qt_path the file path to the sql database
#' 
#' @return the existing expid(s) in a given database 
#' 
#' @import DBI
#' @import RSQLite
#' 
#' @author Ahlam Mentag
#' 
#' @export
get_expid <- function(qt_path){
  con_qt <- dbConnect(SQLite(), qt_path)
  
  query <- dbGetQuery(con_qt, sprintf("Select expid from experiment 
                                            "))
  return(query)
}




