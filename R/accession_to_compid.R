#' @title extract the the corresponding compound_id(s)
#'       of a given compound accession
#'
#' @param qt_path the file path to the sql database
#' @param compound_accession complete or partial compound accession identifiers
#' 
#' @return the corresponding compound_id(s) 
#' 
#' @import DBI
#' @import RSQLite
#' @import CompoundDb
#' 
#' @author Ahlam Mentag, Rebecca Dauwe
#' 
#' @export
accession_to_compid <- function(qt_path, compound_accession) {
  
  con_qt <- DBI::dbConnect(RSQLite::SQLite(), qt_path)
  on.exit(DBI::dbDisconnect(con_qt), add = TRUE)
  
  # Convert to character vector
  compound_accession <- as.character(unlist(compound_accession))
  
  # Create LIKE conditions for partial matching
  conditions <- paste(
    rep("compound_accession LIKE ?", length(compound_accession)),
    collapse = " OR "
  )
  
  query <- sprintf(
    "SELECT compound_id, compound_accession
     FROM ms_compound
     WHERE %s",
    conditions
  )
  
  # Add wildcards for partial matching
  params <- as.list(paste0("%", compound_accession, "%"))
  
  compounds <- DBI::dbGetQuery(
    con_qt,
    query,
    params = params
  )
  
  return(compounds)
}