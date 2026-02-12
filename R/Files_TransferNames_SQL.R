#' @tiltle Load compound metadata from a SQLite database and association file
#'
#' @description Files_TransferNames_SQL() function
#' reads compound metadata from a SQLite database and prepares
#' subdatabase specific compound tables based on the database names defined
#' in a compound association file. 
#'   
#' @param sqlite_path Character string. Path to the SQLite database containing
#'   the \code{ms_compound} table.
#' @param assoc_path Character string. Path to the association file defining
#'   reference and target databases and compound ID mappings.
#'
#' @return A list with four elements:
#' \describe{
#'   \item{ftng}{Data frame of compounds from the reference database.}
#'   \item{ftps}{Data frame of compounds from the target database.}
#'   \item{syn}{Data frame containing the union of reference and target compounds.}
#'   \item{Assoc}{The association data frame as read from \code{assoc_path}.}
#' }
#'
#' @author Ahlam Mentag
#' 
#' @export
Files_TransferNames_SQL <- function(sqlite_path, assoc_path) {
  
  library(DBI)
  library(RSQLite)
  

  #Read association file
  Assoc <- read.table(
    assoc_path,
    header = TRUE,
    sep = "\t",
    stringsAsFactors = FALSE
  )
  
  ref_db    <- unique(Assoc$ref_database)
  target_db <- unique(Assoc$target_database)
  
  if (length(ref_db) != 1 || length(target_db) != 1) {
    stop("Association file must contain exactly one ref and one target database")
  }
  
  ref_db    <- ref_db[1]
  target_db <- target_db[1]
  

  #Connect to SQLite
  con <- dbConnect(SQLite(), sqlite_path)
  on.exit(dbDisconnect(con), add = TRUE)
  

  #Helper to read compounds
  read_compounds <- function(db_name) {
    
    df <- dbGetQuery(
      con,
      paste0(
        "SELECT
           compound_id ,
           name,
           subsid,
           smiles
         FROM ms_compound
         WHERE source_database = '", db_name, "'"
      )
    )
    
    df[] <- lapply(df, as.character)
    df$name[is.na(df$name)] <- "NULL"
    
    df
  }
  

  ## Load sub-databases
  ftng <- read_compounds(ref_db)
  ftps <- read_compounds(target_db)
  
  syn <- unique(rbind(ftng, ftps))
  
  return(list(ftng, ftps, syn, Assoc))
}
