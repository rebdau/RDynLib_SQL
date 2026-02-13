#' @tiltle Load compound metadata from a SQLite database and association file
#'
#' @description Files_TransferNames_SQL() function
#' reads compound metadata from a SQLite database and prepares
#' subdatabase specific compound tables based on the database names defined
#' in a compound association file. 
#'   
#' @param sqlite_dir Character string. Directory containing SQLite databases.
#'   Each database file must be named <database_name>.sqlite.
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
Files_TransferNames_SQL <- function(sqlite_dir, assoc_path) {
  
  library(DBI)
  library(RSQLite)
  
  # Read association file
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
  
  # Helper function that connects to the correct database
  read_compounds <- function(db_name) {
    
    sqlite_path <- file.path(sqlite_dir, paste0(db_name, ".sqlite"))
    
    if (!file.exists(sqlite_path)) {
      stop(paste("Database file not found:", sqlite_path))
    }
    
    con <- dbConnect(SQLite(), sqlite_path)
    on.exit(dbDisconnect(con), add = TRUE)
    
    df <- dbGetQuery(
      con,
      "SELECT compound_id, name, subsid, smiles FROM ms_compound"
    )
    
    df[] <- lapply(df, as.character)
    df$name[is.na(df$name)] <- "NULL"
    
    df
  }
  
  # Load both databases
  ftng <- read_compounds(ref_db)
  ftps <- read_compounds(target_db)
  
  syn <- unique(rbind(ftng, ftps))
  
  return(list(
    ftng = ftng,
    ftps = ftps,
    syn  = syn,
    Assoc = Assoc
  ))
}
