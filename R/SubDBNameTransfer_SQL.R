#' @title Transfer compound names across SQLite sub-databases.
#' 
#' @description SubDBNameTransfer_SQL() propagates compound names,
#' identifiers, and SMILES strings between associated compounds stored
#' in separate SQLite databases (one per sub-database).
#'
#' @param sqlite_dir Character string. Directory containing SQLite databases.
#'   Each database file must be named <database_name>.sqlite.
#'   
#' @param assoc_path Character string. Path to the association file defining
#'   reference and target compound_ids.
#'
#' @return A data frame of compounds having multiple mappings.
#'
#' @author Ahlam Mentag
#' 
#' @export
SubDBNameTransfer_SQL <- function(sqlite_dir, assoc_path) {
  
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
  
  # Load compound metadata from both databases
  Names_lst <- Files_TransferNames_SQL(
    sqlite_dir = sqlite_dir,
    assoc_path = assoc_path
  )
  
  # Perform name transfer logic
  transfertNam <- TransferNames_SQL(Names_lst)
  
  Updated_Names <- transfertNam[[1]]
  Mult.frame    <- transfertNam[[2]]
  
  # Helper function to update one database
  update_database <- function(db_name, df) {
    
    sqlite_path <- file.path(sqlite_dir, paste0(db_name, ".sqlite"))
    
    if (!file.exists(sqlite_path)) {
      stop(paste("Database file not found:", sqlite_path))
    }
    
    con <- dbConnect(SQLite(), sqlite_path)
    on.exit(dbDisconnect(con), add = TRUE)
    
    for (i in seq_len(nrow(df))) {
      
      dbExecute(
        con,
        "UPDATE ms_compound
         SET name = ?, subsid = ?, smiles = ?
         WHERE compound_id = ?",
        params = list(
          df$COMPNAME[i],
          df$SUBSID[i],
          df$SMILES[i],
          as.integer(df$COMPID[i])
        )
      )
    }
  }

  
  update_database(ref_db,    Updated_Names[[1]])
  update_database(target_db, Updated_Names[[2]])
  
  return(Mult.frame)
}

