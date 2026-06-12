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
  
  ref_db <- unique(Assoc$ref_database)
  target_db <- unique(Assoc$target_database)
  
  if (length(ref_db) != 1 || length(target_db) != 1) {
    stop(
      "Association file must contain exactly one ref_database and one target_database"
    )
  }
  
  ref_db <- ref_db[1]
  target_db <- target_db[1]
  
  cat("Reference database :", ref_db, "\n")
  cat("Target database    :", target_db, "\n")
  
  # Load metadata
  Names_lst <- Files_TransferNames_SQL(
    sqlite_dir = sqlite_dir,
    assoc_path = assoc_path
  )
  
  transfertNam <- TransferNames_SQL(Names_lst)
  
  Updated_Names <- transfertNam[[1]]
  Mult.frame <- transfertNam[[2]]
  
  update_database <- function(db_name, df) {
    
    sqlite_path <- file.path(sqlite_dir, db_name)
    
    cat("\nOpening :", sqlite_path, "\n")
    
    if (!file.exists(sqlite_path)) {
      stop(
        paste(
          "Database file not found:",
          normalizePath(sqlite_path, winslash = "/",
                        mustWork = FALSE)
        )
      )
    }
    
    con <- tryCatch(
      dbConnect(SQLite(), sqlite_path),
      error = function(e) {
        stop(
          paste(
            "Could not connect to database:",
            sqlite_path,
            "\n",
            conditionMessage(e)
          )
        )
      }
    )
    
    on.exit(dbDisconnect(con), add = TRUE)
    
    if (!"ms_compound" %in% dbListTables(con)) {
      stop(
        paste(
          "Table 'ms_compound' not found in",
          db_name
        )
      )
    }
    
    required_cols <- c(
      "compound_id",
      "name",
      "subsid",
      "smiles"
    )
    
    missing_cols <- setdiff(required_cols, colnames(df))
    
    if (length(missing_cols) > 0) {
      stop(
        paste(
          "Missing columns:",
          paste(missing_cols, collapse = ", ")
        )
      )
    }
    
    for (i in seq_len(nrow(df))) {
      
      dbExecute(
        con,
        "UPDATE ms_compound
         SET name = ?,
             subsid = ?,
             smiles = ?
         WHERE compound_id = ?",
        params = list(
          as.character(df$COMPNAME[i]),
          as.character(df$SUBSID[i]),
          as.character(df$SMILES[i]),
          as.integer(df$COMPID[i])
        )
      )
    }
    
    cat("Updated", nrow(df), "records in", db_name, "\n")
  }
  
  update_database(ref_db, Updated_Names[[1]])
  update_database(target_db, Updated_Names[[2]])
  
  return(Mult.frame)
}