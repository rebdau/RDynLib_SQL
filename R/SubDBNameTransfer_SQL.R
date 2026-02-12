#' @title Transfer compound names across SQLite sub-databases.
#' 
#' @description SubDBNameTransfer_SQL() function propagates compound names, 
#' identifiers, and SMILES strings between associated compounds stored
#' in a SQLite database. Associations are provided via a single association file
#' defining reference and target compound IDs.
#'
#' @param sqlite_path Character string. Path to the SQLite database containing
#'   the \code{ms_compound} table.
#' @param assoc_path Character string. Path to the association file defining
#'   reference and target compound_ids of two aligned compounds.
#'
#' @return Updates compound names, SUBSID, and SMILES directly in the SQLite
#'   database using compound IDs.
#'
#' @author Ahlam Mentag
#' 
#' @export
SubDBNameTransfer_SQL <- function(sqlite_path, assoc_path) {
  
  library(DBI)
  library(RSQLite)
  
  # Load names and associations
  Names_lst <- Files_TransferNames_SQL(
    sqlite_path = sqlite_path,
    assoc_path  = assoc_path
  )
  
  #Transfer names
  transfertNam <- TransferNames_SQL(Names_lst)
  
  Updated_Names <- transfertNam[[1]]
  Mult.frame    <- transfertNam[[2]]

  # Write back to SQLite
  con <- dbConnect(SQLite(), sqlite_path)
  on.exit(dbDisconnect(con), add = TRUE)
  
  update_table <- function(df) {
    
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
  
  update_table(Updated_Names[[1]])  # ref
  update_table(Updated_Names[[2]])  # target
  update_table(Updated_Names[[3]])  # syn (if needed)
  
  return(Mult.frame)
}
