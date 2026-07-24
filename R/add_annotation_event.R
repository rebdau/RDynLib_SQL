#' @title Add annotations from a given CSV file to the SQL database.
#' 
#' @description
#' The 'add_annotation_event()' function creates the annotation_history table
#' in the SQL database if it does not already exist and populates it from a
#' given CSV file.
#' The CSV file should contain the following columns: expid, compound_id,
#' date, type, description, and reference.
#' The function then creates an incremental primary key, 'annotation_id'.
#' If the same compound already exists in the annotation_history table and in
#' the CSV file with the same date, the corresponding row in the SQL table will
#' be updated. If they have different dates, a new row will be inserted.
#' 
#' @param sql_path Path to the SQL database.
#' @param csv_file Path to the CSV file containing the annotations.
#' 
#' @returns The new or updated 'annotation_history' table in the SQL database.
#' 
#' @import DBI
#' @import RSQLite
#' @import readr
#' 
#' @author Ahlam Mentag
#' 
#' @export
add_annotation_event <- function(sql_path, csv_file) {
  
  con <- dbConnect(SQLite(), sql_path)
  on.exit(dbDisconnect(con), add = TRUE)
  
  # Create table if it does not exist
  dbExecute(con, "
    CREATE TABLE IF NOT EXISTS annotation_history (
      annotation_id INTEGER PRIMARY KEY AUTOINCREMENT,
      expid INTEGER,
      compound_id INTEGER NOT NULL,
      date TEXT NOT NULL,
      type TEXT,
      description TEXT,
      reference TEXT,

      FOREIGN KEY(compound_id)
        REFERENCES ms_compound(compound_id),

      UNIQUE(compound_id, date)
    );
  ")
  
  # Read csv (semicolon separated)
  annot <- read_csv2(
    csv_file,
    show_col_types = FALSE,
    trim_ws = TRUE
  )
  
  annot <- annot[, names(annot) != ""]
  
  ## Check required columns
  required <- c(
    "expid",
    "compound_id",
    "date",
    "type",
    "description",
    "reference"
  )
  
  miss <- setdiff(required, names(annot))
  
  if (length(miss) > 0)
    stop(
      "Missing columns: ",
      paste(miss, collapse = ", "),
      "\nColumns found:\n",
      paste(names(annot), collapse = ", ")
    )
  
  ## Keep only valid compounds
  valid_ids <- dbGetQuery(
    con,
    "SELECT compound_id FROM ms_compound"
  )$compound_id
  
  annot <- annot[annot$compound_id %in% valid_ids, ]
  
  sql <- "
    INSERT INTO annotation_history
      (expid, compound_id, date, type, description, reference)
    VALUES
      (?, ?, ?, ?, ?, ?)

    ON CONFLICT(compound_id, date)
    DO UPDATE SET
      expid = excluded.expid,
      type = excluded.type,
      description = excluded.description,
      reference = excluded.reference;
  "
  
  dbBegin(con)
  
  for(i in seq_len(nrow(annot))) {
    
    dbExecute(
      con,
      sql,
      params = list(
        annot$expid[i],
        annot$compound_id[i],
        annot$date[i],
        annot$type[i],
        annot$description[i],
        annot$reference[i]
      )
    )
    
  }
  
  dbCommit(con)
  
  invisible(nrow(annot))
}