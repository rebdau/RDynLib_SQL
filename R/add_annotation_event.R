#' @title Add annotations from a given CSV file to the SQL database.
#' 
#' @description
#' The 'add_annotation_event()' function creates the annotation_history table
#' in the SQL database if it does not already exist and populates it from a
#' given CSV file.
#' The CSV file should contain the following columns: expid, compound_id,
#' date, type, description, and references.
#' Both `reference` and `references` are accepted as CSV column names.
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
  
  # Connect to SQLite database
  con <- DBI::dbConnect(RSQLite::SQLite(), sql_path)
  on.exit(DBI::dbDisconnect(con), add = TRUE)
  
  # Read CSV (semicolon separated)
  annot <- readr::read_csv2(
    csv_file,
    show_col_types = FALSE,
    trim_ws = TRUE,
    col_types = readr::cols(
      .default = readr::col_character()
    )
  )
  
  # Remove unnamed columns
  annot <- annot[
    ,
    !grepl("^\\.\\.\\.[0-9]+$", names(annot)),
    drop = FALSE
  ]
  
  # Accept both reference and references in CSV
  if ("reference" %in% names(annot) &&
      !"references" %in% names(annot)) {
    
    names(annot)[names(annot) == "reference"] <- "references"
  }
  
  # Check required columns
  required <- c(
    "expid",
    "compound_id",
    "date",
    "type",
    "description",
    "references"
  )
  
  miss <- setdiff(required, names(annot))
  
  if (length(miss) > 0) {
    stop(
      "Missing columns: ",
      paste(miss, collapse = ", "),
      "\nColumns found:\n",
      paste(names(annot), collapse = ", ")
    )
  }
  
  # Validate required values
  if (any(is.na(annot$compound_id) |
          trimws(annot$compound_id) == "" |
          is.na(annot$date) |
          trimws(annot$date) == "")) {
    stop("compound_id and date must not contain missing values.")
  }
  
  # Validate compound IDs
  if (any(!grepl("^[0-9]+$", annot$compound_id))) {
    stop("compound_id must contain valid integer identifiers.")
  }
  
  # Create table and migrate existing schema
  DBI::dbWithTransaction(con, {
    
    # Create annotation_history if missing
    DBI::dbExecute(con, '
      CREATE TABLE IF NOT EXISTS annotation_history (
        annotation_id INTEGER PRIMARY KEY AUTOINCREMENT,
        expid INTEGER,
        compound_id INTEGER NOT NULL,
        date TEXT NOT NULL,
        type TEXT,
        description TEXT,
        "references" TEXT,

        FOREIGN KEY(compound_id)
          REFERENCES ms_compound(compound_id),

        UNIQUE(compound_id, date)
      );
    ')
    
    # Check existing columns
    existing_cols <- DBI::dbListFields(
      con,
      "annotation_history"
    )
    
    # Rename reference to references if needed
    if ("reference" %in% existing_cols &&
        !"references" %in% existing_cols) {
      
      DBI::dbExecute(
        con,
        'ALTER TABLE annotation_history
         RENAME COLUMN "reference" TO "references"'
      )
    }
    
    # Add references column if missing
    existing_cols <- DBI::dbListFields(
      con,
      "annotation_history"
    )
    
    if (!"references" %in% existing_cols) {
      
      DBI::dbExecute(
        con,
        'ALTER TABLE annotation_history
         ADD COLUMN "references" TEXT'
      )
    }
    
    # Check for existing duplicate annotations
    duplicates <- DBI::dbGetQuery(con, '
      SELECT compound_id, date, COUNT(*) AS n
      FROM annotation_history
      GROUP BY compound_id, date
      HAVING COUNT(*) > 1
    ')
    
    if (nrow(duplicates) > 0) {
      stop(
        "Cannot create UNIQUE constraint: ",
        "duplicate (compound_id, date) pairs exist ",
        "in annotation_history."
      )
    }
    
    # Create UNIQUE index for existing databases
    DBI::dbExecute(con, '
      CREATE UNIQUE INDEX IF NOT EXISTS
        idx_annotation_compound_date
      ON annotation_history(compound_id, date)
    ')
    
    # Get valid compound IDs
    valid_ids <- DBI::dbGetQuery(
      con,
      "SELECT compound_id FROM ms_compound"
    )$compound_id
    
    # Keep only valid compounds
    annot <- annot[
      annot$compound_id %in% as.character(valid_ids),
      ,
      drop = FALSE
    ]
    
    # SQL query for insert or update
    sql <- '
      INSERT INTO annotation_history
        (
          expid,
          compound_id,
          date,
          type,
          description,
          "references"
        )
      VALUES (?, ?, ?, ?, ?, ?)

      ON CONFLICT(compound_id, date)
      DO UPDATE SET
        expid = excluded.expid,
        type = excluded.type,
        description = excluded.description,
        "references" = excluded."references";
    '
    
    # Insert or update each annotation
    for (i in seq_len(nrow(annot))) {
      
      DBI::dbExecute(
        con,
        sql,
        params = list(
          annot$expid[i],
          annot$compound_id[i],
          annot$date[i],
          annot$type[i],
          annot$description[i],
          annot$references[i]
        )
      )
    }
  })
  
  # Return number of processed annotations
  invisible(nrow(annot))
}