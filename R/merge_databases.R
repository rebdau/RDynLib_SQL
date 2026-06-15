#' @title Merge two SQLite spectral databases.
#'
#' @description
#' Merges two SQLite databases containing mass spectrometry data into a
#' single output database. The function aligns tables present in both
#' databases, shifts identifier columns to avoid ID conflicts,
#' and combines the records. Tables that exist in only
#' one database are copied directly to the output database.
#'
#' Identifier columns (`expid`, `compound_id`, `spectrum_id`, `peak_id`)
#' are automatically shifted to ensure uniqueness across the merged data.
#' Retention times (`rtime`, `retention_time`) from the additional database
#' are converted from seconds to minutes before merging.
#'
#' @param main_db 'Character(1)' Path to the primary SQLite database.
#' 
#' @param add_db 'Character(1)' Path to the second SQLite database that would be 
#' added to the main database.
#' 
#' @param output_db Character. Path where the merged SQLite database will
#'   be written.
#'
#' @return No value is returned. The function writes the merged database
#'   to `output_db` and prints progress messages during execution.
#'
#' @import DBI
#' @import RSQLite
#' @import dplyr
#' 
#' @author Ahlam Mentag
#' 
#' @export
merge_databases <- function(main_db, add_db, output_db) {
  
  con_main <- dbConnect(SQLite(), main_db)
  con_add  <- dbConnect(SQLite(), add_db)
  con_out  <- dbConnect(SQLite(), output_db)
  
  on.exit({
    dbDisconnect(con_main)
    dbDisconnect(con_add)
    dbDisconnect(con_out)
  })
  
  
  # TABLE LISTS 
  
  tables_main <- setdiff(dbListTables(con_main), "sqlite_sequence")
  tables_add  <- setdiff(dbListTables(con_add), "sqlite_sequence")
  
  tables_main_l <- tolower(tables_main)
  tables_add_l  <- tolower(tables_add)
  
  all_tables <- union(tables_main_l, tables_add_l)
  
  required_tables <- c(
    "synonym",
    "experiment",
    "ms_compound",
    "msms_spectrum",
    "msms_spectrum_peak"
  )
  
  all_tables <- union(all_tables, required_tables)
  
  
  # SAFE READ
  
  safe_read <- function(con, tables_vec, tbl_name) {
    hit <- tables_vec[tolower(tables_vec) == tbl_name]
    if (length(hit) == 0) return(NULL)
    dbReadTable(con, hit[1])
  }
  
  
  # SHIFT IDS
  
  shift_ids <- function(df, exp_shift, spectrum_shift, peak_shift, compound_shift) {
    
    if ("expid" %in% names(df))
      df$expid <- as.integer(df$expid) + exp_shift
    
    if ("spectrum_id" %in% names(df))
      df$spectrum_id <- as.integer(df$spectrum_id) + spectrum_shift
    
    if ("peak_id" %in% names(df))
      df$peak_id <- as.integer(df$peak_id) + peak_shift
    
    if ("compound_id" %in% names(df))
      df$compound_id <- as.character(as.integer(df$compound_id) + compound_shift)
    
    if ("subsid" %in% names(df))
      df$subsid <- as.integer(df$subsid)
    
    df
  }
  
  
  # SHIFT CALCULATION
  
  compute_shift <- function(con_main, con_add, table, key) {
    
    main_tbl <- safe_read(con_main, tables_main, table)
    add_tbl  <- safe_read(con_add, tables_add, table)
    
    if (is.null(main_tbl) || is.null(add_tbl)) return(0)
    if (!key %in% names(main_tbl) || !key %in% names(add_tbl)) return(0)
    
    main_vals <- suppressWarnings(as.integer(main_tbl[[key]]))
    add_vals  <- suppressWarnings(as.integer(add_tbl[[key]]))
    
    max_main <- if (any(!is.na(main_vals))) max(main_vals, na.rm = TRUE) else 0
    min_add  <- if (any(!is.na(add_vals)))  min(add_vals, na.rm = TRUE) else 0
    
    max_main + 1 - min_add
  }
  
  exp_shift      <- compute_shift(con_main, con_add, "experiment", "expid")
  compound_shift <- compute_shift(con_main, con_add, "ms_compound", "compound_id")
  spectrum_shift <- compute_shift(con_main, con_add, "msms_spectrum", "spectrum_id")
  peak_shift     <- compute_shift(con_main, con_add, "msms_spectrum_peak", "peak_id")
  
  
  # SCHEMA 
  
  FTMS_SCHEMA <- list(
    expid = "integer",
    spectrum_id = "integer",
    peak_id = "integer",
    compound_id = "character",
    subsid = "integer",
    rtime = "numeric",
    retention_time = "numeric",
    ppm_deviation = "numeric"   
  )
  
  
  # TYPE FIX
  
  force_schema <- function(df, schema) {
    for (col in names(schema)) {
      if (!col %in% names(df)) next
      
      target <- schema[[col]]
      
      df[[col]] <- switch(
        target,
        character = as.character(df[[col]]),
        numeric   = as.numeric(df[[col]]),
        integer   = as.integer(df[[col]]),
        df[[col]]
      )
    }
    df
  }
  
  
  # STRONG TYPE NORMALIZER (IMPORTANT FIX)
  
  normalize_types <- function(df1, df2) {
    
    common_cols <- intersect(names(df1), names(df2))
    
    for (col in common_cols) {
      
      # if either side is numeric-like, force numeric
      if (is.numeric(df1[[col]]) || is.numeric(df2[[col]])) {
        df1[[col]] <- suppressWarnings(as.numeric(df1[[col]]))
        df2[[col]] <- suppressWarnings(as.numeric(df2[[col]]))
      } else {
        df1[[col]] <- as.character(df1[[col]])
        df2[[col]] <- as.character(df2[[col]])
      }
    }
    
    list(df1 = df1, df2 = df2)
  }
  
  
  # TABLE LOOP
  
  for (tbl_name in all_tables) {
    
    tbl_main <- safe_read(con_main, tables_main, tbl_name)
    tbl_add  <- safe_read(con_add, tables_add, tbl_name)
    
    if (is.null(tbl_main) && is.null(tbl_add)) next
    
    if (is.null(tbl_main)) tbl_main <- tbl_add[0, , drop = FALSE]
    if (is.null(tbl_add))  tbl_add  <- tbl_main[0, , drop = FALSE]
    
    # align columns
    all_cols <- union(names(tbl_main), names(tbl_add))
    tbl_main <- tbl_main[, all_cols, drop = FALSE]
    tbl_add  <- tbl_add[, all_cols, drop = FALSE]
    
    # schema enforcement
    tbl_main <- force_schema(tbl_main, FTMS_SCHEMA)
    tbl_add  <- force_schema(tbl_add, FTMS_SCHEMA)
    
    # CRITICAL FIX: harmonize types before bind_rows
    fixed <- normalize_types(tbl_main, tbl_add)
    tbl_main <- fixed$df1
    tbl_add  <- fixed$df2
    
    # shift add DB only
    if (nrow(tbl_add) > 0) {
      
      tbl_add <- shift_ids(tbl_add, exp_shift, spectrum_shift, peak_shift, compound_shift)
      
      if ("rtime" %in% names(tbl_add))
        tbl_add$rtime <- as.numeric(tbl_add$rtime) / 60
      
      if ("retention_time" %in% names(tbl_add))
        tbl_add$retention_time <- as.numeric(tbl_add$retention_time) / 60
    }
    
    merged_tbl <- bind_rows(tbl_main, tbl_add)
    
    dbWriteTable(con_out, tbl_name, merged_tbl, overwrite = TRUE)
    
    cat("Merged:", tbl_name, "\n")
  }
  
  
  # UNIQUE TABLES FROM MAIN
  
  for (tbl_name in setdiff(tables_main_l, tables_add_l)) {
    tbl <- dbReadTable(con_main, tables_main[tolower(tables_main) == tbl_name])
    dbWriteTable(con_out, tbl_name, tbl, overwrite = TRUE)
  }
  
  
  # UNIQUE TABLES FROM ADD
  
  for (tbl_name in setdiff(tables_add_l, tables_main_l)) {
    
    tbl <- dbReadTable(con_add, tables_add[tolower(tables_add) == tbl_name])
    
    if (nrow(tbl) > 0)
      tbl <- shift_ids(tbl, exp_shift, spectrum_shift, peak_shift, compound_shift)
    
    dbWriteTable(con_out, tbl_name, tbl, overwrite = TRUE)
  }
  
  cat("Database merge completed. Output saved to:", output_db, "\n")
}