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
  
  library(DBI)
  library(RSQLite)
  library(dplyr)
  
  if (file.exists(output_db))
    file.remove(output_db)
  
  con_main <- dbConnect(SQLite(), main_db)
  con_add  <- dbConnect(SQLite(), add_db)
  con_out  <- dbConnect(SQLite(), output_db)
  
  on.exit({
    dbDisconnect(con_main)
    dbDisconnect(con_add)
    dbDisconnect(con_out)
  })
  
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
  
  safe_read <- function(con, tables_vec, tbl_name) {
    
    hit <- tables_vec[tolower(tables_vec) == tbl_name]
    
    if (length(hit) == 0)
      return(NULL)
    
    dbReadTable(con, hit[1])
  }
  
  shift_ids <- function(df,
                        exp_shift,
                        spectrum_shift,
                        peak_shift,
                        compound_shift) {
    
    if ("expid" %in% names(df))
      df$expid <- as.integer(df$expid) + exp_shift
    
    if ("spectrum_id" %in% names(df))
      df$spectrum_id <- as.integer(df$spectrum_id) + spectrum_shift
    
    if ("peak_id" %in% names(df))
      df$peak_id <- as.integer(df$peak_id) + peak_shift
    
    if ("compound_id" %in% names(df))
      df$compound_id <- as.character(
        as.integer(df$compound_id) + compound_shift
      )
    
    if ("subsid" %in% names(df))
      df$subsid <- as.integer(df$subsid)
    
    df
  }
  
  compute_shift <- function(con_main,
                            con_add,
                            table,
                            key) {
    
    main_tbl <- safe_read(con_main, tables_main, table)
    add_tbl  <- safe_read(con_add, tables_add, table)
    
    if (is.null(main_tbl) || is.null(add_tbl))
      return(0)
    
    if (!key %in% names(main_tbl))
      return(0)
    
    if (!key %in% names(add_tbl))
      return(0)
    
    main_vals <- suppressWarnings(
      as.integer(main_tbl[[key]])
    )
    
    add_vals <- suppressWarnings(
      as.integer(add_tbl[[key]])
    )
    
    max_main <- if (length(main_vals))
      max(main_vals, na.rm = TRUE) else 0
    
    min_add <- if (length(add_vals))
      min(add_vals, na.rm = TRUE) else 0
    
    max_main + 1 - min_add
  }
  
  exp_shift <- compute_shift(
    con_main, con_add,
    "experiment", "expid"
  )
  
  compound_shift <- compute_shift(
    con_main, con_add,
    "ms_compound", "compound_id"
  )
  
  spectrum_shift <- compute_shift(
    con_main, con_add,
    "msms_spectrum", "spectrum_id"
  )
  
  peak_shift <- compute_shift(
    con_main, con_add,
    "msms_spectrum_peak", "peak_id"
  )
  
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
  
  force_schema <- function(df, schema) {
    
    if (is.null(df))
      return(df)
    
    for (col in names(schema)) {
      
      if (!col %in% names(df))
        next
      
      df[[col]] <- switch(
        schema[[col]],
        character = as.character(df[[col]]),
        integer   = as.integer(df[[col]]),
        numeric   = as.numeric(df[[col]]),
        df[[col]]
      )
    }
    
    df
  }
  
  normalize_types <- function(df1, df2) {
    
    common_cols <- intersect(
      names(df1),
      names(df2)
    )
    
    for (col in common_cols) {
      
      cls1 <- class(df1[[col]])[1]
      cls2 <- class(df2[[col]])[1]
      
      if (cls1 %in% c("numeric", "integer") ||
          cls2 %in% c("numeric", "integer")) {
        
        df1[[col]] <- suppressWarnings(
          as.numeric(df1[[col]])
        )
        
        df2[[col]] <- suppressWarnings(
          as.numeric(df2[[col]])
        )
        
      } else {
        
        df1[[col]] <- as.character(df1[[col]])
        df2[[col]] <- as.character(df2[[col]])
      }
    }
    
    list(df1 = df1, df2 = df2)
  }
  
  for (tbl_name in all_tables) {
    
    cat("Processing:", tbl_name, "\n")
    
    tbl_main <- safe_read(
      con_main,
      tables_main,
      tbl_name
    )
    
    tbl_add <- safe_read(
      con_add,
      tables_add,
      tbl_name
    )
    
    if (is.null(tbl_main) && is.null(tbl_add))
      next
    
    if (is.null(tbl_main))
      tbl_main <- tbl_add[0, , drop = FALSE]
    
    if (is.null(tbl_add))
      tbl_add <- tbl_main[0, , drop = FALSE]
    
    tbl_main <- force_schema(
      tbl_main,
      FTMS_SCHEMA
    )
    
    tbl_add <- force_schema(
      tbl_add,
      FTMS_SCHEMA
    )
    
    fixed <- normalize_types(
      tbl_main,
      tbl_add
    )
    
    tbl_main <- fixed$df1
    tbl_add  <- fixed$df2
    
    if (nrow(tbl_add) > 0) {
      
      tbl_add <- shift_ids(
        tbl_add,
        exp_shift,
        spectrum_shift,
        peak_shift,
        compound_shift
      )
      
      if ("rtime" %in% names(tbl_add))
        tbl_add$rtime <- as.numeric(tbl_add$rtime) / 60
      
      if ("retention_time" %in% names(tbl_add))
        tbl_add$retention_time <- as.numeric(tbl_add$retention_time) / 60
    }
    
    merged_tbl <- bind_rows(
      tbl_main,
      tbl_add
    )
    
    dbWriteTable(
      con_out,
      tbl_name,
      merged_tbl,
      overwrite = TRUE
    )
    
    cat("Merged:", tbl_name, "\n")
  }
  
  cat(
    "\nDatabase merge completed.\nOutput:",
    output_db,
    "\n"
  )
}