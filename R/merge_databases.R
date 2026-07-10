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
  
  if (file.exists(output_db)) {
    
    con_test <- tryCatch(
      DBI::dbConnect(RSQLite::SQLite(), output_db),
      error = function(e) NULL
    )
    
    if (!is.null(con_test)) {
      DBI::dbDisconnect(con_test)
      file.remove(output_db)
    } else {
      stop("Output database is locked. Close all connections first.")
    }
  }
  
  
  con_main <- DBI::dbConnect(RSQLite::SQLite(), main_db)
  con_add  <- DBI::dbConnect(RSQLite::SQLite(), add_db)
  con_out  <- DBI::dbConnect(RSQLite::SQLite(), output_db)
  
  
  on.exit({
    DBI::dbDisconnect(con_main)
    DBI::dbDisconnect(con_add)
    DBI::dbDisconnect(con_out)
  })
  
  
  tables_main <- setdiff(DBI::dbListTables(con_main), "sqlite_sequence")
  tables_add  <- setdiff(DBI::dbListTables(con_add), "sqlite_sequence")
  
  all_tables <- union(tables_main, tables_add)
  
  
  safe_read <- function(con, tables_vec, tbl_name) {
    
    hit <- tables_vec[
      tolower(tables_vec) == tolower(tbl_name)
    ]
    
    if (length(hit) == 0)
      return(NULL)
    
    DBI::dbReadTable(con, hit[1])
  }
  
  
  get_schema <- function(con, tbl_name) {
    
    res <- DBI::dbGetQuery(
      con,
      sprintf(
        "SELECT sql
         FROM sqlite_master
         WHERE type='table'
         AND lower(name)=lower('%s')",
        tbl_name
      )
    )
    
    if (nrow(res) != 1)
      return(NULL)
    
    res$sql[[1]]
  }
  
  
  normalize_types <- function(main_df, add_df) {
    
    common <- intersect(
      names(main_df),
      names(add_df)
    )
    
    for (col in common) {
      
      if (col == "expid") {
        main_df[[col]] <- as.character(main_df[[col]])
        add_df[[col]]  <- as.character(add_df[[col]])
        next
      }
      
      
      cl <- class(main_df[[col]])[1]
      
      if (cl == "character") {
        
        main_df[[col]] <- as.character(main_df[[col]])
        add_df[[col]]  <- as.character(add_df[[col]])
        
      } else if (cl == "integer") {
        
        main_df[[col]] <- as.integer(main_df[[col]])
        add_df[[col]]  <- as.integer(add_df[[col]])
        
      } else if (cl == "numeric") {
        
        main_df[[col]] <- as.numeric(main_df[[col]])
        add_df[[col]]  <- as.numeric(add_df[[col]])
      }
    }
    
    list(
      main = main_df,
      add = add_df
    )
  }
  
  
  compute_shift <- function(table, key) {
    
    main_tbl <- safe_read(con_main, tables_main, table)
    add_tbl  <- safe_read(con_add, tables_add, table)
    
    if (is.null(main_tbl) || is.null(add_tbl))
      return(0)
    
    if (!key %in% names(main_tbl) ||
        !key %in% names(add_tbl))
      return(0)
    
    main_val <- suppressWarnings(as.integer(main_tbl[[key]]))
    add_val  <- suppressWarnings(as.integer(add_tbl[[key]]))
    
    if (length(main_val) == 0 ||
        length(add_val) == 0)
      return(0)
    
    max(main_val, na.rm = TRUE) + 1 -
      min(add_val, na.rm = TRUE)
  }
  
  
  spectrum_shift <- compute_shift(
    "msms_spectrum",
    "spectrum_id"
  )
  
  peak_shift <- compute_shift(
    "msms_spectrum_peak",
    "peak_id"
  )
  
  compound_shift <- compute_shift(
    "ms_compound",
    "compound_id"
  )
  
  experiment_shift <- compute_shift(
    "experiment",
    "expid"
  )
  
  shift_ids <- function(df) {
    
    if ("spectrum_id" %in% names(df))
      df$spectrum_id <- as.integer(df$spectrum_id) + spectrum_shift
    
    if ("peak_id" %in% names(df))
      df$peak_id <- as.integer(df$peak_id) + peak_shift
    
    if ("compound_id" %in% names(df))
      df$compound_id <- as.character(
        as.integer(df$compound_id) + compound_shift
      )
    
    if ("expid" %in% names(df))
      df$expid <- as.character(
        as.integer(df$expid) + experiment_shift
      )
    
    df
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
    
    
    if (!is.null(tbl_main) &&
        is.null(tbl_add)) {
      
      schema <- get_schema(
        con_main,
        tbl_name
      )
      
      DBI::dbExecute(con_out, schema)
      
      DBI::dbAppendTable(
        con_out,
        tbl_name,
        tbl_main
      )
      
      cat("Copied main:", tbl_name, "\n")
      next
    }
    
    
    if (is.null(tbl_main) &&
        !is.null(tbl_add)) {
      
      schema <- get_schema(
        con_add,
        tbl_name
      )
      
      DBI::dbExecute(con_out, schema)
      
      DBI::dbAppendTable(
        con_out,
        tbl_name,
        tbl_add
      )
      
      cat("Copied add:", tbl_name, "\n")
      next
    }
    
    
    fixed <- normalize_types(
      tbl_main,
      tbl_add
    )
    
    tbl_main <- fixed$main
    tbl_add  <- fixed$add
    
    
    tbl_add <- shift_ids(tbl_add)
    
    
    if ("rtime" %in% names(tbl_add))
      tbl_add$rtime <-
      as.numeric(tbl_add$rtime) / 60
    
    
    if ("retention_time" %in% names(tbl_add))
      tbl_add$retention_time <-
      as.numeric(tbl_add$retention_time) / 60
    
    
    merged <- dplyr::bind_rows(
      tbl_main,
      tbl_add
    )
    
    
    schema <- get_schema(
      con_main,
      tbl_name
    )
    
    
    if (tbl_name == "experiment") {
      
      schema <- gsub(
        "`expid` INTEGER",
        "`expid` TEXT",
        schema
      )
    }
    
    
    DBI::dbExecute(
      con_out,
      schema
    )
    
    
    DBI::dbAppendTable(
      con_out,
      tbl_name,
      merged
    )
    
    
    cat("Merged:", tbl_name, "\n")
  }
  
  
  cat(
    "\nDatabase merge completed.\nOutput:",
    output_db,
    "\n"
  )
}