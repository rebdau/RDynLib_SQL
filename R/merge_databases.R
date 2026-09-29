#' @title Merge two SQLite spectral databases
#'
#' @description
#' Merges two SQLite databases containing mass spectrometry data into a
#' single output database. Tables present in both databases are merged,
#' while tables present in only one database are copied to the output.
#' If the same table contains different columns between databases, all
#' columns are retained and missing values are filled with `NA`.
#'
#' Identifier columns (`expid`, `compound_id`, `spectrum_id`, and `peak_id`)
#' from the additional database are shifted when necessary to avoid ID
#' conflicts.
#'
#' Retention times (`rtime` and `retention_time`) from the additional
#' database can optionally be converted from seconds to minutes.
#'
#' @param main_db `Character(1)`. Path to the primary SQLite database.
#'
#' @param add_db `Character(1)`. Path to the SQLite database to add to
#'   the primary database.
#'
#' @param output_db `Character(1)`. Path where the merged SQLite database
#'   will be written.
#'
#' @param convert_add_rt `Logical(1)`. If `TRUE`, retention times
#'   (`rtime` and `retention_time`) from `add_db` are converted from
#'   seconds to minutes by dividing them by 60. Default is `FALSE`.
#'
#' @return No value is returned. The merged SQLite database is written
#'   to `output_db`, and progress messages are printed during execution.
#'
#' @import DBI
#' @import RSQLite
#' @import dplyr
#'
#' @author Ahlam Mentag
#'
#' @export
merge_databases <- function(
    main_db,
    add_db,
    output_db,
    convert_add_rt = FALSE
) {
  

  
  if (file.exists(output_db)) {
    
    con_test <- tryCatch(
      DBI::dbConnect(RSQLite::SQLite(), output_db),
      error = function(e) NULL
    )
    
    if (!is.null(con_test)) {
      
      DBI::dbDisconnect(con_test)
      
      if (!file.remove(output_db)) {
        stop("Could not remove existing output database.")
      }
      
    } else {
      stop("Output database is locked. Close all connections first.")
    }
  }
  
  

  
  con_main <- DBI::dbConnect(
    RSQLite::SQLite(),
    main_db
  )
  
  con_add <- DBI::dbConnect(
    RSQLite::SQLite(),
    add_db
  )
  
  con_out <- DBI::dbConnect(
    RSQLite::SQLite(),
    output_db
  )
  
  
  on.exit({
    
    try(DBI::dbDisconnect(con_main), silent = TRUE)
    try(DBI::dbDisconnect(con_add),  silent = TRUE)
    try(DBI::dbDisconnect(con_out),  silent = TRUE)
    
  }, add = TRUE)
  
  

  
  tables_main <- setdiff(
    DBI::dbListTables(con_main),
    "sqlite_sequence"
  )
  
  tables_add <- setdiff(
    DBI::dbListTables(con_add),
    "sqlite_sequence"
  )
  
  
  all_table_names <- unique(
    tolower(c(
      tables_main,
      tables_add
    ))
  )
  
  

  
  find_table <- function(tables_vec, tbl_name) {
    
    hit <- tables_vec[
      tolower(tables_vec) == tolower(tbl_name)
    ]
    
    if (length(hit) == 0) {
      return(NULL)
    }
    
    hit[[1]]
  }
  
  

  
  safe_read <- function(
    con,
    tables_vec,
    tbl_name
  ) {
    
    real_name <- find_table(
      tables_vec,
      tbl_name
    )
    
    if (is.null(real_name)) {
      return(NULL)
    }
    
    DBI::dbReadTable(
      con,
      real_name
    )
  }
  
  
  # Get SQLite table information
  
  get_table_info <- function(
    con,
    tables_vec,
    tbl_name
  ) {
    
    real_name <- find_table(
      tables_vec,
      tbl_name
    )
    
    if (is.null(real_name)) {
      return(NULL)
    }
    
    DBI::dbGetQuery(
      con,
      paste0(
        "PRAGMA table_info(",
        DBI::dbQuoteIdentifier(
          con,
          real_name
        ),
        ")"
      )
    )
  }
  
  
  # Determine SQLite column type
  
  get_column_type <- function(
    col,
    info_main,
    info_add
  ) {
    
    # expid is stored as TEXT in merged database
    if (tolower(col) == "expid") {
      return("TEXT")
    }
    
    
    # Prefer type from MAIN database
    if (!is.null(info_main)) {
      
      idx <- which(
        tolower(info_main$name) ==
          tolower(col)
      )
      
      if (length(idx) > 0) {
        
        tp <- info_main$type[idx[1]]
        
        if (
          !is.na(tp) &&
          nzchar(tp)
        ) {
          return(tp)
        }
      }
    }
    
    
    # Otherwise use ADD database
    if (!is.null(info_add)) {
      
      idx <- which(
        tolower(info_add$name) ==
          tolower(col)
      )
      
      if (length(idx) > 0) {
        
        tp <- info_add$type[idx[1]]
        
        if (
          !is.na(tp) &&
          nzchar(tp)
        ) {
          return(tp)
        }
      }
    }
    
    
    # Fallback
    "TEXT"
  }
  
  

  # Create output table using UNION of columns
  
  create_union_table <- function(
    con_out,
    tbl_name,
    tbl_main,
    tbl_add,
    info_main,
    info_add
  ) {
    
    columns <- union(
      if (!is.null(tbl_main)) {
        names(tbl_main)
      } else {
        character()
      },
      if (!is.null(tbl_add)) {
        names(tbl_add)
      } else {
        character()
      }
    )
    
    
    if (length(columns) == 0) {
      stop(
        "Table ",
        tbl_name,
        " has no columns."
      )
    }
    
    
    definitions <- vapply(
      columns,
      function(col) {
        
        type <- get_column_type(
          col,
          info_main,
          info_add
        )
        
        paste(
          DBI::dbQuoteIdentifier(
            con_out,
            col
          ),
          type
        )
      },
      character(1)
    )
    
    
    sql <- paste0(
      "CREATE TABLE ",
      DBI::dbQuoteIdentifier(
        con_out,
        tbl_name
      ),
      " (\n",
      paste(
        definitions,
        collapse = ",\n"
      ),
      "\n)"
    )
    
    
    DBI::dbExecute(
      con_out,
      sql
    )
  }
  
  
  #  Harmonize columns
  
  harmonize_columns <- function(
    main_df,
    add_df
  ) {
    
    all_cols <- union(
      if (!is.null(main_df)) {
        names(main_df)
      } else {
        character()
      },
      if (!is.null(add_df)) {
        names(add_df)
      } else {
        character()
      }
    )
    
    
    # Add missing columns to MAIN
    if (!is.null(main_df)) {
      
      missing_main <- setdiff(
        all_cols,
        names(main_df)
      )
      
      for (col in missing_main) {
        main_df[[col]] <- NA
      }
      
      main_df <- main_df[
        ,
        all_cols,
        drop = FALSE
      ]
    }
    
    
    # Add missing columns to ADD
    if (!is.null(add_df)) {
      
      missing_add <- setdiff(
        all_cols,
        names(add_df)
      )
      
      for (col in missing_add) {
        add_df[[col]] <- NA
      }
      
      add_df <- add_df[
        ,
        all_cols,
        drop = FALSE
      ]
    }
    
    
    list(
      main = main_df,
      add = add_df
    )
  }
  
  
  #  Normalize R column types
  
  normalize_types <- function(
    main_df,
    add_df
  ) {
    
    if (
      is.null(main_df) ||
      is.null(add_df)
    ) {
      
      return(
        list(
          main = main_df,
          add = add_df
        )
      )
    }
    
    
    common <- intersect(
      names(main_df),
      names(add_df)
    )
    
    
    for (col in common) {
      
      # IDs stored as character
      if (
        tolower(col) %in%
        c(
          "expid",
          "compound_id"
        )
      ) {
        
        main_df[[col]] <-
          as.character(main_df[[col]])
        
        add_df[[col]] <-
          as.character(add_df[[col]])
        
        next
      }
      
      
      # Character has priority
      if (
        is.character(main_df[[col]]) ||
        is.character(add_df[[col]])
      ) {
        
        main_df[[col]] <-
          as.character(main_df[[col]])
        
        add_df[[col]] <-
          as.character(add_df[[col]])
        
      } else if (
        is.numeric(main_df[[col]]) ||
        is.numeric(add_df[[col]])
      ) {
        
        main_df[[col]] <-
          as.numeric(main_df[[col]])
        
        add_df[[col]] <-
          as.numeric(add_df[[col]])
      }
    }
    
    
    list(
      main = main_df,
      add = add_df
    )
  }
  
  
  #  Compute ID shifts
  
  compute_shift <- function(
    table,
    key
  ) {
    
    main_tbl <- safe_read(
      con_main,
      tables_main,
      table
    )
    
    add_tbl <- safe_read(
      con_add,
      tables_add,
      table
    )
    
    
    if (
      is.null(main_tbl) ||
      is.null(add_tbl)
    ) {
      return(0L)
    }
    
    
    if (
      !key %in% names(main_tbl) ||
      !key %in% names(add_tbl)
    ) {
      return(0L)
    }
    
    
    main_val <- suppressWarnings(
      as.integer(
        main_tbl[[key]]
      )
    )
    
    add_val <- suppressWarnings(
      as.integer(
        add_tbl[[key]]
      )
    )
    
    
    main_val <- main_val[
      !is.na(main_val)
    ]
    
    add_val <- add_val[
      !is.na(add_val)
    ]
    
    
    if (
      length(main_val) == 0 ||
      length(add_val) == 0
    ) {
      return(0L)
    }
    
    
    as.integer(
      max(main_val) +
        1L -
        min(add_val)
    )
  }
  
  
  #  Calculate shifts
  
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
  
  
  cat(
    "\nID shifts:\n"
  )
  
  cat(
    " spectrum_id :",
    spectrum_shift,
    "\n"
  )
  
  cat(
    " peak_id     :",
    peak_shift,
    "\n"
  )
  
  cat(
    " compound_id :",
    compound_shift,
    "\n"
  )
  
  cat(
    " expid       :",
    experiment_shift,
    "\n\n"
  )
  
  
  #Shift IDs from ADD database
  
  shift_ids <- function(df) {
    
    if (is.null(df)) {
      return(NULL)
    }
    
    
    # spectrum_id
    if ("spectrum_id" %in% names(df)) {
      
      x <- suppressWarnings(
        as.integer(
          df$spectrum_id
        )
      )
      
      df$spectrum_id <- ifelse(
        is.na(x),
        NA_integer_,
        x + spectrum_shift
      )
    }
    
    
    # peak_id
    if ("peak_id" %in% names(df)) {
      
      x <- suppressWarnings(
        as.integer(
          df$peak_id
        )
      )
      
      df$peak_id <- ifelse(
        is.na(x),
        NA_integer_,
        x + peak_shift
      )
    }
    
    
    # compound_id
    if ("compound_id" %in% names(df)) {
      
      x <- suppressWarnings(
        as.integer(
          df$compound_id
        )
      )
      
      df$compound_id <- ifelse(
        is.na(x),
        NA_character_,
        as.character(
          x + compound_shift
        )
      )
    }
    
    
    # expid
    if ("expid" %in% names(df)) {
      
      x <- suppressWarnings(
        as.integer(
          df$expid
        )
      )
      
      df$expid <- ifelse(
        is.na(x),
        NA_character_,
        as.character(
          x + experiment_shift
        )
      )
    }
    
    
    df
  }
  
  
  # Process every table
  
  for (tbl_lower in all_table_names) {
    
    main_name <- find_table(
      tables_main,
      tbl_lower
    )
    
    add_name <- find_table(
      tables_add,
      tbl_lower
    )
    
    
    # Prefer MAIN table name
    if (!is.null(main_name)) {
      tbl_name <- main_name
    } else {
      tbl_name <- add_name
    }
    
    
    cat(
      "\nProcessing:",
      tbl_name,
      "\n"
    )
    
    
    # Read tables
    
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
    
    
    info_main <- get_table_info(
      con_main,
      tables_main,
      tbl_name
    )
    
    info_add <- get_table_info(
      con_add,
      tables_add,
      tbl_name
    )
    
    
    # Shift IDs from ADD database
    
    tbl_add <- shift_ids(
      tbl_add
    )
    
    
    # Convert ADD retention time from seconds -> minutes

    if (
      convert_add_rt &&
      !is.null(tbl_add)
    ) {
      
      if ("rtime" %in% names(tbl_add)) {
        
        tbl_add$rtime <-
          as.numeric(
            tbl_add$rtime
          ) / 60
      }
      
      
      if ("retention_time" %in% names(tbl_add)) {
        
        tbl_add$retention_time <-
          as.numeric(
            tbl_add$retention_time
          ) / 60
      }
    }

    
    fixed_cols <- harmonize_columns(
      tbl_main,
      tbl_add
    )
    
    tbl_main <- fixed_cols$main
    tbl_add  <- fixed_cols$add
    

    
    fixed_types <- normalize_types(
      tbl_main,
      tbl_add
    )
    
    tbl_main <- fixed_types$main
    tbl_add  <- fixed_types$add
    

    
    create_union_table(
      con_out,
      tbl_name,
      tbl_main,
      tbl_add,
      info_main,
      info_add
    )

    
    if (
      !is.null(tbl_main) &&
      nrow(tbl_main) > 0
    ) {
      
      DBI::dbAppendTable(
        con_out,
        tbl_name,
        tbl_main
      )
    }
    
    
    
    if (
      !is.null(tbl_add) &&
      nrow(tbl_add) > 0
    ) {
      
      DBI::dbAppendTable(
        con_out,
        tbl_name,
        tbl_add
      )
    }
    

    
    if (
      !is.null(tbl_main) &&
      !is.null(tbl_add)
    ) {
      
      cat(
        "Merged:",
        tbl_name,
        "| main:",
        nrow(tbl_main),
        "| add:",
        nrow(tbl_add),
        "| total:",
        nrow(tbl_main) + nrow(tbl_add),
        "\n"
      )
      
    } else if (!is.null(tbl_main)) {
      
      cat(
        "Only in MAIN:",
        tbl_name,
        "| rows:",
        nrow(tbl_main),
        "\n"
      )
      
    } else {
      
      cat(
        "Only in ADD:",
        tbl_name,
        "| rows:",
        nrow(tbl_add),
        "\n"
      )
    }
  }
  
  
  cat(
    "\n========================================\n",
    "Database merge completed.\n",
    "Output: ",
    output_db,
    "\n",
    "ADD RT conversion: ",
    ifelse(
      convert_add_rt,
      "seconds -> minutes",
      "none"
    ),
    "\n",
    "========================================\n",
    sep = ""
  )
}