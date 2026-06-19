#' @title Filter and create a SQL subdatabase from one or more experiment IDs
#'
#' @description
#' Creates a new SQLite database containing only records associated with the
#' specified experiment IDs.
#'
#' @param expid Numeric vector of experiment IDs to keep.
#' @param input_db Character. Path to the input SQLite database.
#' @param output_db Character. Path to the output SQLite database.
#' @param overwrite Logical. Whether to overwrite an existing output database.
#'
#' @return Invisibly returns the path to the output database.
#'
#' @author Ahlam Mentag
#'
#' @import DBI
#' @import RSQLite
#'
#' @export
filterDynLibDB <- function(expid,
                           input_db,
                           output_db,
                           overwrite = TRUE) {
  
  if (!file.exists(input_db))
    stop("Input database does not exist.")
  
  if (!length(expid))
    stop("'expid' must contain at least one experiment ID.")
  
  expid <- unique(as.integer(expid))
  
  if (file.exists(output_db)) {
    if (overwrite) {
      file.remove(output_db)
    } else {
      stop("Output database already exists.")
    }
  }
  
  ## Copy the entire database first
  file.copy(input_db, output_db, overwrite = overwrite)
  
  con <- DBI::dbConnect(RSQLite::SQLite(), output_db)
  
  on.exit(DBI::dbDisconnect(con))
  
  expid_sql <- paste(expid, collapse = ",")
  
  ## Verify requested experiments exist
  existing <- DBI::dbGetQuery(
    con,
    sprintf(
      "SELECT expid FROM experiment WHERE expid IN (%s)",
      expid_sql
    )
  )
  
  if (!nrow(existing))
    stop("None of the provided expid values were found.")
  
  ## Get compound IDs to keep
  compounds_keep <- DBI::dbGetQuery(
    con,
    sprintf(
      "SELECT DISTINCT compound_id
             FROM ms_compound
             WHERE expid IN (%s)",
      expid_sql
    )
  )
  
  compound_ids <- compounds_keep$compound_id
  
  ## Get spectrum IDs to keep
  if (length(compound_ids)) {
    
    compound_ids_sql <- paste(compound_ids, collapse = ",")
    
    spectra_keep <- DBI::dbGetQuery(
      con,
      sprintf(
        "SELECT DISTINCT spectrum_id
                 FROM msms_spectrum
                 WHERE compound_id IN (%s)",
        compound_ids_sql
      )
    )
    
    spectrum_ids <- spectra_keep$spectrum_id
    
  } else {
    
    spectrum_ids <- integer()
  }
  
  ## Remove unwanted experiment rows
  DBI::dbExecute(
    con,
    sprintf(
      "DELETE FROM experiment
             WHERE expid NOT IN (%s)",
      expid_sql
    )
  )
  
  ## Remove unwanted compounds
  DBI::dbExecute(
    con,
    sprintf(
      "DELETE FROM ms_compound
             WHERE expid NOT IN (%s)",
      expid_sql
    )
  )
  
  ## Remove unwanted feature matrix rows
  fm_table <- if ("feature_matrix" %in% DBI::dbListTables(con))
    "feature_matrix"
  else
    "Feature_matrix"
  
  DBI::dbExecute(
    con,
    sprintf(
      "DELETE FROM %s
             WHERE expid NOT IN (%s)",
      fm_table,
      expid_sql
    )
  )
  
  ## Remove unwanted spectra
  if (length(compound_ids)) {
    
    DBI::dbExecute(
      con,
      sprintf(
        "DELETE FROM msms_spectrum
                 WHERE compound_id NOT IN (%s)",
        compound_ids_sql
      )
    )
    
  } else {
    
    DBI::dbExecute(
      con,
      "DELETE FROM msms_spectrum"
    )
  }
  
  ## Remove unwanted peaks
  if (length(spectrum_ids)) {
    
    spectrum_ids_sql <- paste(spectrum_ids, collapse = ",")
    
    DBI::dbExecute(
      con,
      sprintf(
        "DELETE FROM msms_spectrum_peak
                 WHERE Spectrum_id NOT IN (%s)",
        spectrum_ids_sql
      )
    )
    
  } else {
    
    DBI::dbExecute(
      con,
      "DELETE FROM msms_spectrum_peak"
    )
  }
  
  invisible(output_db)
}