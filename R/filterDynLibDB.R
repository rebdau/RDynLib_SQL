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
  
  con_in <- DBI::dbConnect(RSQLite::SQLite(), input_db)
  con_out <- DBI::dbConnect(RSQLite::SQLite(), output_db)
  
  on.exit({
    DBI::dbDisconnect(con_in)
    DBI::dbDisconnect(con_out)
  })
  
  expid_sql <- paste(expid, collapse = ",")
  
  ## experiment
  experiment <- DBI::dbGetQuery(
    con_in,
    sprintf(
      "SELECT * FROM experiment WHERE expid IN (%s)",
      expid_sql
    )
  )
  
  if (!nrow(experiment))
    stop("None of the provided expid values were found.")
  
  DBI::dbWriteTable(
    con_out,
    "experiment",
    experiment,
    row.names = FALSE
  )
  
  ## ms_compound
  compounds <- DBI::dbGetQuery(
    con_in,
    sprintf(
      "SELECT * FROM ms_compound WHERE expid IN (%s)",
      expid_sql
    )
  )
  
  if (!nrow(compounds)) {
    compounds <- DBI::dbGetQuery(
      con_in,
      "SELECT * FROM ms_compound WHERE 1 = 0"
    )
  }
  
  DBI::dbWriteTable(
    con_out,
    "ms_compound",
    compounds,
    row.names = FALSE
  )
  
  compound_ids <- unique(compounds$compound_id)
  
  ## Feature_matrix
  feature_matrix <- DBI::dbGetQuery(
    con_in,
    sprintf(
      "SELECT * FROM Feature_matrix WHERE expid IN (%s)",
      expid_sql
    )
  )
  
  if (!nrow(feature_matrix)) {
    feature_matrix <- DBI::dbGetQuery(
      con_in,
      "SELECT * FROM Feature_matrix WHERE 1 = 0"
    )
  }
  
  DBI::dbWriteTable(
    con_out,
    "Feature_matrix",
    feature_matrix,
    row.names = FALSE
  )
  
  ## msms_spectrum
  if (length(compound_ids)) {
    
    compound_ids_sql <- paste(compound_ids, collapse = ",")
    
    spectra <- DBI::dbGetQuery(
      con_in,
      sprintf(
        paste0(
          "SELECT * FROM msms_spectrum ",
          "WHERE compound_id IN (%s)"
        ),
        compound_ids_sql
      )
    )
    
  } else {
    
    spectra <- DBI::dbGetQuery(
      con_in,
      "SELECT * FROM msms_spectrum WHERE 1 = 0"
    )
  }
  
  DBI::dbWriteTable(
    con_out,
    "msms_spectrum",
    spectra,
    row.names = FALSE
  )
  
  ## msms_spectrum_peak
  if (nrow(spectra)) {
    
    spectrum_ids_sql <- paste(
      unique(spectra$spectrum_id),
      collapse = ","
    )
    
    peaks <- DBI::dbGetQuery(
      con_in,
      sprintf(
        paste0(
          "SELECT * FROM msms_spectrum_peak ",
          "WHERE Spectrum_id IN (%s)"
        ),
        spectrum_ids_sql
      )
    )
    
  } else {
    
    peaks <- DBI::dbGetQuery(
      con_in,
      "SELECT * FROM msms_spectrum_peak WHERE 1 = 0"
    )
  }
  
  DBI::dbWriteTable(
    con_out,
    "msms_spectrum_peak",
    peaks,
    row.names = FALSE
  )
  
  invisible(output_db)
}