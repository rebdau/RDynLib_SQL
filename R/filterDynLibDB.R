#' @title Filter and Create an sql subdatabase of a given experiment id  
#'
#' 
#' @return a new sql database containing just the given experiment id
#' 
#' @import DBI
#' @import RSQLite
#' 
#' @param expid 'numeric(0) experiment id to keep. 
#' @param input_db 'character' file path to the global sql database to filter
#' @param output_db 'character' output file path to the filtered sql database
#' @author Ahlam Mentag
#' 
#' @export
filterDynLibDB <- function(expid,
                           input_db,
                           output_db,
                           overwrite = TRUE) {
  
  if (!file.exists(input_db))
    stop("Input database does not exist.")
  
  if (file.exists(output_db)) {
    if (overwrite) {
      file.remove(output_db)
    } else {
      stop("Output database already exists.")
    }
  }
  
  con_in <- dbConnect(SQLite(), input_db)
  con_out <- dbConnect(SQLite(), output_db)
  
  on.exit({
    dbDisconnect(con_in)
    dbDisconnect(con_out)
  })
  

  ## experiment
  experiment <- dbGetQuery(
    con_in,
    sprintf(
      "SELECT * FROM experiment WHERE expid = '%s'",
      expid
    )
  )
  
  if (nrow(experiment) == 0)
    stop("expid not found.")
  
  dbWriteTable(con_out, "experiment", experiment)
  
  
  ## ms_compound
  compounds <- dbGetQuery(
    con_in,
    sprintf(
      "SELECT * FROM ms_compound WHERE expid = '%s'",
      expid
    )
  )
  
  dbWriteTable(con_out, "ms_compound", compounds)
  
  compound_ids <- compounds$compound_id
  

  ## Feature_matrix
  feature_matrix <- dbGetQuery(
    con_in,
    sprintf(
      "SELECT * FROM Feature_matrix WHERE expid = '%s'",
      expid
    )
  )
  
  dbWriteTable(con_out, "Feature_matrix", feature_matrix)
  

  ## msms_spectrum
  if (length(compound_ids) > 0) {
    
    compound_ids_sql <- paste(
      sprintf("'%s'", compound_ids),
      collapse = ","
    )
    
    spectra <- dbGetQuery(
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
    spectra <- data.frame()
  }
  
  dbWriteTable(con_out, "msms_spectrum", spectra)
  

  ## msms_spectrum_peak
  if (nrow(spectra) > 0) {
    
    spectrum_ids_sql <- paste(
      spectra$spectrum_id,
      collapse = ","
    )
    
    peaks <- dbGetQuery(
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
    peaks <- data.frame()
  }
  
  dbWriteTable(con_out, "msms_spectrum_peak", peaks)
  
  invisible(output_db)
}



