#' @title convert a sql database to an mgf file
#'
#' @description sql_to_mgf() function first convert a sqlite database to a CompDb
#'  object and then to a spectra object after that we extract the spectra 
#'  variables needed to use with MsBackendMgf package and finally we get that
#'  mgf file of a given experiment with a spectrum_type in a sqlite database.
#'  
#' @param sqlite_path 'character()' the file path to a sql database.
#' 
#' @param expid 'numeric()' experiment number of the data to convert to  
#'        an mgf file.
#' 
#' @param spectrum_type 'character()' For flax ftms data it could be 
#'       "pseudo_MS2" or "individual".
#'       For Qtof flax data it could be "same_energy" 
#'       or "all_energies".
#'       
#' @param msLevel 'numeric()' msLevel number of the data to convert to  
#'        an mgf file.
#'               
#' @param output_file 'character()' the output directory of the final mgf file.
#' 
#' @return an mgf file with the data from the sqlite database.
#' 
#' @import Spectra
#' @import CompoundDb
#' @import MsBackendMgf
#' @import DBI
#' @import RSQLite 
#' 
#' @author Ahlam Mentag
#' 
#' @export                
sql_to_mgf <- function(sqlite_path,
                       expid,
                       msLevel = NULL,
                       spectrum_type = NULL,
                       output_file) {

  
  ## GNPS helper functions
  source(
    "https://raw.githubusercontent.com/jorainer/xcms-gnps-tools/master/customFunctions.R"
  )
  
  ##  Create temporary database
  tmp_db <- tempfile(fileext = ".sqlite")
  
  # Remove temporary database when function exits
  on.exit({
    if (file.exists(tmp_db)) unlink(tmp_db)
  }, add = TRUE)
  
  # Connect to original and temporary databases
  con_src <- dbConnect(SQLite(), sqlite_path)
  con_tmp <- dbConnect(SQLite(), tmp_db)
  
  # Ensure connections are closed on error
  on.exit({
    if (dbIsValid(con_src)) dbDisconnect(con_src)
    if (dbIsValid(con_tmp)) dbDisconnect(con_tmp)
  }, add = TRUE)
  
  ##  Keep only required tables
  required_tables <- c(
    "metadata",
    "ms_compound",
    "msms_spectrum",
    "msms_spectrum_peak"
  )
  
  # Check that all required tables exist
  missing_tables <- setdiff(
    required_tables,
    dbListTables(con_src)
  )
  
  if (length(missing_tables) > 0) {
    stop(
      "Missing required tables: ",
      paste(missing_tables, collapse = ", ")
    )
  }
  
  # Copy tables with their original SQL schemas
  for (tbl in required_tables) {
    
    schema <- dbGetQuery(
      con_src,
      paste0(
        "SELECT sql FROM sqlite_master WHERE type = 'table' ",
        "AND name = ",
        as.character(dbQuoteString(con_src, tbl))
      )
    )$sql
    
    dbExecute(con_tmp, schema)
    
    # Copy table content
    data <- dbReadTable(con_src, tbl)
    
    dbAppendTable(con_tmp, tbl, data)
  }
  
  # Close connections before loading with CompoundDb
  dbDisconnect(con_src)
  dbDisconnect(con_tmp)
  
  ##  Load temporary database
  cdb <- CompDb(tmp_db)
  
  ##  Convert database to Spectra
  sps <- Spectra(cdb)
  
  ##  Filtering
  sps <- sps[sps$expid %in% expid]
  
  if (!is.null(msLevel)) {
    sps <- sps[sps$msLevel %in% msLevel]
  }
  
  if (!is.null(spectrum_type)) {
    sps <- sps[sps$spectrum_type %in% spectrum_type]
  }
  
  if (length(sps) == 0) {
    stop("No spectra left after filtering")
  }
  
  ## 6. Keep variables needed for export
  sps <- selectSpectraVariables(
    sps,
    c(
      "spectrum_id",
      "compound_id",
      "compound_accession",
      "precursorMz",
      "rtime",
      "nodename",
      "dataOrigin",
      "msLevel",
      "precScanNum"
    )
  )
  
  ## 7. Nodename becomes GNPS feature ID
  sps$feature_id <- sps$nodename
  
  ## 8. Save metadata before formatting
  spectrum_ids <- sps$spectrum_id
  compound_ids <- sps$compound_id
  compound_accessions <- sps$compound_accession
  ms_levels <- sps$msLevel
  prec_scan_nums <- sps$precScanNum
  data_origins <- sps$dataOrigin
  
  ## Format spectra for GNPS
  sps_gnps <- formatSpectraForGNPS(sps)
  
  ## Create unique scan identifiers
  sps_gnps$scan_accession <- paste0(
    compound_accessions,
    "_",
    spectrum_ids
  )
  
  ##  Additional MGF metadata
  sps_gnps$COMPOUND_ID <- compound_ids
  sps_gnps$COMPOUND_ACCESSION <- compound_accessions
  sps_gnps$MSLEVEL <- ms_levels
  sps_gnps$PRECSANNUM <- prec_scan_nums
  sps_gnps$DATAORIGIN <- data_origins
  
  ##  MGF mapping
  mgf_mapping <- spectraVariableMapping(MsBackendMgf())
  
  # Remove acquisitionNum -> SCANS mapping
  mgf_mapping <- mgf_mapping[
    names(mgf_mapping) != "acquisitionNum"
  ]
  
  # Map character identifier to SCANS
  mgf_mapping["scan_accession"] <- "SCANS"
  
  ##Export MGF
  export(
    sps_gnps,
    backend = MsBackendMgf(),
    file = output_file,
    mapping = mgf_mapping
  )
  
  message(
    "MGF successfully exported: ",
    output_file,
    " (", length(sps_gnps), " spectra)"
  )
  
  invisible(TRUE)
}