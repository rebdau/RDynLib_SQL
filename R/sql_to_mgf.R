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
  
  ## Load database
  cdb <- CompDb(sqlite_path)
  
  cdb <- addJoinDefinition(
    cdb,
    table_a  = "ms_compound",
    table_b  = "feature_matrix",
    column_a = "nodename",
    column_b = "nodename"
  )
  
  ## Convert database to Spectra
  sps <- Spectra(cdb)
  
  ## Filtering
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
  
  ## Keep variables needed for export
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
  
  ## nodename becomes GNPS feature ID
  sps$feature_id <- sps$nodename
  
  ## Save metadata before formatSpectraForGNPS()
  spectrum_ids <- sps$spectrum_id
  compound_ids <- sps$compound_id
  compound_accessions <- sps$compound_accession
  ms_levels <- sps$msLevel
  prec_scan_nums <- sps$precScanNum
  data_origins <- sps$dataOrigin
  
  ## Format spectra for GNPS
  sps_gnps <- formatSpectraForGNPS(sps)
  
  
  sps_gnps$scan_accession <- paste0(
    compound_accessions,
    "_",
    spectrum_ids
  )
  
  ## Additional MGF metadata
  sps_gnps$COMPOUND_ID <- compound_ids
  sps_gnps$COMPOUND_ACCESSION <- compound_accessions
  sps_gnps$MSLEVEL <- ms_levels
  sps_gnps$PRECSANNUM <- prec_scan_nums
  sps_gnps$DATAORIGIN <- data_origins
  
  ## MGF mapping
  
  mgf_mapping <- spectraVariableMapping(MsBackendMgf())
  
  ## Remove the normal acquisitionNum -> SCANS mapping
  mgf_mapping <- mgf_mapping[
    !(names(mgf_mapping) == "acquisitionNum")
  ]
  
  ## Map our character identifier to SCANS
  mgf_mapping["scan_accession"] <- "SCANS"
  
  ## Export MGF
  export(
    sps_gnps,
    backend = MsBackendMgf(),
    file = output_file,
    mapping = mgf_mapping
  )
  
  invisible(TRUE)
}