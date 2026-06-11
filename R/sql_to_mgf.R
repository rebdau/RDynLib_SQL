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
#'       "not_assembled" or "assembled".
#'       For Qtof flax data it could be "intersect_single_energy" 
#'       or "merged_MSMS_all_energies".
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
  
  # GNPS helper functions
  source("https://raw.githubusercontent.com/jorainer/xcms-gnps-tools/master/customFunctions.R")
  
  ## Load database as Spectra
  cdb <- CompDb(sqlite_path)
  
  cdb <- addJoinDefinition(
    cdb,
    table_a  = "ms_compound",
    table_b  = "feature_matrix",
    column_a = "nodename",
    column_b = "nodename"
  )
  
  sps <- Spectra(cdb)
  
  ## Filtering
  sps <- sps[sps$expid %in% expid]
  
  if (!is.null(msLevel)) {
    sps <- sps[sps$msLevel %in% msLevel]
  }
  
  if (!is.null(spectrum_type)) {
    sps <- sps[sps$spectrum_type %in% spectrum_type]
  }
  
  if (length(sps) == 0)
    stop("No spectra left after filtering")
  
  # Restrict variables to GNPS relevant ones
  sps <- selectSpectraVariables(
    sps,
    c( "compound_id",
       "precursorMz",
       "rtime",
       "nodename",
       "dataOrigin",
       "msLevel",
       "precScanNum"
    )
  )
  
  
  ## Map nodename to feature_id (GNPS requirement)
  sps$feature_id <- sps$nodename
  
  sps <- selectSpectraVariables(
    sps,
    c( "compound_id",
       "precursorMz",
       "rtime",
       "feature_id",
       "dataOrigin",
       "msLevel",
       "precScanNum"
    )
  )
  
  
  ## Format spectra for GNPS
  sps_gnps <- formatSpectraForGNPS(sps)
  
  sps_gnps$COMPOUND_ID <- sps$compound_id
  sps_gnps$MSLEVEL     <- sps$msLevel
  sps_gnps$PRECSANNUM <- sps$precScanNum
  sps_gnps$DATAORIGIN  <- sps$dataOrigin
  
  
  ## Export MGF using MsBackendMgf
  export(
    sps_gnps,
    backend = MsBackendMgf(),
    file    = output_file
  )
  
  invisible(TRUE)
}



