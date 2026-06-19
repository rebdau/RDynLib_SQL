#' @title Aligning compounds from the same instrument 
#'
#' @description
#'
#' `Aligning_FTn_FTp_SQL()` can align FTneg-FTpos or QTOFneg-QTOFpos, the
#'  reference experiment is always the one with the negative ionization mode 
#'  and adds the aligned compounds to a given text file.
#'
#' @details
#'
#' the function reads the both sql databases FTMS or QTOF, and extracts
#' their modes, create an output file to store the alignment results in a .txt
#' file if it's already given by the user otherwise we create a new one.
#' Then calls a set of defined functions :
#'
#' - [Aligning_General_SQL()]: performs the local alignment, and adjust the mz
#'   values if the alignment is between two different polarities.
#'
#' - `matchFTSyn_SQL()`: removes the matches found by the previous function
#'   if they contain fewer than `minIon` matched peaks.
#'
#' - `RemoveOutliers()`: calculate the standard deviation in the rt difference
#'   found between each matched feature and remove all rows outside a range.
#'
#' - `RegressionPie_LCalign()`: apply a regression model  to reduce the
#'   influence of false matches and outliers.
#'
#' - `FillAssocFTnQTOFn_SQL()`: store the final aligned compounds to a .txt
#'   file.
#'
#' @param FTn_path `character(1)` the path to the reference database with 
#'        negative polarity.
#'        
#' @param FTp_path `character(1)` the path of the SQL database to align with.
#'
#' @param FTn_expnr `numeric(1)` the reference experiment to align.
#'
#' @param FTp_expnr `numeric(1)` the experiment to align with the previous
#'        reference experiment.
#'
#' @param Assoc `character(1)` the file that contains some previous alignment
#'        results done by the user. If `NULL` we create it in the function.
#'
#' @param err `numeric(1)` is the m/z tolerance window used to decide which QTOF
#'        peaks are potential isomeric matches for a given FT peak.
#'
#' @param t.ini `numeric(1)` used in the `Aligning_General_SQL()` to define the
#'        number of neighboring peaks.
#'
#' @param lc.err `numeric(1)` the retention-time alignment tolerance.
#'
#' @param rng `numeric(1)` used in `RemoveOutliers()` refers to the number of
#'        standard  deviations used to define the outlier cutoff.
#'
#' @param minIon `numeric(1)`  
#'    The minimum matching ions between spectra of peak pairs, which are found
#'    after applying the regression model, and which are the final matching 
#'    written in the Assoc table. (default `0.01`).
#'
#' @param minPeaks Numeric. Minimum MS2 matching ratio required.
#'
#' @param startpoint `numeric(1)` used in `RegressionPie_LCalign()` refers to
#'        the minimum retention time for the first knot in piecewise regression.
#'
#' @param save_assoc `logical(1)` If true we add the alignment result to the
#'        Assoc text file.
#'
#' @return A `data.frame` containing the matched FTMSn-FTMSp or QTOFn-QTOFp 
#'         compounds with columns:
#' - `ref_compid`: FTMSneg or QTOFneg compound ID
#' - `target_compid`: FTMSpos or QTOFpos compound ID
#' - `ref_database`: filename of the reference database
#' - `target_database`: filename of the database to align with.
#' This table is updated with new associations found by the alignment procedure.
#'
#'
#' @importFrom DBI dbGetQuery
#'
#' @importFrom DBI dbConnect
#'
#' @importFrom DBI dbDisconnect
#'
#' @importFrom RSQLite SQLite
#'
#' @importFrom utils read.table
#'
#' @importFrom utils write.table
#'
#' @author Ahlam Mentag
#'
#' @export
Aligning_FTn_FTp_SQL <- function(
    FTn_path, FTp_path, FTn_expnr, FTp_expnr,
    Assoc = NULL,
    err = 0.02,
    t.ini = 5,
    lc.err = 1,
    rng = 2,
    minIon = 0.01,
    minPeaks = 0.6,
    startpoint = 1,
    save_assoc = FALSE
) {
  
  # Load / init Assoc
  if (is.character(Assoc)) {
    
    assoc_file <- Assoc
    
    if (file.exists(Assoc)) {
      Assoc_df <- read.table(Assoc, header = TRUE, sep = "\t")
    } else {
      warning("File ", assoc_file, " not found. Proceeding without.")
      Assoc_df <- data.frame(
        ref_compid = integer(0),
        target_compid = integer(0),
        ref_database = character(0),
        target_database = character(0)
      )
    }
    
  } else if (is.data.frame(Assoc)) {
    Assoc_df <- Assoc
    assoc_file <- NULL
    
  } else {
    Assoc_df <- data.frame(
      ref_compid = integer(0),
      target_compid = integer(0),
      ref_database = character(0),
      target_database = character(0)
    )
    assoc_file <- NULL
  }
  
  # DB connection
  FTn_con <- dbConnect(RSQLite::SQLite(), FTn_path)
  FTp_con <- dbConnect(RSQLite::SQLite(), FTp_path)
  
  on.exit({
    dbDisconnect(FTn_con)
    dbDisconnect(FTp_con)
  }, add = TRUE)
  
  
  # DB names 
  ftn_db_name <- basename(FTn_path)
  ftp_db_name <- basename(FTp_path)
  
  # Experiment mode
  ftn_mode <- dbGetQuery(
    FTn_con,
    sprintf("SELECT mode FROM experiment WHERE expid = %d", FTn_expnr)
  )$mode
  
  ftp_mode <- dbGetQuery(
    FTp_con,
    sprintf("SELECT mode FROM experiment WHERE expid = %d", FTp_expnr)
  )$mode
  
  if (!length(ftn_mode)) stop("FTn experiment not found")
  if (!length(ftp_mode)) stop("FTp experiment not found")
  
  polarity_ftn <- ifelse(ftn_mode == "neg", 0, 1)
  polarity_ftp <- ifelse(ftp_mode == "neg", 0, 1)
  
  
  # Alignment
  LCal <- Aligning_General_SQL(
    FT_con = FTn_con,
    QTOF_con = FTp_con,
    FT_expnr = FTn_expnr,
    QTOF_expnr = FTp_expnr,
    err,
    t.ini
  )
  
  LCal <- matchNegPos_SQL(
    LCal = LCal,
    FTn_con = FTn_con,
    FTp_con = FTp_con,
    polarity_ftn = polarity_ftn,
    polarity_ftp = polarity_ftp,
    minIon = minIon
  )
  
  LCal <- RemoveOutliers(LCal, rng)
  
  rg <- RegressionPie_LCalign_SQL(LCal, startpoint)
  PlotPie_LCalign(LCal, rg)
  
  
  # Fill associations
  Assoc_df <- FillAssocFTnFTpn_SQL(
    FTn_con = FTn_con,
    FTp_con = FTp_con,
    Assoc = Assoc_df,
    FTn_expnr = FTn_expnr,
    FTp_expnr = FTp_expnr,
    cutoff = 1,
    rg = rg,
    lc.err = lc.err,
    err = err,
    minIon = minIon,
    polarity_ftn = polarity_ftn,
    polarity_ftp = polarity_ftp,
    FTn_path = ftn_db_name,
    FTp_path = ftp_db_name
  )
  
  
  if (save_assoc && !is.null(assoc_file)) {
    write.table(
      Assoc_df,
      file = assoc_file,
      sep = "\t",
      row.names = FALSE,
      quote = FALSE
    )
  }
  
  return(Assoc_df)
}