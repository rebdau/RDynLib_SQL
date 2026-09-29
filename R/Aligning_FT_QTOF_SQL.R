#' @title Aligning FTMS and QTOF experiments with the same/different polarities
#'
#' @description
#'
#' `Aligning_FT_QTOF_SQL()` aligns FTMS and QTOF experiments with the same
#' or different polarities, by calling a set of defined functions, and adds the
#' aligned compounds to a given text file.
#'
#' @details
#'
#' the function reads the both sql databases FTMS and QTOF, and extracts
#' their mode, create an output file to store the alignment results in a .txt
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
#' @param FT_path `character(1)` the path to the FTMS SQL database.
#'
#' @param QTOF_path `character(1)` the path to the QTOF SQL database.
#'
#' @param FT_expnr `numeric(1)` the FTMS reference experiment to align.
#'
#' @param QTOF_expnr `numeric(1)` the QTOF experiment to align with the previous
#'        FTMS experiment.
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
#' @param spectrum_type_FT `character(1)` Spectrum type to use for the FT-MS
#'        spectra. If `NULL` and the database does not contain a `spectrum_type`
#'        column, all spectra are used.
#'
#' @param spectrum_type_QTOF `character(1)` Spectrum type to use for the QTOF
#'        spectra. If `NULL` and the database does not contain a `spectrum_type`
#'        column, all spectra are used.
#'
#' @return A `data.frame` containing the matched FT–QTOF compounds with columns:
#' - `ref_compid`: FTMS compound ID
#' - `target_compid`: QTOF compound ID
#' - `ref_database`: filename of the FTMS database
#' - `target_database`: filename of the QTOF database.
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
Aligning_FT_QTOF_SQL <- function(
    FT_path,
    QTOF_path,
    FT_expnr,
    QTOF_expnr,
    Assoc = NULL,
    err = 0.02,
    t.ini = 5,
    lc.err = 1,
    rng = 2,
    minIon = 0.01,
    minPeaks = 0.6,
    startpoint = 1,
    save_assoc = FALSE,
    spectrum_type_FT = NULL,
    spectrum_type_QTOF = NULL
) {
  
  # Check experiment IDs
  if (length(FT_expnr) != 1L || is.na(FT_expnr)) {
    stop("`FT_expnr` must be a single experiment ID.")
  }
  
  if (length(QTOF_expnr) != 1L || is.na(QTOF_expnr)) {
    stop("`QTOF_expnr` must be a single experiment ID.")
  }
  
  
  # Load / initialize Assoc
  if (is.character(Assoc)) {
    
    assoc_file <- Assoc
    
    if (file.exists(Assoc)) {
      
      Assoc_df <- utils::read.table(
        Assoc,
        header = TRUE,
        sep = "\t",
        stringsAsFactors = FALSE
      )
      
    } else {
      
      warning(
        "File ",
        assoc_file,
        " not found. Proceeding without."
      )
      
      Assoc_df <- data.frame(
        ref_compid = integer(0),
        target_compid = integer(0),
        ref_database = character(0),
        target_database = character(0),
        stringsAsFactors = FALSE
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
      target_database = character(0),
      stringsAsFactors = FALSE
    )
    
    assoc_file <- NULL
  }
  
  
  # Database connections
  FT_con <- DBI::dbConnect(
    RSQLite::SQLite(),
    FT_path
  )
  
  QTOF_con <- DBI::dbConnect(
    RSQLite::SQLite(),
    QTOF_path
  )
  
  on.exit({
    DBI::dbDisconnect(FT_con)
    DBI::dbDisconnect(QTOF_con)
  }, add = TRUE)
  
  
  # Check spectrum_type columns
  has_type_FT <- "spectrum_type" %in%
    DBI::dbListFields(
      FT_con,
      "msms_spectrum"
    )
  
  has_type_QTOF <- "spectrum_type" %in%
    DBI::dbListFields(
      QTOF_con,
      "msms_spectrum"
    )
  
  
  # Check FT spectrum type
  if (has_type_FT) {
    
    available_types_FT <- DBI::dbGetQuery(
      FT_con,
      "
      SELECT DISTINCT spectrum_type
      FROM msms_spectrum
      WHERE spectrum_type IS NOT NULL
        AND spectrum_type != ''
      "
    )$spectrum_type
    
    
    if (is.null(spectrum_type_FT)) {
      
      stop(
        "`spectrum_type_FT` must be provided because ",
        "the FT database contains a `spectrum_type` column. ",
        "Available types: ",
        paste(available_types_FT, collapse = ", ")
      )
    }
    
    
    if (length(spectrum_type_FT) != 1L ||
        is.na(spectrum_type_FT) ||
        !nzchar(spectrum_type_FT)) {
      
      stop(
        "`spectrum_type_FT` must be a single non-empty character value."
      )
    }
    
    
    if (!spectrum_type_FT %in% available_types_FT) {
      
      stop(
        "spectrum_type_FT = '",
        spectrum_type_FT,
        "' not found. Available types: ",
        paste(available_types_FT, collapse = ", ")
      )
    }
  }
  
  
  # Check QTOF spectrum type
  if (has_type_QTOF) {
    
    available_types_QTOF <- DBI::dbGetQuery(
      QTOF_con,
      "
      SELECT DISTINCT spectrum_type
      FROM msms_spectrum
      WHERE spectrum_type IS NOT NULL
        AND spectrum_type != ''
      "
    )$spectrum_type
    
    
    if (is.null(spectrum_type_QTOF)) {
      
      stop(
        "`spectrum_type_QTOF` must be provided because ",
        "the QTOF database contains a `spectrum_type` column. ",
        "Available types: ",
        paste(available_types_QTOF, collapse = ", ")
      )
    }
    
    
    if (length(spectrum_type_QTOF) != 1L ||
        is.na(spectrum_type_QTOF) ||
        !nzchar(spectrum_type_QTOF)) {
      
      stop(
        "`spectrum_type_QTOF` must be a single non-empty character value."
      )
    }
    
    
    if (!spectrum_type_QTOF %in% available_types_QTOF) {
      
      stop(
        "spectrum_type_QTOF = '",
        spectrum_type_QTOF,
        "' not found. Available types: ",
        paste(available_types_QTOF, collapse = ", ")
      )
    }
  }
  
  
  # Detect polarity
  ft_mode <- DBI::dbGetQuery(
    FT_con,
    sprintf(
      "SELECT mode
       FROM experiment
       WHERE expid = %d",
      FT_expnr
    )
  )$mode
  
  
  qtof_mode <- DBI::dbGetQuery(
    QTOF_con,
    sprintf(
      "SELECT mode
       FROM experiment
       WHERE expid = %d",
      QTOF_expnr
    )
  )$mode
  
  
  if (!length(ft_mode)) {
    stop(
      "Experiment '",
      FT_expnr,
      "' not found in ",
      FT_path
    )
  }
  
  if (!length(qtof_mode)) {
    stop(
      "Experiment '",
      QTOF_expnr,
      "' not found in ",
      QTOF_path
    )
  }
  
  
  ft_mode <- tolower(ft_mode[1])
  qtof_mode <- tolower(qtof_mode[1])
  
  negative_values <- c(
    "neg",
    "negative",
    "-",
    "negativ"
  )
  
  polarity_ft <- ifelse(
    ft_mode %in% negative_values,
    0,
    1
  )
  
  polarity_qtof <- ifelse(
    qtof_mode %in% negative_values,
    0,
    1
  )
  
  
  # Initial LC alignment
  LCal <- Aligning_General_SQL(
    FT_con = FT_con,
    QTOF_con = QTOF_con,
    FT_expnr = FT_expnr,
    QTOF_expnr = QTOF_expnr,
    err = err,
    t.ini = t.ini
  )
  
  
  if (is.null(LCal) || nrow(LCal) == 0) {
    warning("No initial alignment candidates found.")
    return(Assoc_df)
  }
  
  
  # MS2 matching
  LCal <- matchFTSyn_SQL(
    LCal = LCal,
    FT_con = FT_con,
    QTOF_con = QTOF_con,
    FT_expnr = FT_expnr,
    QTOF_expnr = QTOF_expnr,
    polarity_ft = polarity_ft,
    polarity_qtof = polarity_qtof,
    minPeaks = minPeaks,
    spectrum_type_FT = spectrum_type_FT,
    spectrum_type_QTOF = spectrum_type_QTOF
  )
  
  
  if (is.null(LCal) || nrow(LCal) == 0) {
    warning(
      "No candidates remained after MS2 matching."
    )
    return(Assoc_df)
  }
  
  
  # Remove RT outliers
  LCal <- RemoveOutliers(
    LCal,
    rng
  )
  
  
  if (is.null(LCal) || nrow(LCal) == 0) {
    warning(
      "No candidates remained after outlier removal."
    )
    return(Assoc_df)
  }
  
  
  # Regression
  rg <- RegressionPie_LCalign_SQL(
    LCal,
    startpoint
  )
  
  
  if (interactive() && capabilities("cairo")) {
    try(
      PlotPie_LCalign(
        LCal,
        rg
      ),
      silent = TRUE
    )
  }
  
  
  # Fill association table
  Assoc_df <- FillAssocFTnQTOFn_SQL(
    FT_con = FT_con,
    QTOF_con = QTOF_con,
    Assoc = Assoc_df,
    FT_expnr = FT_expnr,
    QTOF_expnr = QTOF_expnr,
    cutoff = 1,
    rg = rg,
    lc.err = lc.err,
    err = err,
    minIon = minIon,
    polarity_ft = polarity_ft,
    polarity_qtof = polarity_qtof,
    FT_path = FT_path,
    QTOF_path = QTOF_path,
    spectrum_type_FT = spectrum_type_FT,
    spectrum_type_QTOF = spectrum_type_QTOF
  )
  
  
  # Save associations
  if (save_assoc && !is.null(assoc_file)) {
    
    utils::write.table(
      Assoc_df,
      file = assoc_file,
      sep = "\t",
      row.names = FALSE,
      quote = FALSE
    )
  }
  
  
  return(Assoc_df)
}