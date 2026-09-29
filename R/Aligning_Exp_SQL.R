#' @title Align two experiments from the same SQL database
#'
#' @description
#' Aligns two experiments with the same polarity stored in the same DynLib
#' SQLite database using m/z, MS2 similarity, and retention-time regression.
#'
#' @param db_path `character(1)` Path to the SQLite database.
#' @param ref_expnr `numeric(1)` Reference experiment ID.
#' @param target_expnr `numeric(1)` Target experiment ID.
#' @param Assoc Association data frame or path to an existing association file.
#' @param err `numeric(1)` m/z tolerance used for the initial alignment.
#' @param t.ini `numeric(1)` Number of neighboring peaks used for alignment.
#' @param lc.err `numeric(1)` Retention-time alignment tolerance.
#' @param rng `numeric(1)` Number of standard deviations used to remove outliers.
#' @param thr1 `numeric(1)` Minimum fraction of matching MS2 ions.
#' @param thr2 `numeric(1)` Minimum MS2 dot-product similarity.
#' @param IntThres `numeric(1)` Minimum MS2 peak intensity.
#' @param startpoint `numeric(1)` Minimum retention time for piecewise regression.
#' @param tol `numeric(1)` m/z tolerance for matching MS2 fragment ions.
#' @param save_assoc `logical(1)` Whether to save the associations to `Assoc`.
#' @param spectrum_type_ref `character(1)` Spectrum type to use for the experiment
#'        to align. If `NULL` and the database does not contain a `spectrum_type`
#'        column, all spectra are used.
#'
#' @param spectrum_type_target `character(1)` Spectrum type to use for the 
#'        experiment to align with. If `NULL` and the database does not contain
#'        a `spectrum_type` column, all spectra are used.
#'        
#' @return A `data.frame` containing the aligned compound IDs and database names.
#'
#' @importFrom DBI dbConnect dbDisconnect dbGetQuery
#' @importFrom RSQLite SQLite
#' @importFrom utils read.table write.table
#'
#' @author Ahlam Mentag
#'
#' @export
Aligning_Exp_SQL <- function(
    db_path,
    ref_expnr,
    target_expnr,
    Assoc = NULL,
    err = 0.01,
    t.ini = 5,
    lc.err = 1,
    rng = 2,
    thr1 = 1,
    thr2 = 0.9,
    IntThres = 100,
    startpoint = 1,
    tol = 0.01,
    save_assoc = FALSE,
    spectrum_type_ref = NULL,
    spectrum_type_target = NULL
) {
  
  # Check experiment IDs
  if (length(ref_expnr) != 1L || is.na(ref_expnr)) {
    stop("`ref_expnr` must be a single experiment ID.")
  }
  
  if (length(target_expnr) != 1L || is.na(target_expnr)) {
    stop("`target_expnr` must be a single experiment ID.")
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
        " not found. A new association table will be created."
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
  
  
  # Database connection
  con <- DBI::dbConnect(
    RSQLite::SQLite(),
    db_path
  )
  
  on.exit(
    DBI::dbDisconnect(con),
    add = TRUE
  )
  
  
  # Check spectrum_type column
  has_spectrum_type <- "spectrum_type" %in%
    DBI::dbListFields(
      con,
      "msms_spectrum"
    )
  
  
  # Check requested spectrum types
  if (has_spectrum_type) {
    
    available_types <- DBI::dbGetQuery(
      con,
      "
      SELECT DISTINCT spectrum_type
      FROM msms_spectrum
      WHERE spectrum_type IS NOT NULL
        AND spectrum_type != ''
      "
    )$spectrum_type
    
    
    if (is.null(spectrum_type_ref)) {
      
      stop(
        "`spectrum_type_ref` must be provided because ",
        "the database contains a `spectrum_type` column. ",
        "Available types: ",
        paste(available_types, collapse = ", ")
      )
    }
    
    
    if (is.null(spectrum_type_target)) {
      
      stop(
        "`spectrum_type_target` must be provided because ",
        "the database contains a `spectrum_type` column. ",
        "Available types: ",
        paste(available_types, collapse = ", ")
      )
    }
    
    
    if (length(spectrum_type_ref) != 1L ||
        is.na(spectrum_type_ref) ||
        !nzchar(spectrum_type_ref)) {
      
      stop(
        "`spectrum_type_ref` must be a single non-empty character value."
      )
    }
    
    
    if (length(spectrum_type_target) != 1L ||
        is.na(spectrum_type_target) ||
        !nzchar(spectrum_type_target)) {
      
      stop(
        "`spectrum_type_target` must be a single non-empty character value."
      )
    }
    
    
    if (!spectrum_type_ref %in% available_types) {
      
      stop(
        "spectrum_type_ref = '",
        spectrum_type_ref,
        "' not found. Available types: ",
        paste(available_types, collapse = ", ")
      )
    }
    
    
    if (!spectrum_type_target %in% available_types) {
      
      stop(
        "spectrum_type_target = '",
        spectrum_type_target,
        "' not found. Available types: ",
        paste(available_types, collapse = ", ")
      )
    }
  }
  
  
  # Check experiments
  experiments <- DBI::dbGetQuery(
    con,
    sprintf(
      "
      SELECT expid, mode
      FROM experiment
      WHERE expid IN (%d, %d)
      ",
      ref_expnr,
      target_expnr
    )
  )
  
  
  if (!ref_expnr %in% experiments$expid) {
    stop(
      "Reference experiment ",
      ref_expnr,
      " not found."
    )
  }
  
  
  if (!target_expnr %in% experiments$expid) {
    stop(
      "Target experiment ",
      target_expnr,
      " not found."
    )
  }
  
  
  ref_mode <- experiments$mode[
    experiments$expid == ref_expnr
  ][1]
  
  target_mode <- experiments$mode[
    experiments$expid == target_expnr
  ][1]
  
  
  # This function is for same-polarity experiments
  if (tolower(ref_mode) != tolower(target_mode)) {
    
    stop(
      "Aligning_Exp_SQL() is intended for experiments ",
      "with the same polarity. Found: ",
      ref_mode,
      " and ",
      target_mode,
      "."
    )
  }
  
  
  # Initial LC alignment
  LCal <- Aligning_General_SQL(
    FT_con = con,
    QTOF_con = con,
    FT_expnr = ref_expnr,
    QTOF_expnr = target_expnr,
    err = err,
    t.ini = t.ini
  )
  
  
  if (is.null(LCal) || nrow(LCal) == 0) {
    warning("No initial alignment candidates found.")
    return(Assoc_df)
  }
  
  
  # MS2 matching
  LCal <- match_FTexp_SQL(
    LCal = LCal,
    con = con,
    expnr1 = ref_expnr,
    expnr2 = target_expnr,
    thr1 = thr1,
    thr2 = thr2,
    IntThres = IntThres,
    tol = tol,
    spectrum_type_ref = spectrum_type_ref,
    spectrum_type_target = spectrum_type_target
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
  
  
  # Piecewise regression
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
  
  
  # Final associations
  Assoc_df <- FillAssocExp_SQL(
    con = con,
    Assoc = Assoc_df,
    ref_expnr = ref_expnr,
    target_expnr = target_expnr,
    rg = rg,
    lc.err = lc.err,
    err = err,
    db_path = db_path
  )
  
  
  # Save association file
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