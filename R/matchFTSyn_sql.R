#' @title Match FT and QTOF MS2 spectra based on local alignments 
#'
#' @description
#' `matchFTSyn_SQL()` filters locally aligned FTMS–QTOF candidate pairs in
#' (`LCal`) by comparing MS2 peak lists stored in SQL databases.
#'
#' @details
#' The function:
#' - Loads FT and QTOF MS2 spectra from SQLite via DBI
#' - Applies polarity and spectrum-type filtering
#' - Computes MS2 peak overlap for each FT–QTOF candidate pair
#' - Removes pairs with less matched peaks based on a threshold
#'
#' @param LCal 
#' `data.frame` A data.frame returned by `Aligning_General_SQL()` function.
#' 
#' @param FT_con `character(1)` A DBI connection object 
#' to the FTMS SQLite database (the reference database).
#' 
#' @param QTOF_con `character(1)` A DBI connection object to the QTOF
#' SQLite database.
#' 
#' @param polarity_ft `numeric(1)` Integer (0 or 1), refers to the polarity 
#' of the ftms  experiment to align.
#' 
#' @param polarity_qtof `numeric(1)`Integer (0 or 1), refers to the polarity 
#' of the qtof experiment to align.
#' 
#' @param minPeaks `Numeric(1)` Minimum MS2 matching ratio required.
#' 
#' @param spectrum_type_FT `character(1)` Spectrum type to use for the FT-MS
#'        spectra. If `NULL` and the database does not contain a `spectrum_type`
#'        column, all spectra are used.
#'
#' @param spectrum_type_QTOF `character(1)` Spectrum type to use for the QTOF
#'        spectra. If `NULL` and the database does not contain a `spectrum_type`
#'        column, all spectra are used.
#'
#' @return A filtered version of `LCal`, keeping only MS2 supported matches.
#'
#' @importFrom DBI dbConnect
#'
#' @importFrom DBI dbDisconnect
#'
#' @importFrom RSQLite SQLite
#' 
#' @importFrom DBI dbGetQuery
#' 
#' @importMethodsFrom DBI dbListFields
#' 
#' @author Ahlam Mentag
#'
#' @export
matchFTSyn_SQL <- function(
    LCal,
    FT_con,
    QTOF_con,
    FT_expnr,
    QTOF_expnr,
    polarity_ft = 0,
    polarity_qtof = 0,
    minPeaks = 0.6,
    spectrum_type_FT = NULL,
    spectrum_type_QTOF = NULL
) {
  
  if (is.null(LCal) || nrow(LCal) == 0) {
    message(
      "No candidate alignments found; returning empty LCal."
    )
    return(LCal)
  }
  
  
  # Check spectrum_type columns
  ft_cols <- DBI::dbListFields(
    FT_con,
    "msms_spectrum"
  )
  
  qtof_cols <- DBI::dbListFields(
    QTOF_con,
    "msms_spectrum"
  )
  
  has_type_FT <- "spectrum_type" %in% ft_cols
  has_type_QTOF <- "spectrum_type" %in% qtof_cols
  
  
  # Check spectrum type arguments
  if (has_type_FT && is.null(spectrum_type_FT)) {
    stop(
      "`spectrum_type_FT` must be provided because ",
      "`spectrum_type` exists in the FT database."
    )
  }
  
  if (has_type_QTOF && is.null(spectrum_type_QTOF)) {
    stop(
      "`spectrum_type_QTOF` must be provided because ",
      "`spectrum_type` exists in the QTOF database."
    )
  }
  
  
  # Load compounds from selected FT experiment
  subdb_ft <- DBI::dbGetQuery(
    FT_con,
    "
    SELECT
      retention_time,
      mass_measured,
      compound_id
    FROM ms_compound
    WHERE expid = ?
    ORDER BY compound_id
    ",
    params = list(
      FT_expnr
    )
  )
  
  
  # Load compounds from selected QTOF experiment
  subdb_syn <- DBI::dbGetQuery(
    QTOF_con,
    "
    SELECT
      retention_time,
      mass_measured,
      compound_id
    FROM ms_compound
    WHERE expid = ?
    ORDER BY compound_id
    ",
    params = list(
      QTOF_expnr
    )
  )
  
  
  # Build FT query
  ft_query <- "
    SELECT
      s.compound_id,
      p.mz
    FROM msms_spectrum s
    JOIN msms_spectrum_peak p
      USING(spectrum_id)
    JOIN ms_compound c
      USING(compound_id)
    WHERE s.ms_level = ?
      AND s.polarity = ?
      AND c.expid = ?
  "
  
  ft_params <- list(
    2,
    polarity_ft,
    FT_expnr
  )
  
  
  # Add spectrum_type filter only if column exists
  if (has_type_FT) {
    
    ft_query <- paste0(
      ft_query,
      " AND s.spectrum_type = ?"
    )
    
    ft_params <- c(
      ft_params,
      list(spectrum_type_FT)
    )
  }
  
  
  ft_query <- paste0(
    ft_query,
    " ORDER BY s.compound_id, p.mz"
  )
  
  
  # Retrieve FT MS2 peaks
  ms2_ft_df <- DBI::dbGetQuery(
    FT_con,
    ft_query,
    params = ft_params
  )
  
  
  # Build QTOF query
  qtof_query <- "
    SELECT
      s.compound_id,
      p.mz
    FROM msms_spectrum s
    JOIN msms_spectrum_peak p
      USING(spectrum_id)
    JOIN ms_compound c
      USING(compound_id)
    WHERE s.ms_level = ?
      AND s.polarity = ?
      AND c.expid = ?
  "
  
  qtof_params <- list(
    2,
    polarity_qtof,
    QTOF_expnr
  )
  
  
  # Add spectrum_type filter only if column exists
  if (has_type_QTOF) {
    
    qtof_query <- paste0(
      qtof_query,
      " AND s.spectrum_type = ?"
    )
    
    qtof_params <- c(
      qtof_params,
      list(spectrum_type_QTOF)
    )
  }
  
  
  qtof_query <- paste0(
    qtof_query,
    " ORDER BY s.compound_id, p.mz"
  )
  
  
  # Retrieve QTOF MS2 peaks
  ms2_syn_df <- DBI::dbGetQuery(
    QTOF_con,
    qtof_query,
    params = qtof_params
  )
  
  
  # Check spectra
  if (nrow(ms2_ft_df) == 0) {
    
    warning(
      "No FT MS2 spectra found for experiment ",
      FT_expnr,
      "."
    )
    
    return(
      LCal[0, , drop = FALSE]
    )
  }
  
  
  if (nrow(ms2_syn_df) == 0) {
    
    warning(
      "No QTOF MS2 spectra found for experiment ",
      QTOF_expnr,
      "."
    )
    
    return(
      LCal[0, , drop = FALSE]
    )
  }
  
  
  # Round and split FT spectra by compound
  ms2_ft_list <- split(
    round(ms2_ft_df$mz),
    ms2_ft_df$compound_id
  )
  
  ms2_ft_list <- lapply(
    ms2_ft_list,
    unique
  )
  
  
  # Round and split QTOF spectra by compound
  ms2_syn_list <- split(
    round(ms2_syn_df$mz),
    ms2_syn_df$compound_id
  )
  
  ms2_syn_list <- lapply(
    ms2_syn_list,
    unique
  )
  
  
  # Align MS2 lists by compound ID
  ms2.ft <- ms2_ft_list[
    as.character(
      subdb_ft$compound_id
    )
  ]
  
  ms2.syn <- ms2_syn_list[
    as.character(
      subdb_syn$compound_id
    )
  ]
  
  
  # Replace missing spectra with empty vectors
  ms2.ft[
    sapply(ms2.ft, is.null)
  ] <- list(integer(0))
  
  ms2.syn[
    sapply(ms2.syn, is.null)
  ] <- list(integer(0))
  
  
  # Filter candidate pairs
  i <- 1
  
  while (i <= nrow(LCal)) {
    
    ft.compid <- LCal[i, 1]
    syn.compid <- LCal[i, 7]
    
    
    ft.row <- which(
      subdb_ft$compound_id == ft.compid
    )
    
    syn.row <- which(
      subdb_syn$compound_id == syn.compid
    )
    
    
    # Compound not found
    if (length(ft.row) == 0 ||
        length(syn.row) == 0) {
      
      LCal <- LCal[
        -i,
        ,
        drop = FALSE
      ]
      
      next
    }
    
    
    # Retrieve spectra
    ms2ion <- ms2.ft[[ft.row[1]]]
    msmsion <- ms2.syn[[syn.row[1]]]
    
    
    # Remove pairs without spectra
    if (length(ms2ion) == 0 ||
        length(msmsion) == 0) {
      
      LCal <- LCal[
        -i,
        ,
        drop = FALSE
      ]
      
      next
    }
    
    
    # Number of matching fragment ions
    same_ion <- sum(
      ms2ion %in% msmsion
    )
    
    
    # Fraction of FT ions found in QTOF
    matching_ratio <- same_ion /
      length(ms2ion)
    
    
    if (matching_ratio < minPeaks) {
      
      LCal <- LCal[
        -i,
        ,
        drop = FALSE
      ]
      
      next
    }
    
    
    i <- i + 1
  }
  
  
  unique(LCal)
}