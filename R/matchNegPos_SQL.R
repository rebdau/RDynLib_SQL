#' @title Match  MS2 spectra between the same instruments type 
#'        using SQL databases. 
#'
#' @description
#' 
#' `matchNegPos_SQL()` Removes false positive negFTm/z,posFTm/z peak pairs by 
#                      checking whether a certain percentage of the neg MS2  
#                      spectral peaks can be traced in the pos MS2 spectrum. 
#                      By default, 20% of the neg MS2 peaks are traced in the
#                      pos MS2 spectrum.
#'
#' @param LCal 
#' `data.frame` returned by `Aligning_General_SQL()` function.
#' 
#' @param FTn_con `DBIConnection`  
#'   A `DBI connection` object to the reference SQLite database.
#'
#' @param FTp_con `DBIConnection`  
#'   A `DBI connection` object to the SQLite database to align with.
#'
#' @param polarity_ftn 
#'        `Integer (0/1)` Polarity of negative FTMS or QTOF experiment.
#'        
#' @param polarity_ftp 
#'        `Integer (0/1)` Polarity of positive FTMS or QTOF experiment.
#' 
#' @param minPeaks `Numeric(1)` Minimum matching ratio (default 0.2).
#' 
#'
#' @return Filtered `LCal` with MS2-supported neg–pos pairs.
#'
#' @import DBI RSQLite
#' @export
matchNegPos_SQL <- function(
    LCal,
    FTn_con,
    FTp_con,
    polarity_ftn = 0,
    polarity_ftp = 1,
    minIon = 0.2,
    tol = 0.01,
    spectrum_type_FTn = NULL,
    spectrum_type_FTp = NULL
) {
  
  if (is.null(LCal) || nrow(LCal) == 0) {
    return(LCal)
  }
  
  
  # Check whether spectrum_type exists
  has_type_FTn <- "spectrum_type" %in%
    DBI::dbListFields(FTn_con, "msms_spectrum")
  
  has_type_FTp <- "spectrum_type" %in%
    DBI::dbListFields(FTp_con, "msms_spectrum")
  
  
  # If spectrum_type exists, user must provide a value
  if (has_type_FTn && is.null(spectrum_type_FTn)) {
    stop(
      "`spectrum_type_FTn` must be provided because ",
      "`spectrum_type` exists in the FTn database."
    )
  }
  
  if (has_type_FTp && is.null(spectrum_type_FTp)) {
    stop(
      "`spectrum_type_FTp` must be provided because ",
      "`spectrum_type` exists in the FTp database."
    )
  }
  
  
  # Build FTn query
  query_FTn <- "
    SELECT
      s.compound_id,
      p.mz
    FROM msms_spectrum s
    JOIN msms_spectrum_peak p
      USING(spectrum_id)
    WHERE s.ms_level = ?
      AND s.polarity = ?
  "
  
  params_FTn <- list(
    2,
    polarity_ftn
  )
  
  
  # Add spectrum type filter only if column exists
  if (has_type_FTn) {
    
    query_FTn <- paste0(
      query_FTn,
      " AND s.spectrum_type = ?"
    )
    
    params_FTn <- c(
      params_FTn,
      list(spectrum_type_FTn)
    )
  }
  
  
  # Build FTp query
  query_FTp <- "
    SELECT
      s.compound_id,
      p.mz
    FROM msms_spectrum s
    JOIN msms_spectrum_peak p
      USING(spectrum_id)
    WHERE s.ms_level = ?
      AND s.polarity = ?
  "
  
  params_FTp <- list(
    2,
    polarity_ftp
  )
  
  
  # Add spectrum type filter only if column exists
  if (has_type_FTp) {
    
    query_FTp <- paste0(
      query_FTp,
      " AND s.spectrum_type = ?"
    )
    
    params_FTp <- c(
      params_FTp,
      list(spectrum_type_FTp)
    )
  }
  
  
  # Extract MS2 peaks
  ms2_neg_df <- DBI::dbGetQuery(
    FTn_con,
    query_FTn,
    params = params_FTn
  )
  
  ms2_pos_df <- DBI::dbGetQuery(
    FTp_con,
    query_FTp,
    params = params_FTp
  )
  
  
  # Check spectra
  if (nrow(ms2_neg_df) == 0) {
    warning("No matching MS2 spectra found in the FTn database.")
    return(LCal[0, , drop = FALSE])
  }
  
  if (nrow(ms2_pos_df) == 0) {
    warning("No matching MS2 spectra found in the FTp database.")
    return(LCal[0, , drop = FALSE])
  }
  
  
  # Split peaks by compound
  ms2_neg <- split(
    ms2_neg_df$mz,
    ms2_neg_df$compound_id
  )
  
  ms2_pos <- split(
    ms2_pos_df$mz,
    ms2_pos_df$compound_id
  )
  
  
  # Filter candidate pairs
  i <- 1
  
  while (i <= nrow(LCal)) {
    
    neg_id <- LCal[i, 1]
    pos_id <- LCal[i, 7]
    
    neg_peaks <- ms2_neg[[as.character(neg_id)]]
    pos_peaks <- ms2_pos[[as.character(pos_id)]]
    
    
    # Remove pairs without MS2 spectra
    if (is.null(neg_peaks) || is.null(pos_peaks)) {
      LCal <- LCal[-i, , drop = FALSE]
      next
    }
    
    
    # Remove missing values
    neg_peaks <- neg_peaks[!is.na(neg_peaks)]
    pos_peaks <- pos_peaks[!is.na(pos_peaks)]
    
    if (length(neg_peaks) == 0 ||
        length(pos_peaks) == 0) {
      
      LCal <- LCal[-i, , drop = FALSE]
      next
    }
    
    
    # Adjust positive-mode fragments
    pos_peaks <- pos_peaks - 2
    
    
    # Count matching ions
    same_ion <- sum(
      vapply(
        neg_peaks,
        function(x) {
          any(abs(pos_peaks - x) <= tol)
        },
        logical(1)
      )
    )
    
    
    # Fraction of FTn ions found in FTp
    ratio <- same_ion / length(neg_peaks)
    
    
    if (ratio < minIon) {
      LCal <- LCal[-i, , drop = FALSE]
      next
    }
    
    
    i <- i + 1
  }
  
  
  unique(LCal)
}