#' @title Compare MS/MS spectra of two compounds and compute CSPP similarity metrics
#'
#' @description
#' targMS2comp_SQL() retrieves MS/MS spectra for two compounds from a SQLite
#' database and computes CSPP similarity metrics based on fragment ion and
#' neutral loss comparisons.
#'
#' @details
#' The function extracts one MS/MS spectrum (MS level 2) for each compound from
#' the \code{msms_spectrum} table and retrieves the corresponding fragment peaks
#' from the \code{msms_spectrum_peak} table.
#'
#' Fragment ions are filtered using an intensity threshold and converted to
#' relative intensities. CSPP similarity metrics are then computed using dot
#' product comparisons between shared fragment ions and shared neutral losses.
#'
#' Neutral losses are calculated as the difference between the precursor mass
#' and fragment ion mass. Only peaks exceeding the specified intensity threshold
#' are included in similarity calculations.
#'
#' @param compid1 `Integer(1)` compound identifier corresponding to the
#'   substrate (precursor) compound.
#'
#' @param compid2 `Integer(1)` compound identifier corresponding to the
#'   product compound.
#'
#' @param data_con `DBIConnection` to an open SQLite database containing
#'   MS/MS spectra and fragment peak data.
#'
#' @param IntThres `Numeric(1)` minimum fragment ion intensity threshold.
#'   Only peaks with intensity greater than or equal to this value are considered
#'   during similarity calculations. Default values typically depend on the
#'   instrument type (e.g. \code{100} for FTMS or \code{5} for QTOF data).
#'
#' @return
#' A `data.frame` containing CSPP similarity metrics for the compound pair.
#' The returned columns include precursor and product compound identifiers,
#' precursor masses, number of filtered fragment ions, counts of common ions
#' and neutral losses, dot product similarity scores, and forward/reverse match
#' ratios.
#'
#' If no MS/MS spectra or valid fragment peaks are found for either compound,
#' the function returns \code{NULL}.
#'
#' @author Ahlam Mentag
#'
#' @export
targMS2comp_SQL <- function(compid1,
                            compid2,
                            ms2_split,
                            IntThres = 5) {
  
  ms2_sub  <- ms2_split[[as.character(compid1)]]
  ms2_prod <- ms2_split[[as.character(compid2)]]
  
  if (is.null(ms2_sub) || is.null(ms2_prod))
    return(NULL)
  
  if (nrow(ms2_sub) == 0 || nrow(ms2_prod) == 0)
    return(NULL)
  
  precursor1 <- suppressWarnings(as.numeric(ms2_sub$precursor_mz[1]))
  precursor2 <- suppressWarnings(as.numeric(ms2_prod$precursor_mz[1]))
  
  if (is.na(precursor1) || is.na(precursor2))
    return(NULL)
  
  # Force numeric
  ms2_sub$mz        <- as.numeric(ms2_sub$mz)
  ms2_sub$intensity <- as.numeric(ms2_sub$intensity)
  ms2_prod$mz        <- as.numeric(ms2_prod$mz)
  ms2_prod$intensity <- as.numeric(ms2_prod$intensity)
  
  # Relative intensities
  ms2_sub$rint  <- as.integer(ms2_sub$intensity / max(ms2_sub$intensity, na.rm = TRUE) * 100)
  ms2_prod$rint <- as.integer(ms2_prod$intensity / max(ms2_prod$intensity, na.rm = TRUE) * 100)
  
  
  # IONS 
  
  
  ac <- ms2_sub[, c("mz", "rint", "intensity")]
  ac <- ac[ac$intensity >= IntThres & !is.na(ac$intensity), ]
  
  bd <- ms2_prod[, c("mz", "rint", "intensity")]
  bd <- bd[bd$intensity >= IntThres & !is.na(bd$intensity), ]
  
  if (nrow(ac) == 0 || nrow(bd) == 0)
    return(NULL)
  
  DotIons <- CommonDotProd(ac, bd)
  
  
  # NEUTRAL LOSSES
  
  ac$nloss <- as.integer(round(precursor1) - ac$mz)
  bd$nloss <- as.integer(round(precursor2) - bd$mz)
  
  ac_loss <- ac[, c("nloss", "rint", "intensity")]
  bd_loss <- bd[, c("nloss", "rint", "intensity")]
  
  DotLoss <- CommonDotProd(ac_loss, bd_loss)
  
  # Final output (same structure as original)
  data.frame(
    COMPID.sub   = compid1,
    MZ.sub       = precursor1,
    IONS.sub     = nrow(ac),
    COMPID.prod  = compid2,
    MZ.prod      = precursor2,
    IONS.prod    = nrow(bd),
    COMMON_IONS  = DotIons[[1]],
    DOT_IONS     = round(DotIons[[2]], 2),
    COMMON_LOSS  = DotLoss[[1]],
    DOT_LOSS     = round(DotLoss[[2]], 2),
    FORW_IONS    = round(DotIons[[1]] / nrow(ac), 2),
    REV_IONS     = round(DotIons[[1]] / nrow(bd), 2),
    FORW_LOSS    = round(DotLoss[[1]] / nrow(ac), 2),
    REV_LOSS     = round(DotLoss[[1]] / nrow(bd), 2),
    stringsAsFactors = FALSE
  )
}