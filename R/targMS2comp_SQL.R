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
                            ms2_split) {
  
  ms2_sub  <- ms2_split[[as.character(compid1)]]
  ms2_prod <- ms2_split[[as.character(compid2)]]
  
  if (is.null(ms2_sub) || is.null(ms2_prod))
    return(NULL)
  
  if (nrow(ms2_sub) == 0 || nrow(ms2_prod) == 0)
    return(NULL)
  
  precursor1 <- ms2_sub$precursor_mz[1]
  precursor2 <- ms2_prod$precursor_mz[1]
  
  # Relative intensities
  ms2_sub$rint  <- as.integer(ms2_sub$intensity /
                                max(ms2_sub$intensity) * 100)
  
  ms2_prod$rint <- as.integer(ms2_prod$intensity /
                                max(ms2_prod$intensity) * 100)
  
  ac <- ms2_sub[, c("mz", "rint", "intensity")]
  bd <- ms2_prod[, c("mz", "rint", "intensity")]
  
  if (nrow(ac) == 0 || nrow(bd) == 0)
    return(NULL)
  
  DotIons <- CommonDotProd(ac, bd)
  
  # Neutral losses
  ac$nloss <- round(precursor1) - ac$mz
  bd$nloss <- round(precursor2) - bd$mz
  
  DotLoss <- CommonDotProd(
    ac[, c("nloss", "rint", "intensity")],
    bd[, c("nloss", "rint", "intensity")]
  )
  
  data.frame(
    COMPID.sub   = compid1,
    MZ.sub       = precursor1,
    IONS.sub     = nrow(ac),
    COMPID.prod  = compid2,
    MZ.prod      = precursor2,
    IONS.prod    = nrow(bd),
    COMMON_IONS  = DotIons[[1]],
    DOT_IONS     = round(DotIons[[2]], 3),
    COMMON_LOSS  = DotLoss[[1]],
    DOT_LOSS     = round(DotLoss[[2]], 3),
    FORW_IONS    = round(DotIons[[1]] / nrow(ac), 3),
    REV_IONS     = round(DotIons[[1]] / nrow(bd), 3),
    FORW_LOSS    = round(DotLoss[[1]] / nrow(ac), 3),
    REV_LOSS     = round(DotLoss[[1]] / nrow(bd), 3),
    stringsAsFactors = FALSE
  )
}