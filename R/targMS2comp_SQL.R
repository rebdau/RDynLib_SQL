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
                            IntThres = 0,
                            mz_tol = 0.3) {
  
  compid1 <- as.character(compid1)
  compid2 <- as.character(compid2)
  
  ms2_sub  <- ms2_split[[compid1]]
  ms2_prod <- ms2_split[[compid2]]
  
  if (is.null(ms2_sub) || is.null(ms2_prod)) return(NULL)
  if (nrow(ms2_sub) == 0 || nrow(ms2_prod) == 0) return(NULL)
  
  required <- c("mz", "intensity", "precursor_mz")
  if (!all(required %in% names(ms2_sub)) ||
      !all(required %in% names(ms2_prod))) return(NULL)
  
  precursor1 <- as.numeric(ms2_sub$precursor_mz[1])
  precursor2 <- as.numeric(ms2_prod$precursor_mz[1])
  
  if (is.na(precursor1) || is.na(precursor2)) return(NULL)
  
  # Keep only relevant columns
  ac <- ms2_sub[, c("mz", "intensity")]
  bd <- ms2_prod[, c("mz", "intensity")]
  
  # FIX: explicit column reference (avoids "closure" bug)
  ac <- ac[!is.na(ac$intensity) & ac$intensity >= IntThres, , drop = FALSE]
  bd <- bd[!is.na(bd$intensity) & bd$intensity >= IntThres, , drop = FALSE]
  
  if (nrow(ac) == 0 || nrow(bd) == 0) return(NULL)
  
  # Normalize intensities
  ac$rint <- ac$intensity / max(ac$intensity) * 100
  bd$rint <- bd$intensity / max(bd$intensity) * 100
  
  # Dot product on ions
  DotIons <- CommonDotProd(
    ac[, c("mz", "rint")],
    bd[, c("mz", "rint")],
    tol = mz_tol
  )
  
  # Loss calculation
  ac$nloss <- precursor1 - ac$mz
  bd$nloss <- precursor2 - bd$mz
  
  DotLoss <- CommonDotProd(
    ac[, c("nloss", "rint")],
    bd[, c("nloss", "rint")],
    tol = mz_tol
  )
  
  # Output
  data.frame(
    COMPID.sub  = compid1,
    COMPID.prod = compid2,
    DOT_IONS    = DotIons[[2]],
    DOT_LOSS    = DotLoss[[2]],
    COMMON_IONS = DotIons[[1]],
    COMMON_LOSS = DotLoss[[1]],
    FORW_IONS   = DotIons[[1]] / nrow(ac),
    REV_IONS    = DotIons[[1]] / nrow(bd),
    FORW_LOSS   = DotLoss[[1]] / nrow(ac),
    REV_LOSS    = DotLoss[[1]] / nrow(bd)
  )
}