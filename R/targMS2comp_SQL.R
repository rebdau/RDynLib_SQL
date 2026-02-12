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
targMS2comp_SQL <- function(compid1, compid2,
                            data_con, IntThres) {

    
    
    ## get one MS2 spectrum per compound 
    spec_q <- "
    SELECT spectrum_id, compound_id, precursor_mz
    FROM msms_spectrum
    WHERE compound_id = ?
      AND ms_level = 2
    ORDER BY spectrum_id
    LIMIT 1
  "
    
    s1 <- dbGetQuery(data_con, spec_q, params = list(compid1))
    s2 <- dbGetQuery(data_con, spec_q, params = list(compid2))
    
    if (nrow(s1) == 0 || nrow(s2) == 0) {
      return(NULL)
    }
    
    ##get fragment peaks
    peak_q <- "
    SELECT mz, intensity
    FROM msms_spectrum_peak
    WHERE spectrum_id = ?
  "
    
    p1 <- dbGetQuery(data_con, peak_q, params = list(s1$spectrum_id))
    p2 <- dbGetQuery(data_con, peak_q, params = list(s2$spectrum_id))
    
    if (nrow(p1) == 0 || nrow(p2) == 0) {
      return(NULL)
    }
    
    ## relative intensities
    p1$rint <- as.integer(p1$intensity / max(p1$intensity) * 100)
    p2$rint <- as.integer(p2$intensity / max(p2$intensity) * 100)
    
    ##  filter by intensity 
    ac <- p1[p1$intensity >= IntThres, c("mz", "rint", "intensity")]
    bd <- p2[p2$intensity >= IntThres, c("mz", "rint", "intensity")]
    
    if (nrow(ac) == 0 || nrow(bd) == 0) {
      return(NULL)
    }
    
    ##dot product on ions
    DotIons <- CommonDotProd(ac, bd)
    
    ## neutral losses 
    ac$nloss <- round(s1$precursor_mz) - ac$mz
    bd$nloss <- round(s2$precursor_mz) - bd$mz
    
    ac_n <- ac[, c("nloss", "rint", "intensity")]
    bd_n <- bd[, c("nloss", "rint", "intensity")]
    
    DotLoss <- CommonDotProd(ac_n, bd_n)
    
    ## final output
    fin <- data.frame(
      COMPID.sub   = compid1,
      MZ.sub       = s1$precursor_mz,
      IONS.sub     = nrow(ac),
      COMPID.prod  = compid2,
      MZ.prod      = s2$precursor_mz,
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
    
    return(fin)
  }
  