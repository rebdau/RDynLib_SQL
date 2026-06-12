
#' Extract full MSn fragmentation trees for features
#'
#' This function extracts complete MSn fragmentation trees (MS2, MS3 and
#' MS4 spectra) associated with features from an `XcmsExperiment` object.
#' Starting from feature-associated MS2 spectra, higher-level MSn spectra
#' are retrieved recursively using the `precScanNum` spectra variable.
#'
#' All spectra belonging to the same fragmentation tree are linked to the
#' originating feature through a propagated `feature_id`.
#'
#' @param ftms An \linkS4class{XcmsExperiment} object containing detected
#'   features and multi-level MSn data of a metabolomics experiment. 
#'   The data is expected to include at least MS2 and MS3 spectra and optionally 
#'   higher-order fragmentation levels (e.g. MS4 and above).
#'
#' @param ms2 A \linkS4class{Spectra} object containing all MS2 spectra whose
#'   retention time and precursor m/z fall within the ranges of any
#'   chromatographic peak of a feature in `ftms`, as returned for example by
#'   \code{xcms::featureSpectra(ftms, msLevel = 2L)}.
#'   
#' @return A \linkS4class{Spectra} object containing all MS2, MS3 and MS4
#'   spectra associated with the input features. Each spectrum includes a
#'   `feature_id` spectra variable identifying the originating feature.
#'
#' @details
#' The function reconstructs MSn fragmentation trees by linking spectra
#' across MS levels within each data origin. For each distinct `dataOrigin`
#' in the input `ms2` object, MS2 spectra are matched to MS3 spectra using
#' the `precScanNum` and `acquisitionNum` spectra variables. The same
#' matching strategy is then applied to link MS3 spectra to MS4 spectra.
#'
#' The `feature_id` from the MS2 spectra is propagated to all matched MS3 and
#' MS4 spectra, ensuring that all spectra belonging to the same fragmentation
#' tree can be traced back to the originating feature.
#'
#' Spectra are processed independently for each `dataOrigin` and subsequently
#' combined into a single \linkS4class{Spectra} object.
#'
#' @section Processing steps:
#' \enumerate{
#'   \item Split MS2 spectra by `dataOrigin`.
#'   \item For each subset, retrieve MS3 spectra from `ftms` with matching
#'         `dataOrigin`.
#'   \item Match MS2 to MS3 spectra using `acquisitionNum` and `precScanNum`
#'         via \code{findMatches()}.
#'   \item Propagate `feature_id` from MS2 to matched MS3 spectra.
#'   \item Retrieve MS4 spectra and match them to MS3 spectra using the same
#'         approach.
#'   \item Propagate `feature_id` from MS3 to MS4 spectra.
#'   \item Combine MS2, MS3, and MS4 spectra for each `dataOrigin`.
#'   \item Concatenate results across all data origins into a single object.
#' }
#'
#' @examples
#' ## Assume 'ftms' is an XcmsExperiment object
#' ## Extract MS2 spectra associated with features
#' ms2 <- xcms::featureSpectra(ftms, msLevel = 2L)
#'
#' ## Retrieve full MSn trees
#' msn <- ftms_all_levels(ftms, ms2)
#'
#' length(msn)
#'
#' @importFrom xcms featureSpectra
#' @import Spectra
#' @importFrom S4Vectors findMatches
#' @export
ftms_all_levels <- function(ftms, ms2) {
  
  cat("Number of MS2 spectra found:", length(ms2), "\n")
  
  res <- lapply(unique(ms2$dataOrigin), function(origin) {
    ms2_subset <- filterDataOrigin(ms2, origin)
    
    ms3_filtered <- filterDataOrigin(filterMsLevel(spectra(ftms), 3), origin)
    ## To support n:m matches
    m <- findMatches(ms2_subset$acquisitionNum, ms3_filtered$precScanNum)
    cat("Number of MS3 matched to MS2:", length(m), "\n")
    ms3_filtered <- ms3_filtered[to(m)]
    ms3_filtered$feature_id <- ms2_subset$feature_id[from(m)]
    
    ms4_filtered <- filterDataOrigin(filterMsLevel(spectra(ftms), 4), origin)
    m <- findMatches(ms3_filtered$acquisitionNum, ms4_filtered$precScanNum)
    cat("Number of MS4 matched to MS3:", length(m), "\n")
    ms4_filtered <- ms4_filtered[to(m)]
    ms4_filtered$feature_id <- ms3_filtered$feature_id[from(m)]
    
    # Combine spectra for this origin
    c(ms2_subset, ms3_filtered, ms4_filtered)
  })
  do.call(c, res)
}
