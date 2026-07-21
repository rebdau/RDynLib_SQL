#' @tiltle filter out all peaks that have intensity below 1% (or 0.5%) of 
#' the base peak.
#' 
#' @param sps spectra object to filter
#' @param threshold 'numeric()' filtering threshold
#' 
#' @author Ahlam Mentag
#' 
#' @export
filter_basepeak <- function(sps, threshold = 0.01) {
  
  addProcessing(
    sps,
    function(x) {
      
      if (nrow(x) == 0)
        return(x)
      
      base_peak <- max(x[, 2], na.rm = TRUE)
      
      x[x[, 2] >= threshold * base_peak, , drop = FALSE]
    }
  )
  
}