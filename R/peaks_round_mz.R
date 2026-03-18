#' Round m/z values in a spectrum and keep only highest intensity peaks
#'
#' Given a matrix of spectral data (m/z and intensity), this function rounds
#' the m/z values to integers, and for any duplicated m/z values, retains only
#' the peak with the highest intensity.
#'
#' @param x A numeric matrix with at least two columns: m/z (column 1) and intensity (column 2).
#' @param spectrumMsLevel Integer. The MS level of the spectrum (currently unused, placeholder for compatibility).
#' @param ... Additional arguments (not currently used).
#' @return A numeric matrix of two columns: rounded m/z and corresponding intensity.
#'   Only the highest intensity peak is kept for duplicated m/z values.
#' @examples
#' mat <- matrix(c(100.2, 200, 100.7, 150, 101.4, 50), ncol = 2, byrow = TRUE)
#' peaks_round_mz(mat, spectrumMsLevel = 1)
#' @importFrom dplyr group_by slice_max arrange
#' @export
peaks_round_mz <- function(x, spectrumMsLevel, ...) {
  x[, 1L] <- round_perl(x[, 1L])
  df <- as.data.frame(x)
  colnames(df) <- c("mz", "intensity")
  
  # Keep only highest intensity for duplicated mz
  df <- df %>%
    group_by(mz) %>%
    slice_max(order_by = intensity, n = 1, with_ties = FALSE) %>%
    arrange(mz)
  
  # Return as matrix
  as.matrix(df)
}