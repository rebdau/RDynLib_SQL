#' Combine spectra across a tree
#'
#' Given a list of spectra (matrices with columns `mz` and `intensity`), 
#' this function sums the intensities of matching m/z values across all spectra
#' and returns a single merged spectrum matrix.
#'
#' @param spectra_list A list of numeric matrices, each with columns `mz` and `intensity`.
#'
#' @return A numeric matrix with columns `mz` and `intensity`, 
#'         containing combined intensities for matching m/z values, sorted by m/z.
#' @importFrom dplyr group_by summarise arrange
#' @export
combine_spectra <- function(spectra_list) {
  do.call(rbind, spectra_list) %>%
    as.data.frame() %>%
    group_by(mz) %>%
    summarise(intensity = sum(intensity), .groups = "drop") %>%
    arrange(mz) %>%
    as.matrix()
}