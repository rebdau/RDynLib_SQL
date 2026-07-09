## Helper function: calculate rounded neutral losses
neutral_loss <- function(mz_values, precursorMz) {
  

  
  round_perl(precursorMz- mz_values)
}

## Helper function: round the mz values depending on the resolution type chosen by the user.
round_perl <- function(number) {
  

    floor(number + 0.5)
}

## Helper function: remove the duplicate mz values after the rounding.
remove_duplicates <- function(mz, intensity) {
  

  mz <- round_perl(mz)
  
  unique_mz <- unique(mz)
  
  list(
    mz = unique_mz,
    intensity = vapply(
      unique_mz,
      function(m)
        max(intensity[mz == m], na.rm = TRUE),
      numeric(1)
    )
  )
}

#' Calculate spectral similarity using the DynLib symmetric dot product
#'
#' Function to calculate the symmetric dot product between two spectra. As
#' input **matched** peak matrices are expected, such as those returned by
#' `dynlib_map()`. These input matrices are expected to have the same number
#' of rows, with the same row in both matrices `x` and `y` representing a
#' matching peak. The matrices can also contain neutral loss peaks. Importantly
#' each matrix has to also contain the m/z-weighted intensity sum as an
#' *attribute* to the matrix (i.e. `attributes(x)$wintensity_sum)`).
#'
#' @param x `numeric` 2-column `matrix` with the m/z and intensity values of
#'     **matched** peaks.
#'
#' @param y `numeric` 2-column `matrix` with the m/z and intensity values of
#'     **matched** peaks.
#'
#' @param n  `numeric`; exponent applied to intensities in the weighted dot-product calculation.
#'
#' @param m  `numeric`; exponent applied to m/z values in the weighted dot-product calculation.
#'
#' @return a `numeric(1)` with the similarity between the matched peaks.
#'
#' @author Ahlam Mentag
dynlib_symmetric_dotproduct <- function(x, y, n = 3, m = 0.6, ...) {
  
  if (!nrow(x) || !nrow(y))
    return(0.0)
  
  mz1 <- x[, 1]
  intensity1 <- x[, 2]
  mz2 <- y[, 1]
  intensity2 <- y[, 2]
  
  intensity1 <- (intensity1^n) * (mz1^m)
  intensity2 <- (intensity2^n) * (mz2^m)
  
  sum(intensity1 * intensity2)^2 /
    (attr(x, "wintensity_sum") * attr(y, "wintensity_sum"))
}


#' Match fragment and neutral-loss peaks between spectra
#'
#' Matches/maps peaks between the two provided peak matrices.
#'
#' @return
#'
#' a list with numeric matrices containing only the matched peaks between
#' the two input spectra. The m/z weighted intensity sum of each (cleaned) peak
#' matrix is returned as an attribute `"wintensity_sum"` to each table.
#'
#' @author Ahlam Mentag
dynlibmatch2_map <- function(
    x, y,
    xPrecursorMz, yPrecursorMz,
    n = 3, m = 0.6,
    fragments_method = c("unit.resolution", "high.resolution"),
    tolerance = 0, ppm = 0,
    ...
) {
  
  fragments_method <- match.arg(fragments_method)

  ## Return empty aligned matrices if one spectrum has no peaks
  if (!nrow(x) || !nrow(y)) {
    empty <- matrix(numeric(), ncol = 2, nrow = 0)
    return(list(empty, empty))
  }
  

  ## Clean fragment peaks

  cleaned1 <- remove_duplicates(
    x[, 1], x[, 2]
  )
  
  cleaned2 <- remove_duplicates(
    y[, 1], y[, 2]
  )
  
  mz1 <- cleaned1$mz
  int1 <- cleaned1$intensity
  
  mz2 <- cleaned2$mz
  int2 <- cleaned2$intensity
  

  ## Exact fragment matching

  idx2 <- match(mz1, mz2)
  valid_frag <- !is.na(idx2)
  
  matched1_frag <- cbind(
    mz = mz1[valid_frag],
    intensity = int1[valid_frag]
  )
  
  matched2_frag <- cbind(
    mz = mz2[idx2[valid_frag]],
    intensity = int2[idx2[valid_frag]]
  )
  

  ## Neutral loss matching

  remaining1 <- !valid_frag
  used2 <- rep(FALSE, length(mz2))
  used2[idx2[valid_frag]] <- TRUE
  remaining2 <- !used2
  
  matched1_nl <- NULL
  matched2_nl <- NULL
  
  if (any(remaining1) && any(remaining2)) {
    
    nl1 <- neutral_loss(
      mz1[remaining1],
      xPrecursorMz
    )
    
    nl2 <- neutral_loss(
      mz2[remaining2],
      yPrecursorMz
    )
    
    idx_nl2 <- match(nl1, nl2)
    valid_nl <- which(!is.na(idx_nl2))
    
    if (length(valid_nl) > 0) {
      
      matched1_nl <- cbind(
        mz = mz1[remaining1][valid_nl],
        intensity = int1[remaining1][valid_nl]
      )
      
      matched2_nl <- cbind(
        mz = mz2[remaining2][idx_nl2[valid_nl]],
        intensity = int2[remaining2][idx_nl2[valid_nl]]
      )
    }
  }
  

  ## Combine matched peaks

  matched1 <- matched1_frag
  matched2 <- matched2_frag
  
  if (!is.null(matched1_nl)) {
    matched1 <- rbind(matched1, matched1_nl)
    matched2 <- rbind(matched2, matched2_nl)
  }
  
  
  if (nrow(matched1) != nrow(matched2)) {
    stop("dynlibmatch_map(): internal error - unmatched row counts.")
  }
  
  attr(matched1, "wintensity_sum") <-
    sum(((int1^n) * (mz1^m))^2)
  
  attr(matched2, "wintensity_sum") <-
    sum(((int2^n) * (mz2^m))^2)
  
  return(list(matched1, matched2))
}


