library(Spectra)



neutral_loss <- function(mz_values, precursorMz) {
  round_perl(precursorMz) - mz_values
}
round_perl <- function(number) {
  return(round(number, 2))
}



remove_duplicates <- function(mz, intensity) {
  mz <- round_perl(mz)
  unique_mz <- unique(mz)
  list(mz = unique_mz,
       intensity = vapply(unique_mz, function(m) max(intensity[mz == m],
                                                     na.rm = TRUE), NA_real_))
}

#' @description
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
#' @param n
#'
#' @param m
#'
#' @return a `numeric(1)` with the similarity between the matched peaks.
#'
#' @author Ahlam Mentag
dynlib_symmetric_dotproduct <- function(x, y, n = 3, m = 0.6, ...) {
  if (!nrow(x) || !nrow(y))
    return(0.0)
  
  mz1 <- x[, 1L]
  intensity1 <- x[, 2L]
  mz2 <- y[, 1L]
  intensity2 <- y[, 2L]
  
  intensity1 <- (intensity1^n) * (mz1^m)
  intensity2 <- (intensity2^n) * (mz2^m)
  
  sum(intensity1 * intensity2)^2 /
    (attributes(x)$wintensity_sum * attributes(y)$wintensity_sum)
}

#' @description
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
dynlib_map <- function(x, y, xPrecursorMz, yPrecursorMz, n = 3, m = 0.6, ...) {
  if (!nrow(x) || !nrow(y))
    return(list(x = matrix(numeric(), ncol = 2, nrow = 0),
                y = matrix(numeric(), ncol = 2, nrow = 0)))

  
  cleaned1 <- remove_duplicates(x[, 1L], x[, 2L])
  cleaned2 <- remove_duplicates(y[, 1L], y[, 2L])
  
  mz1 <- cleaned1$mz
  intensity1 <- cleaned1$intensity
  mz2 <- cleaned2$mz
  intensity2 <- cleaned2$intensity
  
  # Matching exact m/z
  common_mz <- intersect(mz1, mz2)
  keep_mz1 <- mz1 %in% common_mz
  keep_mz2 <- mz2 %in% common_mz
  matched1 <- cbind(mz = mz1[keep_mz1], intensity = intensity1[keep_mz1])
  matched2 <- cbind(mz = mz2[keep_mz2], intensity = intensity2[keep_mz2])
  
  ## Process the remaining peaks to see if they could match as neutral losses
  remaining_mz1 <- !keep_mz1
  remaining_mz2 <- !keep_mz2
  if (any(remaining_mz1) && any(remaining_mz2)) {
    nl1 <- neutral_loss(mz1[remaining_mz1], xPrecursorMz)
    nl2 <- neutral_loss(mz2[remaining_mz2], yPrecursorMz)
    common_nl <- intersect(nl1, nl2)
    if (length(common_nl)) {
      keep_mz1 <- nl1 %in% common_nl
      keep_mz2 <- nl2 %in% common_nl
      matched_nl1 <- cbind(mz = mz1[remaining_mz1][keep_mz1],
                           intensity = intensity1[remaining_mz1][keep_mz1])
      matched_nl2 <- cbind(mz = mz2[remaining_mz2][keep_mz2],
                           intensity = intensity2[remaining_mz2][keep_mz2])
      matched1 <- rbind(matched1, matched_nl1)
      matched2 <- rbind(matched2, matched_nl2)
    }
  }
  ## Store the m/z weighted intensity sum as an attribute
  attributes(matched1)$wintensity_sum <- sum(((intensity1^n) * (mz1 ^m))^2)
  attributes(matched2)$wintensity_sum <- sum(((intensity2^n) * (mz2 ^m))^2)
  return(list(x = matched1, y = matched2))
}

count_common_peaks <- function(x, y, xPrecursorMz, yPrecursorMz, ...) {
  # Extract m/z and intensity values from the spectra
  mz1 <- mz(x)[[1]]  # Extract the numeric vector from the list
  intensity1 <- intensity(x)[[1]]
  mz2 <- mz(y)[[1]]
  intensity2 <- intensity(y)[[1]]
  
  # Check for empty spectra
  if (length(mz1) == 0 || length(mz2) == 0)
    return(0)
  
  # Remove duplicates
  cleaned1 <- remove_duplicates(mz1, intensity1)
  cleaned2 <- remove_duplicates(mz2, intensity2)
  
  mz1 <- cleaned1$mz
  intensity1 <- cleaned1$intensity
  mz2 <- cleaned2$mz
  intensity2 <- cleaned2$intensity
  
  # Find common m/z values
  common_mz <- intersect(mz1, mz2)
  remaining_idx1 <- rep(TRUE, length(mz1))
  remaining_idx2 <- rep(TRUE, length(mz2))
  
  if (length(common_mz) > 0) {
    idx1 <- match(common_mz, mz1)
    idx2 <- match(common_mz, mz2)
    
    valid_indices <- !is.na(idx1) & !is.na(idx2)
    idx1 <- idx1[valid_indices]
    idx2 <- idx2[valid_indices]
    
    remaining_idx1[idx1] <- FALSE
    remaining_idx2[idx2] <- FALSE
  }
  
  common_peaks <- length(common_mz)
  
  # Calculate neutral losses only on remaining peaks
  nl1 <- neutral_loss(mz1[remaining_idx1], xPrecursorMz)
  nl2 <- neutral_loss(mz2[remaining_idx2], yPrecursorMz)
  
  common_nl <- intersect(nl1, nl2)
  
  # Add common neutral losses to common peaks count
  common_peaks <- common_peaks + length(common_nl)
  
  return(common_peaks)
}



## To use these functions with `matchSpectra()`:
## Use the `dynlib_map` function as `MAPFUN` and the
## `dynlib_symmetric_dotproduct` function as `FUN`.

#csp <- CompareSpectraParam(
# ppm = 0,
#tolerance = 0.005,
#threshold = 0.8,
#requirePrecursor = TRUE,
#MAPFUN = dynlib_map,
#FUN = dynlib_symmetric_dotproduct
#)

#res <- matchSpectra(query = xxx, target = yyy, param = csp)

library(BiocParallel)
library(Spectra)


sim_final_final <- function(st_sps, dy_sps, polarity_query,
                            polarity_target, machine_target,
                            threshold = 0.8, ppm = 0, tolerance = 0.005) {
  
  library(MetaboAnnotation)
  #Compare just with flax seed's compound that match the following criteria
  st_filtered <- st_sps[st_sps$polarity == polarity_query ]
  
  #st_filtered$precursorMz[st_filtered$msLevel %in% c(3, 4, 5)] <- round_perl(st_filtered$precursorMz[st_filtered$msLevel %in% c(3, 4, 5)])
  
  #Compare just with Dynlib's compound that match the following criteria
  
  dy_filtered <- dy_sps[dy_sps$polarity == polarity_target & 
                          dy_sps$machine == machine_target & !is.na(dy_sps$name) 
                        & !(dy_sps$name == "") & !(grepl("^!", dy_sps$name))]
  
  if (length(st_filtered) == 0 || length(dy_filtered) == 0) {
    return(data.frame())
  }
  
  # Paramètres de comparaison
  param <- CompareSpectraParam(
    ppm = ppm,
    tolerance = tolerance,
    threshold = threshold,
    requirePrecursor = TRUE,
    MAPFUN = dynlib_map,
    FUN = dynlib_symmetric_dotproduct,
    matchedPeaksCount  = TRUE
  )
  
  # Calcul de similarité
  matches <- matchSpectra(
    query = st_filtered,
    target = dy_filtered,
    param = param
  )
  
  # Extraction et filtrage des résultats
  df <- matchedData(matches)
  df <- df[!is.na(df$score) & df$score >= threshold, ]
  
  return(df)
}

