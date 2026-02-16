############################################################
# Load required libraries
############################################################
library(Spectra)
library(BiocParallel)
library(dplyr)
library(MetaboAnnotation)

############################################################
# Rounding strategy
############################################################
mz_rounding_strategy <- function(method = c("low.resolution", "high.resolution"),
                                 digits = 2) {
  method <- match.arg(method)
  if (method == "low.resolution") {
    function(x) floor(x + 0.5)
  } else {
    function(x) round(x, digits)
  }
}

############################################################
# Neutral loss
############################################################
neutral_loss <- function(mz_values, precursorMz, round_fun) {
  round_fun(precursorMz) - mz_values
}

############################################################
# Remove duplicates
############################################################
remove_duplicates <- function(mz, intensity, round_fun) {
  mz <- round_fun(mz)
  unique_mz <- unique(mz)
  list(
    mz = unique_mz,
    intensity = vapply(
      unique_mz,
      function(m) max(intensity[mz == m], na.rm = TRUE),
      NA_real_
    )
  )
}

############################################################
# Symmetric dot product
############################################################
dynlib_symmetric_dotproduct <- function(x, y, n = 3, m = 0.6, ...) {
  if (!nrow(x) || !nrow(y)) return(0.0)
  mz1 <- x[, 1L]; intensity1 <- x[, 2L]
  mz2 <- y[, 1L]; intensity2 <- y[, 2L]
  intensity1 <- (intensity1^n) * (mz1^m)
  intensity2 <- (intensity2^n) * (mz2^m)
  sum(intensity1 * intensity2)^2 /
    (attributes(x)$wintensity_sum * attributes(y)$wintensity_sum)
}

############################################################
# dynlib_map
############################################################
dynlib_map <- function(x, y,
                       xPrecursorMz,
                       yPrecursorMz,
                       resolution = "high.resolution",
                       digits = 2,
                       n = 3,
                       m = 0.6,
                       ...) {
  if (!nrow(x) || !nrow(y))
    return(list(x = matrix(numeric(), ncol = 2, nrow = 0),
                y = matrix(numeric(), ncol = 2, nrow = 0)))
  
  round_fun <- mz_rounding_strategy(resolution, digits)
  
  cleaned1 <- remove_duplicates(x[, 1L], x[, 2L], round_fun)
  cleaned2 <- remove_duplicates(y[, 1L], y[, 2L], round_fun)
  
  mz1 <- cleaned1$mz; intensity1 <- cleaned1$intensity
  mz2 <- cleaned2$mz; intensity2 <- cleaned2$intensity
  
  # Matching exact m/z
  common_mz <- intersect(mz1, mz2)
  keep_mz1 <- mz1 %in% common_mz
  keep_mz2 <- mz2 %in% common_mz
  
  matched1 <- cbind(mz = mz1[keep_mz1], intensity = intensity1[keep_mz1])
  matched2 <- cbind(mz = mz2[keep_mz2], intensity = intensity2[keep_mz2])
  
  # Neutral loss matching
  remaining_mz1 <- !keep_mz1
  remaining_mz2 <- !keep_mz2
  if (any(remaining_mz1) && any(remaining_mz2)) {
    nl1 <- neutral_loss(mz1[remaining_mz1], xPrecursorMz, round_fun)
    nl2 <- neutral_loss(mz2[remaining_mz2], yPrecursorMz, round_fun)
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
  
  attributes(matched1)$wintensity_sum <- sum(((intensity1^n) * (mz1^m))^2)
  attributes(matched2)$wintensity_sum <- sum(((intensity2^n) * (mz2^m))^2)
  
  return(list(x = matched1, y = matched2))
}

############################################################
# Count common peaks
############################################################
count_common_peaks <- function(x, y,
                               xPrecursorMz,
                               yPrecursorMz,
                               resolution = "high.resolution",
                               digits = 2,
                               ...) {
  round_fun <- mz_rounding_strategy(resolution, digits)
  
  mz1 <- mz(x)[[1]]; intensity1 <- intensity(x)[[1]]
  mz2 <- mz(y)[[1]]; intensity2 <- intensity(y)[[1]]
  
  if (length(mz1) == 0 || length(mz2) == 0) return(0)
  
  cleaned1 <- remove_duplicates(mz1, intensity1, round_fun)
  cleaned2 <- remove_duplicates(mz2, intensity2, round_fun)
  
  mz1 <- cleaned1$mz; mz2 <- cleaned2$mz
  common_mz <- intersect(mz1, mz2)
  remaining_idx1 <- !mz1 %in% common_mz
  remaining_idx2 <- !mz2 %in% common_mz
  
  common_peaks <- length(common_mz)
  nl1 <- neutral_loss(mz1[remaining_idx1], xPrecursorMz, round_fun)
  nl2 <- neutral_loss(mz2[remaining_idx2], yPrecursorMz, round_fun)
  common_peaks <- common_peaks + length(intersect(nl1, nl2))
  return(common_peaks)
}

############################################################
# Interactive compare_spectra_filtered
############################################################
#' Compare spectra with interactive selection of collision energy
#' and spectrum type based on dataset values.
#'
#' @param st_sps Spectra object of sample spectra
#' @param dy_sps Spectra object of DynLib reference spectra
#' @param polarity_query Polarity of sample spectra to compare
#' @param polarity_target Polarity of reference spectra to compare
#' @param machine_target Machine type of reference spectra
#' @param threshold Similarity threshold
#' @param ppm PPM tolerance
#' @param tolerance Absolute tolerance
#' @param mz_method m/z rounding method
#' @param digits Number of digits if high resolution
#' @param collision_energy Optional user selection of collision energy
#' @param spectrum_type Optional user selection of spectrum type
#' @return Data.frame of matched spectra, slice-maxed
#' @export
compare_spectra_filtered <- function(
    st_sps, dy_sps,
    polarity_query,
    polarity_target,
    machine_target,
    threshold = 0.8,
    ppm = 0,
    tolerance = 0.005,
    mz_method = c("low.resolution", "high.resolution"),
    digits = 2,
    collision_energy = NULL,
    spectrum_type = NULL
) {
  # Rounding function
  round_fun <- mz_rounding_strategy(mz_method, digits)
  
  # Filter query spectra by polarity
  st_filtered <- st_sps[st_sps$polarity == polarity_query]
  
  # If multiple collision energies exist, ask the user
  if ("collision_energy" %in% names(spectraData(st_filtered))) {
    available_ce <- unique(spectraData(st_filtered)$collision_energy)
    if (length(available_ce) > 1 && is.null(collision_energy)) {
      cat("Available collision energies:\n")
      print(available_ce)
      collision_energy <- readline(prompt = "Select collision energy: ")
      collision_energy <- as.numeric(collision_energy)
    }
    st_filtered <- st_filtered[spectraData(st_filtered)$collision_energy == collision_energy]
  }
  
  # If multiple spectrum types exist, ask the user
  if ("spectrum_type" %in% names(spectraData(st_filtered))) {
    available_types <- unique(spectraData(st_filtered)$spectrum_type)
    if (length(available_types) > 1 && is.null(spectrum_type)) {
      cat("Available spectrum types:\n")
      print(available_types)
      spectrum_type <- readline(prompt = "Select spectrum type: ")
    }
    st_filtered <- st_filtered[spectraData(st_filtered)$spectrum_type %in% spectrum_type]
  }
  
  # Filter reference spectra
  dy_filtered <- dy_sps[
    dy_sps$polarity == polarity_target &
      dy_sps$machine == machine_target &
      !is.na(dy_sps$name) &
      dy_sps$name != "" &
      !grepl("^!", dy_sps$name)
  ]
  
  if (length(st_filtered) == 0 || length(dy_filtered) == 0) return(data.frame())
  
  # matchSpectra parameters
  param <- CompareSpectraParam(
    ppm = ppm,
    tolerance = tolerance,
    threshold = threshold,
    requirePrecursor = TRUE,
    MAPFUN = dynlib_map,
    FUN = dynlib_symmetric_dotproduct,
    matchedPeaksCount = TRUE
  )
  
  # Compute similarity
  matches <- matchSpectra(query = st_filtered, target = dy_filtered, param = param)
  
  # Extract results
  df <- as.data.frame(matchedData(matches))
  for (col in names(df)) {
    if (inherits(df[[col]], "Rle")) df[[col]] <- as.vector(df[[col]])
  }
  
  # Apply threshold
  df <- df[!is.na(df$score) & df$score >= threshold, ]
  if (nrow(df) == 0) return(df)
  
  # Slice-max by score then matched peaks
  df <- df %>%
    group_by(acquisitionNum) %>%
    slice_max(order_by = score, n = 1, with_ties = TRUE) %>%
    slice_max(order_by = matched_peaks_count, n = 1, with_ties = FALSE) %>%
    ungroup()
  
  return(df)
}
