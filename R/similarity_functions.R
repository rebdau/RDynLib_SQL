library(Spectra)

neutral_loss <- function(mz_values, precursorMz,
                         method = c("unit.resolution", "high.resolution"),
                         digits = 4) {
  
  method <- match.arg(method)
  
  round_perl(precursorMz, method = method, digits = digits) - mz_values
}


round_perl <- function(number,
                       method = c("unit.resolution", "high.resolution"),
                       digits = 4) {
  
  method <- match.arg(method)
  
  if (method == "unit.resolution") {
    floor(number + 0.5)
  } else {
    round(number, digits)
  }
}


remove_duplicates <- function(mz, intensity,
                              method = c("unit.resolution", "high.resolution"),
                              digits = 4) {
  
  method <- match.arg(method)
  
  mz <- round_perl(mz, method = method, digits = digits)
  
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
  
  mz1 <- x[, 1]
  intensity1 <- x[, 2]
  mz2 <- y[, 1]
  intensity2 <- y[, 2]
  
  intensity1 <- (intensity1^n) * (mz1^m)
  intensity2 <- (intensity2^n) * (mz2^m)
  
  sum(intensity1 * intensity2)^2 /
    (attr(x, "wintensity_sum") * attr(y, "wintensity_sum"))
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
dynlib_map <- function(x, y,
                       xPrecursorMz, yPrecursorMz,
                       n = 3, m = 0.6,
                       fragments_method = c("unit.resolution", "high.resolution"),
                       precursor_method = c("unit.resolution", "high.resolution"),
                       digits = 4, ...) {
  
  fragments_method <- match.arg(fragments_method)
  precursor_method <- match.arg(precursor_method)
  
  if (!nrow(x) || !nrow(y)) {
    return(list(
      x = matrix(numeric(), ncol = 2, nrow = 0),
      y = matrix(numeric(), ncol = 2, nrow = 0)
    ))
  }
  
  ##Fragment resolution (mz peaks)
  cleaned1 <- remove_duplicates(
    x[, 1], x[, 2],
    method = fragments_method,
    digits = digits
  )
  
  cleaned2 <- remove_duplicates(
    y[, 1], y[, 2],
    method = fragments_method,
    digits = digits
  )
  
  mz1 <- cleaned1$mz
  intensity1 <- cleaned1$intensity
  mz2 <- cleaned2$mz
  intensity2 <- cleaned2$intensity
  
  ## Exact fragment matching
  common_mz <- intersect(mz1, mz2)
  
  keep_mz1 <- mz1 %in% common_mz
  keep_mz2 <- mz2 %in% common_mz
  
  matched1 <- cbind(mz = mz1[keep_mz1],
                    intensity = intensity1[keep_mz1])
  
  matched2 <- cbind(mz = mz2[keep_mz2],
                    intensity = intensity2[keep_mz2])
  
  ##Neutral loss (precursor resolution) 
  remaining_mz1 <- !keep_mz1
  remaining_mz2 <- !keep_mz2
  
  if (any(remaining_mz1) && any(remaining_mz2)) {
    
    nl1 <- neutral_loss(
      mz1[remaining_mz1],
      xPrecursorMz,
      method = precursor_method,
      digits = digits
    )
    
    nl2 <- neutral_loss(
      mz2[remaining_mz2],
      yPrecursorMz,
      method = precursor_method,
      digits = digits
    )
    
    common_nl <- intersect(nl1, nl2)
    
    if (length(common_nl)) {
      
      keep_nl1 <- nl1 %in% common_nl
      keep_nl2 <- nl2 %in% common_nl
      
      matched_nl1 <- cbind(
        mz = mz1[remaining_mz1][keep_nl1],
        intensity = intensity1[remaining_mz1][keep_nl1]
      )
      
      matched_nl2 <- cbind(
        mz = mz2[remaining_mz2][keep_nl2],
        intensity = intensity2[remaining_mz2][keep_nl2]
      )
      
      matched1 <- rbind(matched1, matched_nl1)
      matched2 <- rbind(matched2, matched_nl2)
    }
  }
  
  ## Weighted intensity sums (always computed on cleaned peaks)
  attr(matched1, "wintensity_sum") <-
    sum(((intensity1^n) * (mz1^m))^2)
  
  attr(matched2, "wintensity_sum") <-
    sum(((intensity2^n) * (mz2^m))^2)
  
  list(x = matched1, y = matched2)
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



symmetric_dotproduct_combined <- function(
    st_sps, dy_sps,
    polarity_query, polarity_target,
    precursor_resolution = c("unit.resolution", "high.resolution"),
    fragments_resolution = c("unit.resolution", "high.resolution"),
    threshold = 0.8,
    ppm = 0,
    tolerance = 0.005,
    digits = 4
) {
  
  if (!requireNamespace("rstudioapi", quietly = TRUE)) {
    install.packages("rstudioapi", repos = "https://cloud.r-project.org/")
  }
  library(rstudioapi)
  
  precursor_resolution <- match.arg(precursor_resolution)
  fragments_resolution <- match.arg(fragments_resolution)
  
  st_filtered <- st_sps[st_sps$polarity == polarity_query]
  dy_filtered <- dy_sps[
    dy_sps$polarity == polarity_target &
      !is.na(dy_sps$name) & dy_sps$name != "" &
      !grepl("^!", dy_sps$name)
  ]
  
  if (length(st_filtered) == 0 || length(dy_filtered) == 0)
    return(data.frame())
  
  # MS LEVEL prompt
  if ("msLevel" %in% names(spectraData(st_filtered))) {
    available_ms_st <- sort(unique(spectraData(st_filtered)$msLevel))
    cat("\nAvailable MS levels for st_sps:\n")
    print(available_ms_st)
    
    repeat {
      ms_input_st <- rstudioapi::showPrompt(
        "MS level st_sps",
        paste("Available MS levels:", paste(available_ms_st, collapse = ", "), "\nEnter MS level(s) for st_sps (comma separated, e.g. 2 or 2,3):")
      )
      if (is.null(ms_input_st)) {
        cat("Operation cancelled by user.\n")
        return(data.frame())
      }
      chosen_ms_st <- suppressWarnings(as.numeric(unlist(strsplit(ms_input_st, ","))))
      if (!any(is.na(chosen_ms_st)) && all(chosen_ms_st %in% available_ms_st)) break
      cat("Invalid input. Please choose from:", paste(available_ms_st, collapse = ","), "\n")
    }
    
    st_filtered <- st_filtered[spectraData(st_filtered)$msLevel %in% chosen_ms_st]
  }
  
  if ("msLevel" %in% names(spectraData(dy_filtered))) {
    available_ms_dy <- sort(unique(spectraData(dy_filtered)$msLevel))
    cat("\nAvailable MS levels for dy_sps:\n")
    print(available_ms_dy)
    
    repeat {
      ms_input_dy <- rstudioapi::showPrompt(
        "MS level dy_sps",
        paste("Available MS levels:", paste(available_ms_dy, collapse = ", "), "\nEnter MS level(s) for dy_sps (comma separated):")
      )
      if (is.null(ms_input_dy)) {
        cat("Operation cancelled by user.\n")
        return(data.frame())
      }
      chosen_ms_dy <- suppressWarnings(as.numeric(unlist(strsplit(ms_input_dy, ","))))
      if (!any(is.na(chosen_ms_dy)) && all(chosen_ms_dy %in% available_ms_dy)) break
      cat("Invalid input. Please choose from:", paste(available_ms_dy, collapse = ","), "\n")
    }
    
    dy_filtered <- dy_filtered[spectraData(dy_filtered)$msLevel %in% chosen_ms_dy]
  }
  
  if (length(st_filtered) == 0 || length(dy_filtered) == 0)
    return(data.frame())
  
  # Spectrum.type prompt
  if ("spectrum.type" %in% names(spectraData(st_filtered))) {
    available_type_st <- sort(unique(na.omit(spectraData(st_filtered)$spectrum.type)))
    cat("\nAvailable spectrum types for st_sps:\n")
    print(available_type_st)
    
    repeat {
      type_input_st <- rstudioapi::showPrompt(
        "Spectrum type query",
        paste("Available types:", paste(available_type_st, collapse = ", "), "\nEnter spectrum.type(s) for the query (comma separated if multiple):")
      )
      if (is.null(type_input_st)) {
        cat("Operation cancelled by user.\n")
        return(data.frame())
      }
      
      chosen_type_st <- unlist(strsplit(type_input_st, ","))
      chosen_type_st <- trimws(chosen_type_st)
      
      # Debug: Afficher les valeurs saisies et disponibles
      cat("You entered:", paste(chosen_type_st, collapse = ", "), "\n")
      cat("Available types:", paste(available_type_st, collapse = ", "), "\n")
      
      # Vérifier si toutes les entrées sont valides
      if (all(chosen_type_st %in% available_type_st)) break
      
      cat("Invalid input. Please choose from:", paste(available_type_st, collapse = ", "), "\n")
    }
    
    st_filtered <- st_filtered[spectraData(st_filtered)$spectrum.type %in% chosen_type_st]
  }
  
  if ("spectrum.type" %in% names(spectraData(dy_filtered))) {
    available_type_dy <- sort(unique(na.omit(spectraData(dy_filtered)$spectrum.type)))
    cat("\nAvailable spectrum types for the target:\n")
    print(available_type_dy)
    
    repeat {
      type_input_dy <- rstudioapi::showPrompt(
        "Spectrum type target",
        paste("Available types:", paste(available_type_dy, collapse = ", "), "\nEnter spectrum.type(s) for target (comma separated if multiple):")
      )
      if (is.null(type_input_dy)) {
        cat("Operation cancelled by user.\n")
        return(data.frame())
      }
      
      chosen_type_dy <- unlist(strsplit(type_input_dy, ","))
      chosen_type_dy <- trimws(chosen_type_dy)
      
      # Debug: Afficher les valeurs saisies et disponibles
      cat("You entered:", paste(chosen_type_dy, collapse = ", "), "\n")
      cat("Available types:", paste(available_type_dy, collapse = ", "), "\n")
      
      if (all(chosen_type_dy %in% available_type_dy)) break
      
      cat("Invalid input. Please choose from:", paste(available_type_dy, collapse = ", "), "\n")
    }
    
    dy_filtered <- dy_filtered[spectraData(dy_filtered)$spectrum.type %in% chosen_type_dy]
  }
  
  
  if (length(st_filtered) == 0 || length(dy_filtered) == 0)
    return(data.frame())
  
  # Similarity calculation
  param <- MetaboAnnotation::CompareSpectraParam(
    ppm = ppm,
    tolerance = tolerance,
    threshold = threshold,
    requirePrecursor = TRUE,
    MAPFUN = dynlib_map,
    FUN = dynlib_symmetric_dotproduct,
    matchedPeaksCount = TRUE,
    fragments_method = fragments_resolution,
    precursor_method = precursor_resolution,
    digits = digits
  )
  
  matches <- MetaboAnnotation::matchSpectra(
    query = st_filtered,
    target = dy_filtered,
    param = param
  )
  
  df <- as.data.frame(MetaboAnnotation::matchedData(matches))
  df <- df[!is.na(df$score) & df$score >= threshold, ]
  
  if (nrow(df) == 0) return(df)
  
  df <- dplyr::group_by(df, acquisitionNum)
  df <- dplyr::slice_max(df, score, n = 1, with_ties = TRUE)
  df <- dplyr::slice_max(df, matched_peaks_count, n = 1, with_ties = FALSE)
  df <- dplyr::ungroup(df)
  
  return(df)
}







