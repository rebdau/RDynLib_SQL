#' @description
#'
#' Helper function to round the precursorMz values depending on the resolution  
#' type chosen by the user.
#' 
#' @return 
#' a rounded precursorMz value.
#' 
#' @author Ahlam Mentag
neutral_loss <- function(mz_values, precursorMz,
                         method = c("unit.resolution", "high.resolution"),
                         digits = 4) {
  
  method <- match.arg(method)
  
  round_perl(precursorMz, method = method, digits = digits - mz_values)
}

#' @description
#'
#' Helper function to round the mz values depending on the resolution type 
#' chosen by the user.
#' 
#' @return 
#' a rounded mz value.
#' 
#' @author Rebecca Dauwe
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

#' @description
#'
#' Helper function to remove the duplicate mz values after the rounding.
#' 
#' @return 
#' a list of mz and intensity peaks after removing all the duplicates.
#' 
#' @author Ahlam Mentag
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
dynlibmatch2_map <- function(
    x, y,
    xPrecursorMz, yPrecursorMz,
    n = 3, m = 0.6,
    fragments_method = c("unit.resolution", "high.resolution"),
    precursor_method = c("unit.resolution", "high.resolution"),
    digits = 4,
    tolerance = 0, ppm = 0,
    ...
) {
  
  fragments_method <- match.arg(fragments_method)
  precursor_method <- match.arg(precursor_method)
  
  ## Return empty aligned matrices if one spectrum has no peaks
  if (!nrow(x) || !nrow(y)) {
    empty <- matrix(numeric(), ncol = 2, nrow = 0)
    return(list(empty, empty))
  }
  

  ## Clean fragment peaks

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
      xPrecursorMz,
      method = precursor_method,
      digits = digits
    )
    
    nl2 <- neutral_loss(
      mz2[remaining2],
      yPrecursorMz,
      method = precursor_method,
      digits = digits
    )
    
    idx_nl2 <- match(nl1, nl2)
    valid_nl <- !is.na(idx_nl2)
    
    if (any(valid_nl)) {
      
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




#' @description
#'
#' Main function to calculate the spectral similarity between spectra produced 
#' from a high or unit resolution instrument.
#' 
#' @param st_sps spectra object, represents the query data 
#'
#' @param dy_sps spectra object, represents the target data 
#' 
#' @param polarity_query 'numeric' the query polarity, 0  if negative,
#'         and 1 for positive polarity
#'        
#' @param polarity_target 'numeric' the target polarity, 0  if negative,
#'         and 1 for positive polarity
#'        
#' @param precursor_resolution 'character' the resolution type of the precursor  
#'        masses, it is used to define the rounding type of the precursorMz(),
#'        it could be either "unit.resolution" or "high.resolution".
#'        
#' @param fragments_resolution 'character(1)' the resolution type of the fragments  
#'        masses, it is used to define the rounding type of the fragments,
#'        it could be either "unit.resolution" or "high.resolution".
#'        
#' @param threshold 'numeric(1)' the score threshold to select a best match, by 
#'        default it is 0.8.
#'       
#' @param digits 'numeric(1)' rounding number, by default it is 4.  
#'      
#' @return
#'
#' a matrix with the target and the query matches spectraData, the score, and 
#' the number of matched peaks.
#'
#' @import BiocParallel
#' @import Spectra
#' @import MetaboAnnotation
#' @import rstudioapi
#' 
#' @author Ahlam Mentag
#' 
#' @export
similarity_RDynlib_match <- function(
    st_sps, dy_sps,
    polarity_query, polarity_target,
    precursor_resolution = c("unit.resolution", "high.resolution"),
    fragments_resolution = c("unit.resolution", "high.resolution"),
    threshold = 0.8,
    digits = 4
) {
  
  if (!requireNamespace("rstudioapi", quietly = TRUE)) {
    install.packages("rstudioapi", repos = "https://cloud.r-project.org/")
  }
  library(rstudioapi)
  
  precursor_resolution <- match.arg(precursor_resolution)
  fragments_resolution <- match.arg(fragments_resolution)
  
  # polarity filter
  st_filtered <- st_sps[st_sps$polarity == polarity_query]
  
  dy_filtered <- dy_sps[
    dy_sps$polarity == polarity_target &
      !is.na(dy_sps$name) & dy_sps$name != "" &
      !grepl("^!", dy_sps$name)
  ]
  
  if (length(st_filtered) == 0 || length(dy_filtered) == 0)
    return(data.frame())
  
  
  # detect spectrum type column 
  detect_spectrum_type_column <- function(sps) {
    cols <- names(spectraData(sps))
    if ("spectrum.type" %in% cols) return("spectrum.type")
    if ("spectrum_type" %in% cols) return("spectrum_type")
    return(NULL)
  }
  
  

  interactive_filter <- function(sps, label) {
    
    # Spectrum type filtering 
    type_col <- detect_spectrum_type_column(sps)
    
    if (!is.null(type_col)) {
      
      available_types <- sort(unique(na.omit(
        spectraData(sps)[[type_col]]
      )))
      
      if (length(available_types) == 0)
        return(sps)
      
      cat("\nAvailable spectrum types for", label, ":\n")
      print(available_types)
      
      repeat {
        type_input <- rstudioapi::showPrompt(
          paste("Spectrum type -", label),
          paste("Available types:",
                paste(available_types, collapse = ", "))
        )
        
        if (is.null(type_input)) {
          cat("Operation cancelled by user.\n")
          return(NULL)
        }
        
        chosen_types <- trimws(unlist(strsplit(type_input, ",")))
        
        if (all(chosen_types %in% available_types)) break
        
        cat("Invalid input. Please choose from:",
            paste(available_types, collapse = ", "), "\n")
      }
      
      sps <- sps[
        spectraData(sps)[[type_col]] %in% chosen_types
      ]
    }
    
    
    if (length(sps) == 0)
      return(NULL)
    
    
    # MS level filtering after spectrum.type filtering 
    if ("msLevel" %in% names(spectraData(sps))) {
      
      available_ms <- sort(unique(na.omit(
        spectraData(sps)$msLevel
      )))
      
      if (length(available_ms) == 0)
        return(sps)
      
      cat("\nAvailable MS levels for selected spectrum.type (", 
          label, "):\n", sep = "")
      print(available_ms)
      
      repeat {
        ms_input <- rstudioapi::showPrompt(
          paste("MS level -", label),
          paste("Available MS levels:",
                paste(available_ms, collapse = ", "))
        )
        
        if (is.null(ms_input)) {
          cat("Operation cancelled by user.\n")
          return(NULL)
        }
        
        chosen_ms <- suppressWarnings(
          as.numeric(unlist(strsplit(ms_input, ",")))
        )
        
        chosen_ms <- chosen_ms[!is.na(chosen_ms)]
        
        if (length(chosen_ms) > 0 &&
            all(chosen_ms %in% available_ms)) break
        
        cat("Invalid input. Please choose from:",
            paste(available_ms, collapse = ", "), "\n")
      }
      
      sps <- sps[
        spectraData(sps)$msLevel %in% chosen_ms
      ]
    }
    
    return(sps)
  }
  
  
  # Apply interactive filtering 
  st_filtered <- interactive_filter(st_filtered, "query")
  if (is.null(st_filtered) || length(st_filtered) == 0)
    return(data.frame())
  
  dy_filtered <- interactive_filter(dy_filtered, "target")
  if (is.null(dy_filtered) || length(dy_filtered) == 0)
    return(data.frame())
  
  
  #Similarity calculation
  param <- MetaboAnnotation::CompareSpectraParam(

    threshold = threshold,
    requirePrecursor = TRUE,
    MAPFUN = dynlibmatch2_map,
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
  
  df <- as.data.frame(
    MetaboAnnotation::matchedData(matches)
  )
  print(environment(dynlib_map))
  df <- df[!is.na(df$score) & df$score >= threshold, ]
  
  if (nrow(df) == 0)
    return(df)
  
  df <- dplyr::group_by(df, acquisitionNum)
  df <- dplyr::slice_max(df, score, n = 1, with_ties = TRUE)
  df <- dplyr::slice_max(df, matched_peaks_count,
                         n = 1, with_ties = FALSE)
  df <- dplyr::ungroup(df)
  
  return(df)
}







