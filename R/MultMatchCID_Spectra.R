#' @description
#'
#' Helper function to round the precursorMz values depending on the resolution  
#' type chosen by the user.
#' 
#' @return 
#' a rounded precursorMz value.
#' 
#' @author Ahlam Mentag
neutral_loss <- function(mz_values, precursorMz) {
  
  
  
  round_perl(precursorMz- mz_values)
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
round_perl <- function(number) {
  
  
  floor(number + 0.5)
}

#' @description
#'
#' Helper function to remove the duplicate mz values after the rounding.
#' 
#' @return 
#' a list of mz and intensity peaks after removing all the duplicates.
#' 
#' @author Ahlam Mentag
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
#' Main function to calculate the spectral similarity of a compound to other 
#' compounds in the same database (spectra object) produced from a high 
#' or unit resolution instrument.
#' 
#' @param spectra_db spectra object, represents the whole database. 
#'
#' @param compound_id compound_id to compare with the compounds from spectra_db.  
#' 
#' @param polarity_query 'numeric' the query polarity, 0  if negative,
#'         and 1 for positive polarity.
#'        
#' @param polarity_target 'numeric' the target polarity, 0  if negative,
#'         and 1 for positive polarity.
#'        
#'        
#' @param fragments_resolution 'character(1)' the resolution type of the fragments  
#'        masses, it is used to define the rounding type of the fragments,
#'        it could be either "unit.resolution" or "high.resolution".
#'        
#' @param requirePrecursor 'logical(1)' by default true, it allows to pre-filter
#'        the target spectra prior to the actual similarity calculation for 
#'        each individual query spectrum. 
#'        
#' @param threshold 'numeric(1)' the score threshold to select a best match, by 
#'        default it is 0.8.
#'       
#' @param ppm 'numeric(1)' the acceptable difference between m/z values of the
#'        compared peaks:
#'        - For high resolution : ppm is applied for matching the product ions 
#'          and match precursorMz if requirePrecursor is TRUE.
#'        - For unit resolution : ppm is used for matching precursorMz. 
#'        
#' @param tolerance 'numeric(1)' the acceptable difference between m/z values
#'         of the compared peaks.
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
MultMatchCID_Spectra <- function(
    spectra_db,
    compound_id,
    polarity_query,
    polarity_target,
    spectrum_type_query = NULL,
    spectrum_type_target = NULL,
    ms_level_query = NULL,
    ms_level_target = NULL,
    minIons = 3,
    fragments_resolution = c("unit.resolution","high.resolution"),
    requirePrecursor = TRUE,
    threshold = 0.8,
    ppm = 2,
    tolerance = 0.5
) {
  
  fragments_resolution <- match.arg(fragments_resolution)
  
  
  # FILTER QUERY
  
  query_meta_cols <- intersect(
    c("compound_id","polarity","name","spectrum_type","spectrum.type",
      "msLevel","precursorMz","rtime","expid", "retention_time"),
    spectraVariables(spectra_db)
  )
  
  meta_query <- spectraData(spectra_db, query_meta_cols)
  
  query_sps <- spectra_db[
    meta_query$compound_id == compound_id &
      meta_query$polarity == polarity_query
  ]
  
  if (length(query_sps) == 0)
    stop("Compound not found with given polarity.")
  
  
  # Detect spectrum type column
  spectrum_col <- intersect(c("spectrum_type","spectrum.type"), query_meta_cols)
  spectrum_col <- if(length(spectrum_col) > 0) spectrum_col[1] else NULL
  
  
  # Select only needed variables
  vars_keep <- c("compound_id","polarity","msLevel", "name", "rtime","precursorMz","expid", "retention_time")
  if (!is.null(spectrum_col))
    vars_keep <- c(spectrum_col, vars_keep)
  
  query_sps <- selectSpectraVariables(query_sps, vars_keep)
  
  
  # spectrum_type filter
  
  if (!is.null(spectrum_col)) {
    
    available_types <- unique(spectraData(query_sps, spectrum_col)[[1]])
    
    if (is.null(spectrum_type_query)) {
      
      if (length(available_types) > 1) {
        spectrum_type_query <- rstudioapi::showPrompt(
          "Spectrum type selection (query)",
          paste("Available types:", paste(available_types, collapse = ", ")),
          default = available_types[1]
        )
      }
      
    }
    
    if (!is.null(spectrum_type_query) && spectrum_type_query %in% available_types) {
      
      query_sps <- query_sps[
        spectraData(query_sps, spectrum_col)[[1]] == spectrum_type_query
      ]
      
    }
    
  }
  
  
  # MS level query filter
  
  available_ms <- unique(msLevel(query_sps))
  
  if (is.null(ms_level_query)) {
    
    ms_level_query <- as.numeric(
      rstudioapi::showPrompt(
        "MS level selection (query)",
        paste("Available MS levels:", paste(available_ms, collapse = ", ")),
        default = available_ms[1]
      )
    )
    
  }
  
  if (!(ms_level_query %in% available_ms))
    stop("Invalid query MS level selected.")
  
  query_sps <- query_sps[msLevel(query_sps) == ms_level_query]
  
  if (length(query_sps) == 0)
    stop("No query spectra left after filtering.")
  
  
  query_sp <- query_sps[1]
  query_peaks <- peaksData(query_sp)[[1]]
  
  if (nrow(query_peaks) < minIons)
    stop("Not enough product ions in query spectrum.")
  
  
  
  # FILTER TARGET
  
  target_meta_cols <- intersect(
    c("polarity","name","spectrum_type","spectrum.type",
      "msLevel","precursorMz","rtime","expid", "retention_time"),
    spectraVariables(spectra_db)
  )
  
  meta_target <- spectraData(spectra_db, target_meta_cols)
  
  target_sps <- spectra_db[
    meta_target$polarity == polarity_target &
      !is.na(meta_target$name) &
      meta_target$name != "" &
      !grepl("^!", meta_target$name)
  ]
  
  
  vars_keep <- c("compound_id","polarity","msLevel","rtime", "name", "precursorMz","expid", "retention_time")
  if (!is.null(spectrum_col))
    vars_keep <- c(spectrum_col, vars_keep)
  
  target_sps <- selectSpectraVariables(target_sps, vars_keep)
  
  
  
  # spectrum_type filter target
  if (!is.null(spectrum_col)) {
    
    available_types_target <- unique(spectraData(target_sps, spectrum_col)[[1]])
    
    # Convert NA to a string for display
    display_types <- ifelse(is.na(available_types_target), "NA", available_types_target)
    
    if (is.null(spectrum_type_target)) {
      
      cat("\nAvailable target spectrum types:\n")
      print(display_types)
      
      spectrum_type_target <- rstudioapi::showPrompt(
        "Spectrum type selection (target)",
        "Enter one of the available spectrum types:",
        default = display_types[1]
      )
      
    }
    
    # Convert back if user selected "NA"
    if (identical(spectrum_type_target, "NA"))
      spectrum_type_target <- NA
    
    target_sps <- target_sps[
      spectraData(target_sps, spectrum_col)[[1]] == spectrum_type_target |
        (is.na(spectraData(target_sps, spectrum_col)[[1]]) & is.na(spectrum_type_target))
    ]
    
  }
  
  # MS level target filter
  
  available_target_ms <- unique(msLevel(target_sps))
  
  if (is.null(ms_level_target)) {
    
    ms_level_target <- as.numeric(
      rstudioapi::showPrompt(
        "MS level selection (target)",
        paste("Available MS levels:", paste(available_target_ms, collapse = ", ")),
        default = available_target_ms[1]
      )
    )
    
  }
  
  if (!(ms_level_target %in% available_target_ms))
    stop("Invalid target MS level selected.")
  
  target_sps <- target_sps[msLevel(target_sps) == ms_level_target]
  
  
  target_sps <- target_sps[
    sapply(peaksData(target_sps), nrow) >= minIons
  ]
  
  if (length(target_sps) == 0)
    return("No target spectra left after filtering.")
  
  
  
  # MATCHING
  
  if (fragments_resolution == "unit.resolution") {
    
    param <- MetaboAnnotation::CompareSpectraParam(
      ppm = ppm,
      tolerance = tolerance,
      threshold = threshold,
      requirePrecursor = requirePrecursor,
      MAPFUN = dynlibmatch2_map,
      FUN = dynlib_symmetric_dotproduct,
      matchedPeaksCount = TRUE,
      fragments_method = fragments_resolution
    )
    
    matches <- MetaboAnnotation::matchSpectra(
      query = query_sp,
      target = target_sps,
      param = param
    )
    
  } else {
    
    param <- CompareSpectraParam(
      ppm = ppm,
      tolerance = tolerance,
      MAPFUN = joinPeaksGnps,
      FUN = MsCoreUtils::gnps,
      threshold = threshold,
      requirePrecursor = requirePrecursor,
      matchedPeaksCount = TRUE
    )
    
    matches <- matchSpectra(
      query = query_sp,
      target = target_sps,
      param = param
    )
    
  }
  
  
  
  # RESULTS

  
  df <- as.data.frame(MetaboAnnotation::matchedData(matches))
  df <- df[!is.na(df$score) & df$score >= threshold, ]
  
  if (nrow(df) == 0)
    return("No matches found.")
  
  df <- df[df$matched_peaks_count >= minIons, ]
  
  if (nrow(df) == 0)
    return("No matches found after minIons filtering.")
  
  
  if ("target_compound_id" %in% colnames(df)) df$COMPID <- df$target_compound_id
  if ("target_precursorMz" %in% colnames(df)) df$target_precursorMz <- df$target_precursorMz
  if ("precursorMz" %in% colnames(df)) df$query_precursorMz <- df$precursorMz
  
  df$`#FragIon` <- nrow(query_peaks)
  df$ComIons <- df$matched_peaks_count
  df$DotIons <- df$score

  
  if ("target_retention_time" %in% colnames(df)) df$target_retention_time <- df$target_retention_time
  if ("retention_time" %in% colnames(df)) df$query_retention_time <- df$retention_time
  if ("target_expid" %in% colnames(df)) df$target_expid <- df$target_expid
  if ("expid" %in% colnames(df)) df$query_expid <- df$expid
  if ("target_msLevel" %in% colnames(df)) df$target_msLevel <- df$target_msLevel
  if ("msLevel" %in% colnames(df)) df$query_msLevel <- df$msLevel
  if ("target_name" %in% colnames(df)) df$target_name <- df$target_name
  
  cols_order <- c("target_compound_id", "target_msLevel", "query_msLevel", "query_precursorMz", "target_precursorMz","#FragIon",
                  "ComIons","DotIons","target_name","target_retention_time","query_retention_time","target_expid",
                  "query_expid")
  cols_order <- cols_order[cols_order %in% colnames(df)]
  
  df <- df[, cols_order, drop = FALSE]

  # FORMAT OUTPUT
  
  # Round score
  if ("DotIons" %in% colnames(df)) {
    df$DotIons <- round(df$DotIons, 2)
  }
  
  # Round precursor m/z if unit resolution
  if (fragments_resolution == "unit.resolution") {
    
    if ("query_precursorMz" %in% colnames(df)) {
      df$query_precursorMz <- round(df$query_precursorMz,2)
    }
    
    if ("target_precursorMz" %in% colnames(df)) {
      df$target_precursorMz <- round(df$target_precursorMz,2)
    }
    
  }
  
  
  return(df)
  
}