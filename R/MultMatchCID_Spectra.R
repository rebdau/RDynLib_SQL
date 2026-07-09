#' Match a query spectrum against spectra in a DynLib database
#'
#' Calculates spectral similarity between a query spectrum, selected by
#' `compound_id`, and candidate spectra in a DynLib `Spectra` database.
#' The function supports both high-resolution and unit-resolution MS/MS data
#' and applies resolution-specific peak matching and similarity scoring
#' methods.
#' 
#' @param spectra_db `Spectra` object, represents the whole database. 
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
  
  if (!all(ms_level_query %in% available_ms))
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
  
  if (!all(ms_level_target %in% available_target_ms))
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