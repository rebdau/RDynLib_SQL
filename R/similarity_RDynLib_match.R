#' Search a DynLib spectral database for matches to query spectra
#'
#' Compares spectra from a query `Spectra` object, typically representing a
#' newly acquired experiment, against a target `Spectra` object derived from
#' a DynLib SQL database to identify spectral matches and transfer compound
#' annotations. The function supports both high-resolution and unit-resolution
#' mass spectrometry data and applies resolution-specific spectral matching
#' and similarity scoring methods.
#' 
#' @param st_sps `Spectra` object, represents the query data. 
#'
#' @param dy_sps `Spectra` object, represents the target data.
#' 
#' @param polarity_query 'numeric' the query polarity, 0  if negative,
#'         and 1 for positive polarity.
#'        
#' @param polarity_target 'numeric' the target polarity, 0  if negative,
#'         and 1 for positive polarity.
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
similarity_RDynLib_match <- function(
    st_sps, dy_sps,
    polarity_query, polarity_target,
    fragments_resolution = c("unit.resolution", "high.resolution"),
    requirePrecursor = TRUE,
    threshold = 0.8,
    ppm = 2,
    tolerance = 0.5
) {
  
  fragments_resolution <- match.arg(fragments_resolution)
  
  # Polarity filter
  pol_st <- spectraData(st_sps, "polarity")[["polarity"]]
  st_filtered <- st_sps[pol_st == polarity_query]
  
  dy_meta <- spectraData(dy_sps, c("polarity", "name"))
  
  keep_dy <- dy_meta$polarity == polarity_target &
    !is.na(dy_meta$name) &
    dy_meta$name != "" &
    !grepl("^!", dy_meta$name)
  
  dy_filtered <- dy_sps[keep_dy]
  
  if (length(st_filtered) == 0 || length(dy_filtered) == 0)
    return(data.frame())
  
  
  # Detect spectrum type column
  detect_spectrum_type_column <- function(sps) {
    
    cols <- spectraVariables(sps)
    
    if ("spectrum.type" %in% cols)
      return("spectrum.type")
    
    if ("spectrum_type" %in% cols)
      return("spectrum_type")
    
    return(NULL)
  }
  
  
  # Interactive filtering
  interactive_filter <- function(sps, label) {
    
    type_col <- detect_spectrum_type_column(sps)
    
    if (!is.null(type_col)) {
      
      type_data <- spectraData(sps, type_col)[[type_col]]
      
      available_types <- sort(unique(na.omit(type_data)))
      
      if (length(available_types) > 0) {
        
        repeat {
          
          type_input <- rstudioapi::showPrompt(
            paste("Spectrum type -", label),
            paste(
              "Available types:",
              paste(available_types, collapse = ", ")
            )
          )
          
          if (is.null(type_input))
            return(NULL)
          
          chosen_types <- trimws(
            unlist(strsplit(type_input, ","))
          )
          
          if (all(chosen_types %in% available_types))
            break
        }
        
        sps <- sps[type_data %in% chosen_types]
      }
    }
    
    
    if (length(sps) == 0)
      return(NULL)
    
    
    if ("msLevel" %in% spectraVariables(sps)) {
      
      ms_data <- spectraData(
        sps,
        "msLevel"
      )[["msLevel"]]
      
      available_ms <- sort(unique(na.omit(ms_data)))
      
      if (length(available_ms) > 0) {
        
        repeat {
          
          ms_input <- rstudioapi::showPrompt(
            paste("MS level -", label),
            paste(
              "Available MS levels:",
              paste(available_ms, collapse = ", ")
            )
          )
          
          if (is.null(ms_input))
            return(NULL)
          
          chosen_ms <- suppressWarnings(
            as.numeric(
              unlist(strsplit(ms_input, ","))
            )
          )
          
          chosen_ms <- chosen_ms[!is.na(chosen_ms)]
          
          if (length(chosen_ms) > 0 &&
              all(chosen_ms %in% available_ms))
            break
        }
        
        sps <- sps[ms_data %in% chosen_ms]
      }
    }
    
    sps
  }
  
  
  st_filtered <- interactive_filter(
    st_filtered,
    "query"
  )
  
  if (is.null(st_filtered) || length(st_filtered) == 0)
    return(data.frame())
  
  
  dy_filtered <- interactive_filter(
    dy_filtered,
    "target"
  )
  
  if (is.null(dy_filtered) || length(dy_filtered) == 0)
    return(data.frame())
  
  
  # Detect spectrum type columns
  query_type_col <- detect_spectrum_type_column(st_filtered)
  target_type_col <- detect_spectrum_type_column(dy_filtered)
  
  
  # Keep required target variables
  target_vars <- c(
    "msLevel",
    "rtime",
    "precursorMz",
    "name",
    "scanIndex",
    "dataOrigin",
    "compound_id"
  )
  
  
  if (!is.null(target_type_col))
    target_vars <- c(
      target_vars,
      target_type_col
    )
  
  
  dy_filtered <- selectSpectraVariables(
    dy_filtered,
    target_vars
  )
  
  
  # Matching
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
      query = st_filtered,
      target = dy_filtered,
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
      query = st_filtered,
      target = dy_filtered,
      param = param
    )
  }
  
  
  # Results
  df <- as.data.frame(
    MetaboAnnotation::matchedData(matches)
  )
  
  
  # Rename query spectrum type
  if (!is.null(query_type_col) &&
      query_type_col %in% names(df)) {
    
    names(df)[names(df) == query_type_col] <-
      "query_spectrum_type"
  }
  
  
  # Rename target spectrum type
  if (!is.null(target_type_col)) {
    
    target_col <- paste0(
      "target_",
      target_type_col
    )
    
    if (target_col %in% names(df)) {
      
      names(df)[names(df) == target_col] <-
        "target_spectrum_type"
    }
  }
  
  
  # Rename compound id
  if ("compound_id" %in% names(df)) {
    
    names(df)[names(df) == "compound_id"] <-
      "target_compound_id"
  }
  
  
  # Filter score
  df <- df[
    !is.na(df$score) &
      df$score >= threshold,
  ]
  
  
  if (nrow(df) == 0)
    return(df)
  
  
  print(
    df %>%
      count(msLevel, target_msLevel)
  )
  
  
  # Example debug
  cat(
    "Number of Herbacetin matches before filtering:",
    sum(
      grepl(
        "herbac",
        df$target_name,
        ignore.case = TRUE
      )
    ),
    "\n"
  )
  
  
  # Keep best match
  df <- df %>%
    dplyr::group_by(
      dataOrigin,
      acquisitionNum
    ) %>%
    dplyr::slice_max(
      score,
      n = 1,
      with_ties = TRUE
    ) %>%
    dplyr::slice_max(
      matched_peaks_count,
      n = 1,
      with_ties = FALSE
    ) %>%
    dplyr::ungroup()
  
  
  print(
    df %>%
      count(msLevel, target_msLevel)
  )
  
  
  return(df)
}