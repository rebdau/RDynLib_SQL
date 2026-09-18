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
    msLevel_query = NULL,
    msLevel_target = NULL,
    spectrum_type_query = NULL,
    spectrum_type_target = NULL,
    fragments_resolution = c("unit.resolution", "high.resolution"),
    requirePrecursor = TRUE,
    threshold = 0.8,
    ppm = 2,
    tolerance = 0.5
) {
  
  fragments_resolution <- match.arg(fragments_resolution)
  
  
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
  
  

  
  detect_spectrum_type_column <- function(sps) {
    
    cols <- spectraVariables(sps)
    
    if ("spectrum.type" %in% cols)
      return("spectrum.type")
    
    if ("spectrum_type" %in% cols)
      return("spectrum_type")
    
    return(NULL)
  }
  
  
  filter_spectra <- function(
    sps,
    label,
    msLevel = NULL,
    spectrum_type = NULL
  ) {
    
    # Spectrum type
    
    type_col <- detect_spectrum_type_column(sps)
    
    if (!is.null(type_col)) {
      
      type_data <- spectraData(
        sps,
        type_col
      )[[type_col]]
      
      available_types <- sort(
        unique(na.omit(type_data))
      )
      
      
      if (length(available_types) > 0) {
        
        # If spectrum type was not supplied:
        # ask interactively
        if (is.null(spectrum_type)) {
          
          repeat {
            
            type_input <- rstudioapi::showPrompt(
              paste("Spectrum type -", label),
              paste(
                "Available types:",
                paste(
                  available_types,
                  collapse = ", "
                )
              )
            )
            
            # User cancelled
            if (is.null(type_input))
              return(NULL)
            
            chosen_types <- trimws(
              unlist(
                strsplit(type_input, ",")
              )
            )
            
            if (
              length(chosen_types) > 0 &&
              all(chosen_types %in% available_types)
            ) {
              spectrum_type <- chosen_types
              break
            }
          }
          
        } else {
          
          # Check supplied spectrum types
          if (!all(spectrum_type %in% available_types)) {
            
            stop(
              "Invalid spectrum type for ", label,
              ". Available spectrum types: ",
              paste(
                available_types,
                collapse = ", "
              )
            )
          }
        }
        
        
        # Apply spectrum type filter
        sps <- sps[
          type_data %in% spectrum_type
        ]
      }
    }
    
    
    if (length(sps) == 0)
      return(NULL)
    

    
    if ("msLevel" %in% spectraVariables(sps)) {
      
      ms_data <- spectraData(
        sps,
        "msLevel"
      )[["msLevel"]]
      
      available_ms <- sort(
        unique(na.omit(ms_data))
      )
      
      
      if (length(available_ms) > 0) {
        
        # If MS level was not supplied:
        # ask interactively
        if (is.null(msLevel)) {
          
          repeat {
            
            ms_input <- rstudioapi::showPrompt(
              paste("MS level -", label),
              paste(
                "Available MS levels:",
                paste(
                  available_ms,
                  collapse = ", "
                )
              )
            )
            
            # User cancelled
            if (is.null(ms_input))
              return(NULL)
            
            chosen_ms <- suppressWarnings(
              as.numeric(
                unlist(
                  strsplit(ms_input, ",")
                )
              )
            )
            
            chosen_ms <- chosen_ms[
              !is.na(chosen_ms)
            ]
            
            if (
              length(chosen_ms) > 0 &&
              all(chosen_ms %in% available_ms)
            ) {
              msLevel <- chosen_ms
              break
            }
          }
          
        } else {
          
          # Check supplied MS levels
          if (!all(msLevel %in% available_ms)) {
            
            stop(
              "Invalid MS level for ", label,
              ". Available MS levels: ",
              paste(
                available_ms,
                collapse = ", "
              )
            )
          }
        }
        
        
        # Apply MS level filter
        sps <- sps[
          ms_data %in% msLevel
        ]
      }
    }
    
    
    if (length(sps) == 0)
      return(NULL)
    
    
    return(sps)
  }
  
  

  # Query filtering
  
  st_filtered <- filter_spectra(
    st_filtered,
    label = "query",
    msLevel = msLevel_query,
    spectrum_type = spectrum_type_query
  )
  
  if (is.null(st_filtered) || length(st_filtered) == 0)
    return(data.frame())
  
  

  # Target filtering
  dy_filtered <- filter_spectra(
    dy_filtered,
    label = "target",
    msLevel = msLevel_target,
    spectrum_type = spectrum_type_target
  )
  
  if (is.null(dy_filtered) || length(dy_filtered) == 0)
    return(data.frame())
  
  
  query_type_col <- detect_spectrum_type_column(
    st_filtered
  )
  
  target_type_col <- detect_spectrum_type_column(
    dy_filtered
  )
  
  
  
  target_vars <- c(
    "msLevel",
    "rtime",
    "precursorMz",
    "name",
    "scanIndex",
    "dataOrigin",
    "compound_id",
    "compound_accession",
    "spectrum_id"
  )
  
  
  if (!is.null(target_type_col)) {
    
    target_vars <- c(
      target_vars,
      target_type_col
    )
  }
  
  
  dy_filtered <- selectSpectraVariables(
    dy_filtered,
    target_vars
  )

  
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
    
    param <- MetaboAnnotation::CompareSpectraParam(
      ppm = ppm,
      tolerance = tolerance,
      MAPFUN = joinPeaksGnps,
      FUN = MsCoreUtils::gnps,
      threshold = threshold,
      requirePrecursor = requirePrecursor,
      matchedPeaksCount = TRUE
    )
    
    
    matches <- MetaboAnnotation::matchSpectra(
      query = st_filtered,
      target = dy_filtered,
      param = param
    )
  }
  
  
  df <- as.data.frame(
    MetaboAnnotation::matchedData(matches)
  )
  
  
  if (
    !is.null(query_type_col) &&
    query_type_col %in% names(df)
  ) {
    
    names(df)[
      names(df) == query_type_col
    ] <- "query_spectrum_type"
  }
  
  
  if (!is.null(target_type_col)) {
    
    target_col <- paste0(
      "target_",
      target_type_col
    )
    
    if (target_col %in% names(df)) {
      
      names(df)[
        names(df) == target_col
      ] <- "target_spectrum_type"
    }
  }
  
  
  
  if ("compound_id" %in% names(df)) {
    
    names(df)[
      names(df) == "compound_id"
    ] <- "target_compound_id"
  }
  
  if ("compound_accession" %in% names(df)) {
    
    names(df)[
      names(df) == "compound_accession"
    ] <- "target_compound_accession"
  }
  if ("spectrum_id" %in% names(df)) {
    
    names(df)[
      names(df) == "spectrum_id"
    ] <- "target_spectrum_id"
  }
  
  df <- df[
    !is.na(df$score) &
      df$score >= threshold,
  ]
  
  
  if (nrow(df) == 0)
    return(df)
  
  print(
    df %>%
      dplyr::count(
        msLevel,
        target_msLevel
      )
  )
  
  
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
      dplyr::count(
        msLevel,
        target_msLevel
      )
  )
  
  
  return(df)
}