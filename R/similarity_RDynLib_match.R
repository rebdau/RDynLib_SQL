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
  
  
  # Polarity filtering
  
  pol_st <- Spectra::spectraData(
    st_sps,
    "polarity"
  )[["polarity"]]
  
  st_filtered <- st_sps[
    !is.na(pol_st) &
      pol_st == polarity_query
  ]
  
  
  dy_meta <- Spectra::spectraData(
    dy_sps,
    c("polarity", "name")
  )
  
  keep_dy <- !is.na(dy_meta$polarity) &
    dy_meta$polarity == polarity_target &
    !is.na(dy_meta$name) &
    dy_meta$name != "" &
    !grepl("^!", dy_meta$name)
  
  dy_filtered <- dy_sps[
    keep_dy
  ]
  
  
  if (
    length(st_filtered) == 0 ||
    length(dy_filtered) == 0
  ) {
    return(data.frame())
  }
  
  
  # Detect spectrum type column
  
  detect_spectrum_type_column <- function(sps) {
    
    cols <- Spectra::spectraVariables(sps)
    
    if ("spectrum.type" %in% cols) {
      return("spectrum.type")
    }
    
    if ("spectrum_type" %in% cols) {
      return("spectrum_type")
    }
    
    return(NULL)
  }
  
  
  # Filter spectra by spectrum type and MS level
  
  filter_spectra <- function(
    sps,
    label,
    msLevel = NULL,
    spectrum_type = NULL
  ) {
    
    # Spectrum type filtering
    
    type_col <- detect_spectrum_type_column(sps)
    
    if (!is.null(type_col)) {
      
      type_data <- Spectra::spectraData(
        sps,
        type_col
      )[[type_col]]
      
      available_types <- sort(
        unique(
          stats::na.omit(type_data)
        )
      )
      
      
      if (length(available_types) > 0) {
        
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
            
            if (is.null(type_input)) {
              return(NULL)
            }
            
            chosen_types <- trimws(
              unlist(
                strsplit(
                  type_input,
                  ","
                )
              )
            )
            
            if (
              length(chosen_types) > 0 &&
              all(
                chosen_types %in%
                available_types
              )
            ) {
              spectrum_type <- chosen_types
              break
            }
          }
          
        } else {
          
          if (
            !all(
              spectrum_type %in%
              available_types
            )
          ) {
            
            stop(
              "Invalid spectrum type for ",
              label,
              ". Available spectrum types: ",
              paste(
                available_types,
                collapse = ", "
              )
            )
          }
        }
        
        
        sps <- sps[
          !is.na(type_data) &
            type_data %in%
            spectrum_type
        ]
      }
    }
    
    
    if (length(sps) == 0) {
      return(NULL)
    }
    
    
    # MS level filtering
    
    if (
      "msLevel" %in%
      Spectra::spectraVariables(sps)
    ) {
      
      ms_data <- Spectra::spectraData(
        sps,
        "msLevel"
      )[["msLevel"]]
      
      available_ms <- sort(
        unique(
          stats::na.omit(ms_data)
        )
      )
      
      
      if (length(available_ms) > 0) {
        
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
            
            if (is.null(ms_input)) {
              return(NULL)
            }
            
            chosen_ms <- suppressWarnings(
              as.numeric(
                unlist(
                  strsplit(
                    ms_input,
                    ","
                  )
                )
              )
            )
            
            chosen_ms <- chosen_ms[
              !is.na(chosen_ms)
            ]
            
            if (
              length(chosen_ms) > 0 &&
              all(
                chosen_ms %in%
                available_ms
              )
            ) {
              msLevel <- chosen_ms
              break
            }
          }
          
        } else {
          
          if (
            !all(
              msLevel %in%
              available_ms
            )
          ) {
            
            stop(
              "Invalid MS level for ",
              label,
              ". Available MS levels: ",
              paste(
                available_ms,
                collapse = ", "
              )
            )
          }
        }
        
        
        ms_data <- Spectra::spectraData(
          sps,
          "msLevel"
        )[["msLevel"]]
        
        sps <- sps[
          !is.na(ms_data) &
            ms_data %in%
            msLevel
        ]
      }
    }
    
    
    if (length(sps) == 0) {
      return(NULL)
    }
    
    
    return(sps)
  }
  
  
  # Filter query spectra
  
  st_filtered <- filter_spectra(
    st_filtered,
    label = "query",
    msLevel = msLevel_query,
    spectrum_type = spectrum_type_query
  )
  
  if (
    is.null(st_filtered) ||
    length(st_filtered) == 0
  ) {
    return(data.frame())
  }
  
  
  # Filter target spectra
  
  dy_filtered <- filter_spectra(
    dy_filtered,
    label = "target",
    msLevel = msLevel_target,
    spectrum_type = spectrum_type_target
  )
  
  if (
    is.null(dy_filtered) ||
    length(dy_filtered) == 0
  ) {
    return(data.frame())
  }
  
  
  # Detect spectrum type columns
  
  query_type_col <- detect_spectrum_type_column(
    st_filtered
  )
  
  target_type_col <- detect_spectrum_type_column(
    dy_filtered
  )
  
  
  # Define query variables
  
  query_vars <- c(
    "msLevel",
    "rtime",
    "precursorMz",
    "acquisitionNum",
    "scanIndex",
    "dataOrigin"
  )
  
  if (!is.null(query_type_col)) {
    query_vars <- c(
      query_vars,
      query_type_col
    )
  }
  
  query_vars <- intersect(
    query_vars,
    Spectra::spectraVariables(st_filtered)
  )
  
  
  # Define target variables
  
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
  
  target_vars <- intersect(
    target_vars,
    Spectra::spectraVariables(dy_filtered)
  )
  
  
  # Extract query metadata and peaks
  
  query_data <- Spectra::spectraData(
    st_filtered,
    columns = query_vars
  )
  
  query_peaks <- Spectra::peaksData(
    st_filtered
  )
  
  
  # Extract target metadata and peaks
  
  target_data <- Spectra::spectraData(
    dy_filtered,
    columns = target_vars
  )
  
  target_peaks <- Spectra::peaksData(
    dy_filtered
  )
  
  
  # Store query peaks as mz and intensity list-columns
  
  query_data$mz <- S4Vectors::List(
    lapply(
      query_peaks,
      function(x) {
        as.numeric(x[, "mz"])
      }
    )
  )
  
  query_data$intensity <- S4Vectors::List(
    lapply(
      query_peaks,
      function(x) {
        as.numeric(x[, "intensity"])
      }
    )
  )
  
  
  # Store target peaks as mz and intensity list-columns
  
  target_data$mz <- S4Vectors::List(
    lapply(
      target_peaks,
      function(x) {
        as.numeric(x[, "mz"])
      }
    )
  )
  
  target_data$intensity <- S4Vectors::List(
    lapply(
      target_peaks,
      function(x) {
        as.numeric(x[, "intensity"])
      }
    )
  )
  
  
  # Convert query spectra to in-memory backend
  
  st_filtered <- Spectra::Spectra(
    query_data,
    backend = Spectra::MsBackendDataFrame()
  )
  
  
  # Convert target spectra to in-memory backend
  
  dy_filtered <- Spectra::Spectra(
    target_data,
    backend = Spectra::MsBackendDataFrame()
  )
  
  
  # Perform spectral matching
  
  if (
    fragments_resolution ==
    "unit.resolution"
  ) {
    
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
  
  
  # Convert matching results to data.frame
  
  df <- as.data.frame(
    MetaboAnnotation::matchedData(
      matches
    )
  )
  
  
  if (nrow(df) == 0) {
    return(df)
  }
  
  
  # Rename query spectrum type
  
  if (
    !is.null(query_type_col) &&
    query_type_col %in% names(df)
  ) {
    
    names(df)[
      names(df) == query_type_col
    ] <- "query_spectrum_type"
  }
  
  
  # Rename target spectrum type
  
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
  
  
  # Rename target identifiers if necessary
  
  if (
    "compound_id" %in% names(df) &&
    !"target_compound_id" %in% names(df)
  ) {
    
    names(df)[
      names(df) == "compound_id"
    ] <- "target_compound_id"
  }
  
  
  if (
    "compound_accession" %in% names(df) &&
    !"target_compound_accession" %in% names(df)
  ) {
    
    names(df)[
      names(df) == "compound_accession"
    ] <- "target_compound_accession"
  }
  
  
  if (
    "spectrum_id" %in% names(df) &&
    !"target_spectrum_id" %in% names(df)
  ) {
    
    names(df)[
      names(df) == "spectrum_id"
    ] <- "target_spectrum_id"
  }
  
  
  # Filter by similarity score
  
  if ("score" %in% names(df)) {
    
    df <- df[
      !is.na(df$score) &
        df$score >= threshold,
      ,
      drop = FALSE
    ]
  }
  
  # Round similarity score to 3 decimal places
  
  if ("score" %in% names(df)) {
    df$score <- round(df$score, 3)
  }
  if (nrow(df) == 0) {
    return(df)
  }
  
  
  # Show MS level combinations before best-match filtering
  
  if (
    all(
      c(
        "msLevel",
        "target_msLevel"
      ) %in% names(df)
    )
  ) {
    
    print(
      df %>%
        dplyr::count(
          msLevel,
          target_msLevel
        )
    )
  }
  
  
  # Count Herbacetin matches before best-match filtering
  
  if ("target_name" %in% names(df)) {
    
    cat(
      "Number of Herbacetin matches before filtering:",
      sum(
        grepl(
          "herbac",
          df$target_name,
          ignore.case = TRUE
        ),
        na.rm = TRUE
      ),
      "\n"
    )
  }
  
  
  # Keep the best target match for each query spectrum
  
  grouping_vars <- intersect(
    c(
      "dataOrigin",
      "acquisitionNum"
    ),
    names(df)
  )
  
  
  if (length(grouping_vars) > 0) {
    
    df <- df %>%
      dplyr::group_by(
        dplyr::across(
          dplyr::all_of(grouping_vars)
        )
      ) %>%
      dplyr::slice_max(
        score,
        n = 1,
        with_ties = TRUE
      )
    
    
    if (
      "matched_peaks_count" %in%
      names(df)
    ) {
      
      df <- df %>%
        dplyr::slice_max(
          matched_peaks_count,
          n = 1,
          with_ties = FALSE
        )
    }
    
    
    df <- df %>%
      dplyr::ungroup()
  }
  
  
  # Show final MS level combinations
  
  if (
    all(
      c(
        "msLevel",
        "target_msLevel"
      ) %in% names(df)
    )
  ) {
    
    print(
      df %>%
        dplyr::count(
          msLevel,
          target_msLevel
        )
    )
  }
  
  
  return(df)
}