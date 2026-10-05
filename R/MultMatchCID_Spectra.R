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
#' @param polarity_target 'numeric' the target polarity, 0 if negative,
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
#' @param spectrum_type_query 'character(1)' the spectrum type of the query 
#'        spectra.        
#'        
#' @param spectrum_type_target 'character(1)' the spectrum type of the target 
#'        spectra.
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
#' @param spectrum_type_query Optional character vector specifying the
#'   spectrum types to use as query spectra. Spectrum types present
#'   in a `Spectra` object can be retrieved using
#'   [available_spectrum_types()]. If NULL, no filtering by
#'   spectrum type is applied.
#'   
#' @param spectrum_type_target Optional character vector specifying the
#'   spectrum types to use as target spectra. See `spectrum_type_query`
#'   for available values. If NULL, no filtering by spectrum type
#'   is applied.
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
    fragments_resolution = c("unit.resolution", "high.resolution"),
    requirePrecursor = TRUE,
    threshold = 0.8,
    ppm = 2,
    tolerance = 0.5,
    excludeSelf = TRUE
) {
  
  fragments_resolution <- match.arg(fragments_resolution)
  
  
  # Query spectra
  
  query_meta_cols <- intersect(
    c(
      "compound_id", "polarity", "name", "spectrum_type", "spectrum.type",
      "msLevel", "precursorMz", "rtime", "expid", "retention_time"
    ),
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
  
  spectrum_col <- intersect(
    c("spectrum_type", "spectrum.type"),
    query_meta_cols
  )
  
  spectrum_col <- if (length(spectrum_col) > 0) {
    spectrum_col[1]
  } else {
    NULL
  }
  
  
  vars_keep <- c(
    "compound_id", "polarity", "msLevel", "name",
    "rtime", "precursorMz", "expid", "retention_time"
  )
  
  if (!is.null(spectrum_col))
    vars_keep <- c(spectrum_col, vars_keep)
  
  query_sps <- selectSpectraVariables(
    query_sps,
    vars_keep
  )
  
  
  # Spectrum type filter applied to all spectra of the compound,
  # regardless of MS level
  
  if (!is.null(spectrum_col)) {
    
    available_types <- unique(
      spectraData(query_sps, spectrum_col)[[1]]
    )
    
    if (is.null(spectrum_type_query) &&
        length(available_types) > 1) {
      
      spectrum_type_query <- rstudioapi::showPrompt(
        "Spectrum type selection (query)",
        paste(
          "Available types:",
          paste(available_types, collapse = ", ")
        ),
        default = available_types[1]
      )
    }
    
    if (!is.null(spectrum_type_query) &&
        spectrum_type_query %in% available_types) {
      
      query_sps <- query_sps[
        spectraData(query_sps, spectrum_col)[[1]] ==
          spectrum_type_query
      ]
    }
  }
  
  
  if (length(query_sps) == 0)
    stop("No query spectra left after spectrum_type filtering.")
  
  
  # Apply minimum number of fragment ions filter to query spectra
  
  query_peaks_list <- peaksData(query_sps)
  
  keep_idx <- vapply(
    query_peaks_list,
    nrow,
    integer(1)
  ) >= minIons
  
  query_sps <- query_sps[keep_idx]
  
  if (length(query_sps) == 0)
    stop("No query spectra left after minIons filtering.")
  
  
  # Select query MS level
  
  available_ms_query <- sort(
    unique(msLevel(query_sps))
  )
  
  if (is.null(ms_level_query)) {
    
    if (length(available_ms_query) > 1) {
      
      sel <- rstudioapi::showPrompt(
        "MS level selection (query)",
        paste0(
          "Available MS levels for this compound: ",
          paste(available_ms_query, collapse = ", "),
          "\nEnter one value, or multiple values separated by commas, e.g. 3,4"
        ),
        default = as.character(available_ms_query[1])
      )
      
      ms_level_query <- as.numeric(
        trimws(strsplit(sel, ",")[[1]])
      )
      
    } else {
      
      ms_level_query <- available_ms_query
    }
  }
  
  
  if (!all(ms_level_query %in% available_ms_query)) {
    
    stop(sprintf(
      "Invalid query MS level(s). Available MS levels for this compound: %s",
      paste(available_ms_query, collapse = ", ")
    ))
  }
  
  
  query_sps <- query_sps[
    msLevel(query_sps) %in% ms_level_query
  ]
  
  if (length(query_sps) == 0)
    stop("No query spectra left after ms_level_query filtering.")
  
  
  n_query <- length(query_sps)
  query_ms_levels <- msLevel(query_sps)
  
  message(sprintf(
    "Compound %s: %d query spectrum/spectra retained for MS level(s): %s",
    compound_id,
    n_query,
    paste(
      sort(unique(query_ms_levels)),
      collapse = ", "
    )
  ))
  
  
  # Target spectra
  
  target_meta_cols <- intersect(
    c(
      "compound_id", "polarity", "name", "spectrum_type", "spectrum.type",
      "msLevel", "precursorMz", "rtime", "expid", "retention_time"
    ),
    spectraVariables(spectra_db)
  )
  
  meta_target <- spectraData(
    spectra_db,
    target_meta_cols
  )
  
  
  # Keep target spectra with the requested polarity and a valid name
  
  target_filter <-
    meta_target$polarity == polarity_target &
    !is.na(meta_target$name) &
    meta_target$name != "" &
    !grepl("^!", meta_target$name)
  
  
  # Exclude spectra from the query compound itself to avoid
  # trivial self-matches with a score of 1
  
  if (excludeSelf &&
      "compound_id" %in% colnames(meta_target)) {
    
    target_filter <- target_filter &
      (
        is.na(meta_target$compound_id) |
          meta_target$compound_id != compound_id
      )
  }
  
  
  target_sps_all <- spectra_db[target_filter]
  
  
  vars_keep_t <- c(
    "compound_id", "polarity", "msLevel", "rtime",
    "name", "precursorMz", "expid", "retention_time"
  )
  
  if (!is.null(spectrum_col))
    vars_keep_t <- c(spectrum_col, vars_keep_t)
  
  target_sps_all <- selectSpectraVariables(
    target_sps_all,
    vars_keep_t
  )
  
  
  # Apply spectrum type filter to target spectra
  
  if (!is.null(spectrum_col)) {
    
    available_types_target <- unique(
      spectraData(target_sps_all, spectrum_col)[[1]]
    )
    
    display_types <- ifelse(
      is.na(available_types_target),
      "NA",
      available_types_target
    )
    
    if (is.null(spectrum_type_target)) {
      
      cat("\nAvailable target spectrum types:\n")
      print(display_types)
      
      spectrum_type_target <- rstudioapi::showPrompt(
        "Spectrum type selection (target)",
        "Enter one of the available spectrum types:",
        default = display_types[1]
      )
    }
    
    if (identical(spectrum_type_target, "NA"))
      spectrum_type_target <- NA
    
    target_sps_all <- target_sps_all[
      spectraData(target_sps_all, spectrum_col)[[1]] ==
        spectrum_type_target |
        (
          is.na(
            spectraData(target_sps_all, spectrum_col)[[1]]
          ) &
            is.na(spectrum_type_target)
        )
    ]
  }
  
  
  # Apply minimum number of fragment ions filter to target spectra
  
  target_sps_all <- target_sps_all[
    vapply(
      peaksData(target_sps_all),
      nrow,
      integer(1)
    ) >= minIons
  ]
  
  if (length(target_sps_all) == 0)
    return("No target spectra left after filtering.")
  
  
  # Select target MS level
  
  available_ms_target <- sort(
    unique(msLevel(target_sps_all))
  )
  
  if (is.null(ms_level_target)) {
    
    if (
      length(available_ms_target) > 1 &&
      !all(ms_level_query %in% available_ms_target)
    ) {
      
      sel <- rstudioapi::showPrompt(
        "MS level selection (target)",
        paste0(
          "Available target MS levels: ",
          paste(available_ms_target, collapse = ", "),
          "\nEnter one value, or multiple values separated by commas."
        ),
        default = as.character(available_ms_target[1])
      )
      
      ms_level_target <- as.numeric(
        trimws(strsplit(sel, ",")[[1]])
      )
      
    } else {
      
      ms_level_target <- ms_level_query
    }
  }
  
  
  if (!all(ms_level_target %in% available_ms_target)) {
    
    stop(sprintf(
      "Invalid target MS level(s). Available MS levels: %s",
      paste(available_ms_target, collapse = ", ")
    ))
  }
  
  
  target_sps_all <- target_sps_all[
    msLevel(target_sps_all) %in% ms_level_target
  ]
  
  if (length(target_sps_all) == 0)
    return("No target spectra left after ms_level_target filtering.")
  
  
  target_ms_levels_all <- msLevel(target_sps_all)
  
  ms_target_fixed <- !setequal(
    ms_level_target,
    ms_level_query
  )
  
  
  # Spectral matching
  
  all_results <- list()
  
  for (i in seq_len(n_query)) {
    
    q_sp <- query_sps[i]
    q_ms <- query_ms_levels[i]
    
    
    # If target MS levels were explicitly set to different levels,
    # use all selected target spectra. Otherwise, compare spectra
    # belonging to the same MS level.
    
    if (ms_target_fixed) {
      
      target_sps <- target_sps_all
      
    } else {
      
      target_sps <- target_sps_all[
        target_ms_levels_all == q_ms
      ]
    }
    
    
    if (length(target_sps) == 0) {
      
      message(sprintf(
        "  - MS%d: no target spectra available; spectrum skipped.",
        q_ms
      ))
      
      next
    }
    
    
    # Unit-resolution spectral matching
    
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
        query = q_sp,
        target = target_sps,
        param = param
      )
      
      
      # High-resolution spectral matching
      
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
        query = q_sp,
        target = target_sps,
        param = param
      )
    }
    
    
    # Extract and filter matching results
    
    df_i <- as.data.frame(
      MetaboAnnotation::matchedData(matches)
    )
    
    if (nrow(df_i) == 0)
      next
    
    
    # Retain matches above the similarity threshold
    
    df_i <- df_i[
      !is.na(df_i$score) &
        df_i$score >= threshold,
      ,
      drop = FALSE
    ]
    
    if (nrow(df_i) == 0)
      next
    
    
    # Retain matches with the required minimum number of matched ions
    
    df_i <- df_i[
      df_i$matched_peaks_count >= minIons,
      ,
      drop = FALSE
    ]
    
    if (nrow(df_i) == 0)
      next
    
    
    # Add query and matching information
    
    df_i$query_spectrum_index <- i
    df_i$query_msLevel_used <- q_ms
    df_i$`#FragIon` <- nrow(
      peaksData(q_sp)[[1]]
    )
    df_i$ComIons <- df_i$matched_peaks_count
    df_i$DotIons <- df_i$score
    
    all_results[[length(all_results) + 1]] <- df_i
  }
  
  
  if (length(all_results) == 0)
    return("No matches found for any MS level of this compound.")
  
  
  # Combine matching results
  
  df <- do.call(
    rbind,
    all_results
  )
  
  
  # Rename and add output variables
  
  if ("target_compound_id" %in% colnames(df))
    df$COMPID <- df$target_compound_id
  
  if ("precursorMz" %in% colnames(df))
    df$query_precursorMz <- df$precursorMz
  
  if ("retention_time" %in% colnames(df))
    df$query_retention_time <- df$retention_time
  
  if ("expid" %in% colnames(df))
    df$query_expid <- df$expid
  
  if ("msLevel" %in% colnames(df))
    df$query_msLevel <- df$msLevel
  
  
  # Define output column order
  
  cols_order <- c(
    "query_spectrum_index",
    "target_compound_id",
    "query_msLevel_used",
    "target_msLevel",
    "query_precursorMz",
    "target_precursorMz",
    "#FragIon",
    "ComIons",
    "DotIons",
    "target_name",
    "target_retention_time",
    "query_retention_time",
    "target_expid",
    "query_expid"
  )
  
  cols_order <- cols_order[
    cols_order %in% colnames(df)
  ]
  
  df <- df[
    ,
    cols_order,
    drop = FALSE
  ]
  
  
  # Round similarity scores
  
  if ("DotIons" %in% colnames(df))
    df$DotIons <- round(df$DotIons, 2)
  
  
  # Round precursor m/z values for unit-resolution matching
  
  if (fragments_resolution == "unit.resolution") {
    
    if ("query_precursorMz" %in% colnames(df))
      df$query_precursorMz <- round(
        df$query_precursorMz,
        2
      )
    
    if ("target_precursorMz" %in% colnames(df))
      df$target_precursorMz <- round(
        df$target_precursorMz,
        2
      )
  }
  
  
  rownames(df) <- NULL
  
  return(df)
}