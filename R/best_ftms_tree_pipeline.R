#' Filter MS trees by MS2 precursor intensity
#'
#' @description
#' Keeps MS trees whose MS2 precursor intensities are above a quantile threshold.
#'
#' @param sd data.frame containing MS spectra metadata
#' @param quantile_cut numeric(1) quantile cutoff (default 0.25)
#'
#' @return data.frame filtered spectra
#' @author Ahlam Mentag, Rebecca Dauwe
filter_intense_prec <- function(sd, quantile_cut = 0.25) {
  
  ms2 <- sd[sd$msLevel == 2, ]
  threshold <- quantile(ms2$precursorIntensity, quantile_cut, na.rm = TRUE)
  
  strong_ids <- unique(ms2$MSntreeID[ms2$precursorIntensity > threshold])
  
  sd_filtered <- sd[sd$MSntreeID %in% strong_ids, ]
  
  message(length(unique(sd$feature_id)) - length(unique(sd_filtered$feature_id)),
          " features removed by intensity filter")
  
  sd_filtered
}

#' Deduplicate MS3 spectra
#'
#' @description
#' Removes redundant MS3 spectra by keeping the most intense one per group.
#'
#' @param sd data.frame containing MS spectra metadata
#'
#' @return data.frame deduplicated spectra
#' @author Ahlam Mentag, Rebecca Dauwe
deduplicate_ms3 <- function(sd) {
  
  sd$precursorMz_round <- round(sd$precursorMz)
  
  ms3 <- sd[sd$msLevel == 3, ]
  
  ms3_dedup <- ms3 |>
    dplyr::group_by(dataOrigin, precursorMz_round, precursorIntensity) |>
    dplyr::slice_max(peaksCount, n = 1, with_ties = FALSE) |>
    dplyr::ungroup()
  
  sd_new <- rbind(sd[sd$msLevel != 3, ], ms3_dedup)
  
  message("MS3 dedup complete")
  
  sd_new
}

#' Select longest MS tree per feature
#'
#' @description
#' Selects the longest fragmentation tree for each feature.
#'
#' @param sd data.frame containing MS spectra metadata
#'
#' @return data.frame selected trees
#' @author Ahlam Mentag, Rebecca Dauwe
select_longest_trees <- function(sd) {
  
  tree_length <- aggregate(msLevel ~ MSntreeID + feature_id, sd, max)
  max_len <- aggregate(msLevel ~ feature_id, tree_length, max)
  
  longest <- merge(tree_length, max_len,
                   by = c("feature_id", "msLevel"))
  
  longest
}

#' Select best MS tree per feature
#'
#' @description
#' Selects best MS tree based on peak count and precursor intensity.
#'
#' @param sd data.frame containing MS spectra metadata
#' @param longest_trees data.frame output of select_longest_trees()
#'
#' @return data.frame selected trees
#' @author Ahlam Mentag, Rebecca Dauwe
select_best_tree <- function(sd, longest_trees) {
  
  tree_peaks <- aggregate(peaksCount ~ MSntreeID, sd, sum)
  
  candidates <- merge(longest_trees, tree_peaks, by = "MSntreeID")
  
  max_peaks <- aggregate(peaksCount ~ feature_id, candidates, max)
  
  candidates <- merge(candidates, max_peaks,
                      by = c("feature_id", "peaksCount"))
  
  prec_max <- aggregate(precursorIntensity ~ MSntreeID, sd, max)
  
  candidates <- merge(candidates, prec_max, by = "MSntreeID")
  
  candidates <- candidates[order(candidates$feature_id,
                                 -candidates$precursorIntensity), ]
  
  best <- candidates[!duplicated(candidates$feature_id),
                     c("feature_id", "MSntreeID")]
  
  best
}

#' Run full MS tree selection pipeline
#'
#' @description
#' Executes filtering, deduplication and best tree selection.
#'
#' @param ftms_msn_tree Spectra object
#'
#' @return list with filtered data and selected FTMS object
#' @author Ahlam Mentag, Rebecca Dauwe
best_tree_pipeline <- function(ftms_msn_tree) {
  
  sd <- as.data.frame(spectraData(ftms_msn_tree))
  
  sd1 <- filter_intense_prec(sd)
  sd2 <- deduplicate_ms3(sd1)
  
  longest <- select_longest_trees(sd2)
  best <- select_best_tree(sd2, longest)
  
  best_MSntreeIDs <- best$MSntreeID
  
  ftms_best <- ftms_msn_tree[
    spectraData(ftms_msn_tree)$MSntreeID %in% best_MSntreeIDs
  ]
  
  list(
    sd_filtered = sd2,
    best_trees = best,
    ftms_best = ftms_best
  )
}


#' Select best MS3 and optional MS4 spectra
#'
#' @description
#' Selects best MS3 spectrum and attaches MS4 if available.
#'
#' @param ms3_df data.frame MS3 candidates
#' @param sd_filtered data.frame full filtered spectra
#'
#' @return data.frame enriched spectra
#' @author Ahlam Mentag, Rebecca Dauwe
select_best_ms3ms4 <- function(ms3_df, sd_filtered) {
  
  ms3_df$peaksCount <- ifelse(is.null(ms3_df$peaksCount), NA, ms3_df$peaksCount)
  ms3_df$precursorIntensity <- ifelse(is.null(ms3_df$precursorIntensity), NA, ms3_df$precursorIntensity)
  
  best_ms3 <- ms3_df |>
    dplyr::arrange(desc(peaksCount), desc(precursorIntensity)) |>
    dplyr::slice(1)
  
  # no MS4 link
  if (is.na(best_ms3$ms4_acq)) {
    return(best_ms3)
  }
  
  ms4_row <- sd_filtered |>
    dplyr::filter(
      msLevel == 4,
      acquisitionNum %in% best_ms3$ms4_acq,
      feature_id == best_ms3$feature_id
    )
  

  
  all_cols <- union(names(best_ms3), names(ms4_row))
  
  # add missing columns with NA
  for (col in all_cols) {
    if (!col %in% names(best_ms3)) best_ms3[[col]] <- NA
    if (!col %in% names(ms4_row))  ms4_row[[col]]  <- NA
  }
  
  best_ms3 <- best_ms3[, all_cols, drop = FALSE]
  ms4_row  <- ms4_row[, all_cols, drop = FALSE]
  
  dplyr::bind_rows(best_ms3, ms4_row)
}



#' Generate MS3/MS4 candidates per feature
#'
#' @description
#' Expands MS3 candidates and links MS4 spectra when available.
#'
#' @param sd_best data.frame best MS trees
#' @param sd_filtered data.frame filtered spectra
#' @param ms4_map data.frame MS4 mapping table
#'
#' @return data.frame candidate spectra
#' @author Ahlam Mentag, Rebecca Dauwe
add_ms3_ms4_candidates <- function(sd_best, sd_filtered, ms4_map) {
  
  library(dplyr)
  
  sd_best$precMzRound <- round(sd_best$precursorMz)
  sd_filtered$precMzRound <- round(sd_filtered$precursorMz)
  
  out_list <- list()
  
  for (f in unique(sd_best$feature_id)) {
    
    best_tree <- sd_best[sd_best$feature_id == f, ]
    best_prec <- unique(best_tree$precMzRound[best_tree$msLevel == 3])
    
    ms3_candidates <- sd_filtered |>
      filter(
        feature_id == f,
        msLevel == 3,
        !(precMzRound %in% best_prec)
      )
    
    if (nrow(ms3_candidates) == 0) next
    
    ms3_candidates <- left_join(
      ms3_candidates,
      ms4_map,
      by = c("MSntreeID", "acquisitionNum" = "ms3_acq")
    )
    
    ms3_candidates$precMzRound <- round(ms3_candidates$precursorMz)
    
    grouped <- split(ms3_candidates, ms3_candidates$precMzRound)
    
    out_list[[f]] <- dplyr::bind_rows(
      lapply(grouped, select_best_ms3ms4, sd_filtered)
    )
  }
  
  dplyr::bind_rows(out_list)
}

#' Enrich FTMS trees with additional MS3/MS4 spectra
#'
#' @description
#' Enriches the selected representative MSn trees with complementary
#' MS3 and MS4 spectra from other fragmentation trees associated with
#' the same chromatographic feature.
#'
#' @param ftms_best Spectra object containing the selected best trees.
#' @param sd_filtered data.frame containing filtered spectra metadata.
#' @param ftms_msn_tree Original Spectra object containing all spectra.
#'
#' @return Spectra object containing the enriched fragmentation trees.
#'
#' @author Ahlam Mentag, Rebecca Dauwe

enrich_ms3_ms4 <- function(
    ftms_best,
    sd_filtered,
    ftms_msn_tree
) {
  
  # Extract metadata
  
  sd_best <- as.data.frame(
    Spectra::spectraData(ftms_best)
  )
  
  sd_original <- as.data.frame(
    Spectra::spectraData(ftms_msn_tree)
  )
  
  
  # Create unique spectrum keys
  
  make_key <- function(x) {
    
    paste(
      x$dataOrigin,
      x$scanIndex,
      sep = "_"
    )
    
  }
  
  sd_best$.key <- make_key(sd_best)
  
  sd_filtered$.key <- make_key(sd_filtered)
  
  sd_original$.key <- make_key(sd_original)
  
  
  # Helper columns
  
  sd_best$precMzRound <- round(
    sd_best$precursorMz
  )
  
  sd_filtered$precMzRound <- round(
    sd_filtered$precursorMz
  )
  
  
  # Build MS4 mapping table
  
  ms4_map <- sd_filtered |>
    dplyr::filter(msLevel == 4)
  
  
  # Check precursor scan column
  
  if (!"precScanNum" %in% names(ms4_map)) {
    
    if ("precursorScanNum" %in% names(ms4_map)) {
      
      ms4_map$precScanNum <-
        ms4_map$precursorScanNum
      
    } else {
      
      stop(
        "Missing MS4 precursor scan information."
      )
      
    }
    
  }
  
  
  # Create MS3-MS4 mapping
  
  ms4_map <- ms4_map |>
    dplyr::select(
      MSntreeID,
      precScanNum,
      acquisitionNum
    ) |>
    dplyr::rename(
      ms3_acq = precScanNum,
      ms4_acq = acquisitionNum
    )
  
  
  # Identify additional MS3/MS4 spectra
  
  additional <- add_ms3_ms4_candidates(
    sd_best = sd_best,
    sd_filtered = sd_filtered,
    ms4_map = ms4_map
  )
  
  
  # Helper to assign merged tree IDs
  assign_merged_ids <- function(sps) {
    
    # Convert to memory backend
    sps <- Spectra::setBackend(
      sps,
      Spectra::MsBackendMemory()
    )
    
    # Extract feature IDs
    feature_ids <- sps$feature_id
    
    # Assign one merged tree ID per feature
    sps$MergedMSntreeID <- as.integer(
      factor(feature_ids)
    )
    
    sps
  }
  
  
  # No enrichment found
  
  if (is.null(additional) ||
      nrow(additional) == 0) {
    
    message(
      "No additional MS3/MS4 spectra found."
    )
    
    return(
      assign_merged_ids(ftms_best)
    )
    
  }
  
  
  # Remove spectra already present in best trees
  
  additional <- additional |>
    dplyr::filter(
      !.key %in% sd_best$.key
    ) |>
    dplyr::distinct(
      .key,
      .keep_all = TRUE
    )
  
  
  # Check whether additional spectra remain
  
  if (nrow(additional) == 0) {
    
    message(
      "All candidate spectra are already present ",
      "in the selected trees."
    )
    
    return(
      assign_merged_ids(ftms_best)
    )
    
  }
  
  
  # Match additional spectra to original Spectra object
  
  idx <- match(
    additional$.key,
    sd_original$.key
  )
  
  
  # Remove unmatched spectra
  
  valid <- !is.na(idx)
  
  if (any(!valid)) {
    
    warning(
      sum(!valid),
      " additional spectra could not be found ",
      "in the original Spectra object."
    )
    
  }
  
  idx <- idx[valid]
  
  
  # Check whether spectra were retrieved
  
  if (length(idx) == 0) {
    
    message(
      "No additional spectra could be retrieved."
    )
    
    return(
      assign_merged_ids(ftms_best)
    )
    
  }
  
  
  # Retrieve actual spectra, including fragment peaks
  
  additional_spectra <- ftms_msn_tree[idx]
  
  
  # Combine original and additional spectra
  
  ftms_enriched <- c(
    ftms_best,
    additional_spectra
  )
  
  
  # Convert to memory backend before updating metadata
  
  ftms_enriched <- Spectra::setBackend(
    ftms_enriched,
    Spectra::MsBackendMemory()
  )
  
  
  # Assign one merged tree ID per feature
  ftms_enriched <- assign_merged_ids(
    ftms_enriched
  )
  
  
  # Summary
  
  message(
    "Original spectra: ",
    length(ftms_best)
  )
  
  message(
    "Additional spectra: ",
    length(additional_spectra)
  )
  
  message(
    "Final spectra: ",
    length(ftms_enriched)
  )
  
  message(
    "Enriched features: ",
    length(
      unique(
        Spectra::spectraData(
          additional_spectra,
          "feature_id"
        )[["feature_id"]]
      )
    )
  )
  
  
  # Return enriched Spectra object
  
  ftms_enriched
  
}