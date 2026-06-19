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

#' Enrich FTMS trees with MS3/MS4 spectra
#'
#' @description
#' Rebuilds Spectra object with MS3/MS4 enrichment.
#'
#' @param ftms_best Spectra object
#' @param sd_filtered data.frame filtered spectra
#'
#' @return Spectra object enriched
#' @author Ahlam Mentag, Rebecca Dauwe
enrich_ms3_ms4 <- function(ftms_best, sd_filtered) {
  
  library(dplyr)
  
  sd_best <- as.data.frame(spectraData(ftms_best))
  
  # keys
  sd_best$.key <- paste(sd_best$dataOrigin, sd_best$scanIndex, sep = "_")
  sd_filtered$.key <- paste(sd_filtered$dataOrigin, sd_filtered$scanIndex, sep = "_")
  
  #helper column
  sd_best$precMzRound <- round(sd_best$precursorMz)
  sd_filtered$precMzRound <- round(sd_filtered$precursorMz)
  
  # MS4 map 
  ms4_map <- sd_filtered |>
    dplyr::filter(msLevel == 4)
  
  # detect column names
  if (!"precScanNum" %in% names(ms4_map)) {
    if ("precursorScanNum" %in% names(ms4_map)) {
      ms4_map$precScanNum <- ms4_map$precursorScanNum
    } else {
      ms4_map$precScanNum <- NA
    }
  }
  
  if (!"peaksCount" %in% names(ms4_map)) {
    ms4_map$peaksCount <- NA
  }
  
  ms4_map <- ms4_map |>
    dplyr::select(MSntreeID, precScanNum, acquisitionNum, peaksCount) |>
    dplyr::rename(
      ms3_acq = precScanNum,
      ms4_acq = acquisitionNum
    )
  #enrichment 
  additional_ms3 <- add_ms3_ms4_candidates(
    sd_best,
    sd_filtered,
    ms4_map
  )
  
  # if nothing
  if (is.null(additional_ms3) || nrow(additional_ms3) == 0) {
    message("No enrichment found")
    
    spectraData(ftms_best)$MergedMSntreeID <- as.integer(
      factor(spectraData(ftms_best)$feature_id)
    )
    
    return(ftms_best)
  }
  
 
  # keep only real spectral rows
  valid_keys <- sd_best$.key
  
  all_df <- bind_rows(sd_best, additional_ms3)
  all_df <- all_df[all_df$.key %in% valid_keys, ]
  
  # reorder
  all_df <- all_df[match(valid_keys, all_df$.key), ]
  
  # rebuild
  spectraData(ftms_best) <- S4Vectors::DataFrame(all_df)
  
  # final feature ID 
  spectraData(ftms_best)$MergedMSntreeID <- as.integer(
    factor(spectraData(ftms_best)$feature_id)
  )
  
  ftms_best
}



