# In case of duplicated assignments: Resolving multiple assignments of MS2 spectra to features

# We evaluate the cases above to reduce/remove multiple
# assignments and keep only a single MS2-feature association. First, we
# select the duplicated MS2 spectra by identifying those with the same
# *dataOrigin* and *scanIndex* but different *feature_id* to detect
# duplicates assigned to multiple features. From these duplicated MS2
# spectra, we then check if the difference between the precursor m/z of the
# MS2 and the median m/z of the feature (*feature_mzmed*) is within the
# defined threshold. Only those within the threshold are considered valid
# matches. If no valid match is found, the first mapping is kept and all
# feature IDs are stored as ambiguous. If multiple valid matches are found,
# the feature whose retention time (*feature_rtmed*) is closest to the MS2
# retention time is selected, and the remaining ones are marked as ambiguous.
# All unique (non-duplicated) MS2-feature mappings are retained directly.
# The result is a filtered MS2 object containing only unique assignments and
# an *ambiguous* field listing alternative feature IDs for each retained MS2.

#' @param ms2 A `Spectra` object containing MS2 spectra and their associated
#' feature metadata (including `feature_id`, `feature_mzmed`, and `feature_rtmed`).
#' @param threshold Numeric(1). Maximum allowed difference between precursor m/z
#' of the MS2 spectrum and the feature median m/z (*feature_mzmed*) to be
#' considered a valid match (default = 0.002).


filter_unique_MS2 <- function(ms2, threshold = 0.002) {
  
  sp <- spectraData(ms2)
  #'group by dataOrigin and scanIndex to identify
  #'duplicated assignments in each dataOrigin
  sp$ms2_id <- paste0(sp$dataOrigin, "_", sp$scanIndex)
  sp$feature_id <- as.character(sp$feature_id)
  sp$ambiguous <- vector("list", nrow(sp))
  
  # Find duplicated MS2
  ms2_table <- table(sp$ms2_id)
  duplicated_ms2 <- names(ms2_table)[ms2_table > 1]
  # keep unique MS2 assignments directly
  keep_rows <- which(!sp$ms2_id %in% duplicated_ms2)  
  
  if (length(duplicated_ms2) > 0) {
    # Get indices of duplicated MS2
    dup_idx <- which(sp$ms2_id %in% duplicated_ms2)
    
    # Split by ms2_id
    dup_groups <- split(dup_idx, sp$ms2_id[dup_idx])
    
    for (ms2_name in names(dup_groups)) {
      idx <- dup_groups[[ms2_name]]
      feature_ids <- sp$feature_id[idx]
      
      ms2_mz <- precursorMz(ms2)[idx[1]]
      ms2_rt <- rtime(ms2)[idx[1]]
      
      mz_features <- sp$feature_mzmed[idx]
      rt_features <- sp$feature_rtmed[idx]
      
      # m/z filter
      valid <- feature_ids[abs(ms2_mz - mz_features) <= threshold]
      
      if (length(valid) == 0) {
        # Keep first row, all features in ambiguous
        keep_row <- idx[1]
        sp$ambiguous[keep_row] <- list(feature_ids)
        keep_rows <- c(keep_rows, keep_row)
        next
      }
      
      # pick closest by RT if multiple
      if (length(valid) > 1) {
        rt_diff <- abs(ms2_rt - rt_features[match(valid, feature_ids)])
        keep_feat <- valid[which.min(rt_diff)]
      } else {
        keep_feat <- valid
      }
      
      # Keep only row corresponding to selected feature
      keep_row <- idx[which(feature_ids == keep_feat)[1]]
      keep_rows <- c(keep_rows, keep_row)
      
      # Store all other original features as ambiguous
      sp$ambiguous[keep_row] <- list(setdiff(feature_ids, keep_feat))
    }
  }
  
  ms2_filtered <- ms2[keep_rows, ]
  sp_final <- spectraData(ms2_filtered)
  sp_final$ambiguous <- sp$ambiguous[keep_rows]
  spectraData(ms2_filtered) <- sp_final
  
  ms2_filtered
}


# Here we call the function
#qtof_filtered <- filter_unique_MS2(ms2, threshold = 0.002)