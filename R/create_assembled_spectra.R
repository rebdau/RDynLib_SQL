#' Create assembled spectra from an MSn experiment
#'
#' This function takes a \linkS4class{Spectra} object (e.g., MSn data in `ftms`) 
#' and combines spectra belonging to the same MSn tree (`MSntreeID`). 
#' Intensities are summed for identical m/z values, and a new assembled
#' \linkS4class{Spectra} object is returned, with metadata propagated from 
#' the first MS2 spectrum in each tree.
#'
#' @param ftms A \linkS4class{Spectra} object containing MS2–MSn spectra.
#'
#' @return A new \linkS4class{Spectra} object with assembled spectra,
#'         including columns `mz`, `intensity`, `MSntreeID`, and `spectrum.type = "assembled"`.
#'
#' @importFrom Spectra Spectra spectraData peaksData setBackend MsBackendMemory
#' @export
create_assembled_spectra <- function(ftms) {
  # ftms must be a Spectra object
  # Converting the backend to MsBackendMemory
  # so we could modify the peaksData
  ftms <- setBackend(ftms, MsBackendMemory())
  meta <- spectraData(ftms)
  pd_list <- peaksData(ftms)
  tree_ids <- unique(meta$MSntreeID)
  
  assembled_list <- vector("list", length(tree_ids))
  k <- 0
  
  for (tree_id in tree_ids) {
    idx <- which(meta$MSntreeID == tree_id)
    ms_levels <- meta$msLevel[idx]
    
    # Use already preprocessed peaks
    spectra_list <- pd_list[idx]
    
    # Combine all spectra by summing intensities for matching m/z
    merged_spec <- combine_spectra(spectra_list)
    
    # Use metadata from the first MS2 spectrum
    ms2_row <- idx[ms_levels == 2][1]
    ms2_meta <- meta[ms2_row, , drop = FALSE]
    
    # Create new assembled spectrum
    df <- ms2_meta
    # Mark it as assembled
    df$spectrum.type <- "assembled"
    df$MSntreeID <- tree_id
    # Copy the same mz and intensity (already cleaned)
    df$mz <- list(merged_spec[, "mz"])
    df$intensity <- list(merged_spec[, "intensity"])
    # Remove rows where msLevel is NA
    df <- df[!is.na(df$msLevel), , drop = FALSE]
    k <- k + 1
    # Convert the modified metadata+peaks into a new
    # Spectra object and append it to assembled_list
    assembled_list[[k]] <- Spectra(df)
  }
  
  if (length(assembled_list) == 0) {
    warning("No assembled spectra found!")
    return(NULL)
  }
  
  do.call(c, assembled_list)
}