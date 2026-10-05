#' @title Convert MSn trees to pseudo-MS2 spectra
#'
#' @description
#' Converts MSn spectral trees into pseudo-MS2 spectra by combining the
#' fragment ions from all spectra belonging to the same MSn tree. Fragment
#' m/z values are first rounded to unit mass independently for each spectrum,
#' and duplicated rounded m/z values are resolved by retaining the most
#' intense peak. The processed spectra belonging to the same `MSntreeID` are
#' then merged, with intensities of identical m/z values summed. One
#' pseudo-MS2 spectrum is generated for each MSn tree containing an MS2
#' spectrum.
#'
#' Metadata for each resulting pseudo-MS2 spectrum are inherited from the
#' first MS2 spectrum of the corresponding MSn tree, and the spectrum type
#' is set to `"pseudo_MS2"`.
#'
#' @param sps A `Spectra` object containing MSn spectra. The object must
#'   contain an `MSntreeID` variable identifying spectra belonging to the
#'   same MSn tree and an `msLevel` variable identifying the fragmentation
#'   level.
#'
#' @return A `Spectra` object containing one pseudo-MS2 spectrum per MSn tree.
#'   The resulting spectra contain the combined fragment ions from all
#'   available fragmentation levels of the corresponding tree. If no
#'   pseudo-MS2 spectra can be generated, `NULL` is returned with a warning.
#'
#' @author Ahlam Mentag
#'
#' @noRd
msntree_to_pseudoMS2 <- function(sps) {  
  # Ensure memory backend
  sps <- setBackend(sps, MsBackendMemory())
  
  
  # Helper:  rounding
  
  round_perl <- function(number) {
    floor(number + 0.5)
  }
  
  # Helper: round mz and keep most intense duplicated peak
  
  round_spectrum <- function(x) {
    
    if (nrow(x) == 0)
      return(x)
    
    x[, 1] <- round_perl(x[, 1])
    
    df <- data.frame(
      mz = x[, 1],
      intensity = x[, 2]
    )
    
    df <- df |>
      dplyr::group_by(mz) |>
      dplyr::slice_max(
        order_by = intensity,
        n = 1,
        with_ties = FALSE
      ) |>
      dplyr::ungroup() |>
      dplyr::arrange(mz)
    
    as.matrix(df)
  }
  
  
  # Helper: merge spectra from one tree
  
  combine_spectra <- function(spectra_list) {
    
    if (length(spectra_list) == 1)
      return(spectra_list[[1]])
    
    merged <- do.call(rbind, spectra_list)
    
    merged <- data.frame(
      mz = merged[, 1],
      intensity = merged[, 2]
    ) |>
      dplyr::group_by(mz) |>
      dplyr::summarise(
        intensity = sum(intensity),
        .groups = "drop"
      ) |>
      dplyr::arrange(mz)
    
    as.matrix(merged)
  }
  
  
  # Extract data
  
  meta <- spectraData(sps)
  pd_list <- peaksData(sps)
  
  # Round each spectrum independently
  pd_list <- lapply(pd_list, round_spectrum)
  
  tree_ids <- unique(meta$MSntreeID)
  
  assembled_list <- vector("list", length(tree_ids))
  n_out <- 0
  
  
  # Assemble one spectrum per tree
  for (tree_id in tree_ids) {
    
    idx <- which(meta$MSntreeID == tree_id)
    
    if (length(idx) == 0)
      next
    
    merged_spec <- combine_spectra(pd_list[idx])
    
    ms2_rows <- idx[meta$msLevel[idx] == 2]
    
    if (length(ms2_rows) == 0)
      next
    
    ms2_meta <- meta[ms2_rows[1], , drop = FALSE]
    
    ms2_meta$spectrum.type <- "pseudo_MS2"
    ms2_meta$MSntreeID <- tree_id
    
    ms2_meta$mz <- list(merged_spec[, "mz"])
    ms2_meta$intensity <- list(merged_spec[, "intensity"])
    
    rownames(ms2_meta) <- NULL
    
    n_out <- n_out + 1
    assembled_list[[n_out]] <- Spectra(ms2_meta)
  }
  
  assembled_list <- assembled_list[seq_len(n_out)]
  
  if (length(assembled_list) == 0) {
    warning("No pseudo_MS2 spectra could be created.")
    return(NULL)
  }
  
  do.call(base::c, assembled_list)
}