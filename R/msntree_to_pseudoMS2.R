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
    
    ms2_meta$spectrum.type <- "assembled"
    ms2_meta$MSntreeID <- tree_id
    
    ms2_meta$mz <- list(merged_spec[, "mz"])
    ms2_meta$intensity <- list(merged_spec[, "intensity"])
    
    rownames(ms2_meta) <- NULL
    
    n_out <- n_out + 1
    assembled_list[[n_out]] <- Spectra(ms2_meta)
  }
  
  assembled_list <- assembled_list[seq_len(n_out)]
  
  if (length(assembled_list) == 0) {
    warning("No assembled spectra could be created.")
    return(NULL)
  }
  
  do.call(base::c, assembled_list)
}