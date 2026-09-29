#'@title create an mgf file from a spectra object
#'
#'@param sps spectra object
#'@param mgf_file output file
#'
#'@return an mgf file
#'
#'@author Ahlam Mentag
#'
#'@export
spectra_to_mgf <- function(sps, mgf_file) {
  
  # Extract peaks and metadata
  peaks_list <- Spectra::peaksData(sps)
  meta_data  <- Spectra::spectraData(sps)
  
  # Open output file
  con <- file(mgf_file, "w")
  on.exit(close(con), add = TRUE)
  
  # Write spectra
  for (i in seq_along(peaks_list)) {
    
    peak_matrix <- peaks_list[[i]]
    
    # Skip empty spectra
    if (is.null(peak_matrix) || nrow(peak_matrix) == 0) {
      next
    }
    
    cat("BEGIN IONS\n", file = con)
    
    cat(
      paste0("FEATURE_ID=", i, "\n"),
      file = con
    )
    
    # Precursor m/z
    pepmass <- meta_data$precursorMz[i]
    
    if (is.na(pepmass)) {
      pepmass <- ""
    }
    
    cat(
      paste0("PEPMASS=", pepmass, "\n"),
      file = con
    )
    
    cat("CHARGE=1-\n", file = con)
    cat("MSLEVEL=2\n", file = con)
    
    # Write m/z and intensity
    apply(
      peak_matrix,
      1,
      function(x) {
        cat(x[1], x[2], "\n", file = con)
      }
    )
    
    cat("END IONS\n\n", file = con)
  }
  
  message("MGF file written to: ", mgf_file)
  
  invisible(mgf_file)
}