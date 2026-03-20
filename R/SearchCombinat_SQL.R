#' @title Searching for features that represent combinations of other features. 
#'
#' @description
#' The `SearchCombinat_SQL()` function searches for features that represent 
#' combinations of other features present in the proper and in the neighboring 
#' FGs.
#' 
#' @param XCMS4 the resulting `XcmsExperiment` object from the `ConvertXCMS4()`
#' function.
#' 
#' @param adjPeak `numeric()` number of adjacent feature groups in the 
#' chromatogram to search for possible combinations of features, 1 by default.
#' 
#' @param mode `character(1)` is the mode of the experiment, either neg or pos.
#' 
#' @return An updated `XcmsExperiment` object with new columns in 
#'         featureDefinitions, "Comb" and "massMS1":
#'         
#'         - `Comb` : Ion/neutral combination annotation indicating pairs of
#'           features whose masses sum to another observed feature (formatted
#'           as `"mz1/mz2"`).
#'           
#'         - `massMS1` : Integer-encoded MS1 mass (m/z adjusted for protonation/
#'            deprotonation and scaled by 100) used for ion–neutral combination
#'            matching.
#'         
#' @import xcms
#'
#' @import MsExperiment
#'
#' @importFrom DBI dbDisconnect
#'
#' @import Spectra
#'
#' @import dplyr
#'
#' @author Ahlam Mentag
#'
#' @export  
SearchCombinat_SQL <- function(XCMS4, adjPeak, mode) {
  
  ## Proton mass correction
  prtn <- if (mode == "neg") 1.0073 else -1.0073
  
  # ## Ensure feature_id exists
  # if (!"feature_id" %in% colnames(XCMS4)) {
  #   XCMS4$feature_id <- rownames(XCMS4)
  # }
  # 
  ## Compute MS1 neutral mass proxy
  XCMS4$massMS1 <- as.integer(round((XCMS4$mzmed + prtn) * 100, 2))
  
  ## Initialize combination column
  XCMS4$Comb <- "NULL"
  
  ## Peak group boundaries
  EndGrp   <- max(XCMS4$CorrPkGrp)
  StartGrp <- 1 + adjPeak
  
  i <- StartGrp
  while (i <= EndGrp) {
    
    ## Subset nearby peak groups
    subXCMS <- XCMS4[
      XCMS4$CorrPkGrp %in% (i - adjPeak):(i + adjPeak),
    ]
    
    ## Indices of group of interest
    PkGrpOfInt <- which(subXCMS$CorrPkGrp == i)
    
    ## Pairwise sums of masses
    ion_neut <- outer(subXCMS$massMS1, subXCMS$massMS1, "+")
    
    ## Where summed masses exist as MS1 masses
    ion_neut_ind <- which(ion_neut %in% subXCMS$massMS1)
    
    if (length(ion_neut_ind) == 0) {
      i <- i + 1
      next
    }
    
    n <- nrow(subXCMS)
    
    row_ind <- ceiling(ion_neut_ind / n)
    col_ind <- ifelse(
      ion_neut_ind %% n == 0,
      n,
      ion_neut_ind %% n
    )
    
    INcomp <- match(
      ion_neut[cbind(row_ind, col_ind)],
      subXCMS$massMS1
    )
    
    ## Build combination strings
    Combin <- rep(NA_character_, n)
    Combin[INcomp] <- paste0(
      round(subXCMS$mzmed[row_ind], 2),
      "/",
      round(subXCMS$mzmed[col_ind], 2)
    )
    
    ## Restrict to peak group of interest
    CombSel <- intersect(which(!is.na(Combin)), PkGrpOfInt)
    
    ## Assign back to main table
    if (length(CombSel) > 0) {
      XCMS4$Comb[
        match(subXCMS$feature_id[CombSel], XCMS4$feature_id)
      ] <- Combin[CombSel]
    }
    
    i <- i + 1
  }
  
  return(XCMS4)
}










 