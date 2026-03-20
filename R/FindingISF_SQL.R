#' @title finding in-source fragments.
#' 
#' @description
#' finding whether a m/z feature is an in-source fragment (ISF) of another m/z
#' feature in a particular feature group. 
#'fg is the feature group number.
#'
#' @param XCMS the resulting `XcmsExperiment` object from the 
#'        `XCMSgrouping_SQL()` function.
#'        
#' @param fg set within the AutomFindingISF() function, represents the feature 
#'        group.
#'        
#' @param spd spectra Data of the `XcmsExperiment` object.
#' 
#' @param ms2.ft peaks data of the `XcmsExperiment` object.
#' 
#' @return An updated `XcmsExperiment` object with new annotation column 
#'         created inside FindingISF_SQL() function.
#'         
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


FindingISF_SQL <- function(fg, XCMS, spd, ms2.ft) {
  
  ## Features in this group
  XCMS.sel <- XCMS[XCMS$CON.new == fg, ]
  XCMS.sel$isf <- NA_character_
  
  if (nrow(XCMS.sel) == 0)
    return(XCMS.sel)
  
  ## Build MS2 fragments per feature
  ms2_by_feature <- build_ms2_by_feature(spd, ms2.ft)
  
  ## Loop over features
  for (i in seq_len(nrow(XCMS.sel))) {
    
    fid <- XCMS.sel$feature_id[i]
    parent_mz <- round(XCMS.sel$mzmed[i], 0)
    
    ## Skip if no MS2
    if (!fid %in% names(ms2_by_feature))
      next
    
    frags <- ms2_by_feature[[fid]]
    
    ## Remove parent ion if present
    frags <- frags[frags != parent_mz]
    
    ## check if any other feature matches these fragments
    mz_all <- round(XCMS.sel$mzmed, 0)
    
    if (any(mz_all %in% frags)) {
      XCMS.sel$isf[mz_all %in% frags] <- "ISF"
    }
  }
  
  XCMS.sel
}


build_ms2_by_feature <- function(spd, ms2.ft) {
  
  stopifnot("feature_id" %in% colnames(spd))
  stopifnot(length(ms2.ft) == nrow(spd))
  
  ## Only MS2 spectra with a feature_id
  idx <- which(!is.na(spd$feature_id))
  
  res <- split(
    lapply(ms2.ft[idx], function(p) round(p[, 1], 0)),
    spd$feature_id[idx]
  )
  
  ## Flatten: multiple spectra per feature to one vector
  res <- lapply(res, function(x) unique(unlist(x)))
  
  res
}



 

