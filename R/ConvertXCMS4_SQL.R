#' @title Performing a re-indexing and re-ordering of feature groups based  
#'        on retention time.
#' @description
#' The `ConvertXCMS4_SQL()` computes the mean retention time for each 
#' feature group, orders features accordingly, and generates a corrected, 
#' sequential peak-group index for downstream ion–neutral and combinatorial 
#' analyses.
#' 
#' @param XCMS4 the resulting `XcmsExperiment` object from the 
#'        `AutomFindingISF_SQL()` function.
#' 
#' @return An updated `XcmsExperiment` object with new columns in 
#'         featureDefinitions, "MeAb" and "CorrPKGrp":
#'         
#'         - `MeAb` : Mean retention time (rtmed) of the feature’s peak group 
#'                  (CON.new).
#'           
#'         - `massMS1` : Corrected sequential peak-group index after RT-based 
#'                       reordering.
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
ConvertXCMS4 <- function(XCMS4) {
  
  MeAb <- rep(NA_real_, nrow(XCMS4))
  
  # Column names
  pkcol <- "CON.new"   # peak group
  rtcol <- "rtmed"     # median RT
  subcol <- "subgrp"   # subgroup
  
  # Singletons
  Singlet <- which(as.integer(XCMS4[[pkcol]]) == 1)
  MeAb[Singlet] <- XCMS4[[rtcol]][Singlet]
  
  # Other peak groups
  i <- 2
  while (length(which(XCMS4[[pkcol]] %in% i)) != 0) {
    PkGrp <- which(XCMS4[[pkcol]] %in% i)
    RtAve <- mean(XCMS4[[rtcol]][PkGrp])
    MeAb[PkGrp] <- RtAve
    i <- i + 1
  }
  
  # Add MeAb column
  XCMS4$MeAb <- MeAb
  
  # Order by MeAb
  XCMS4 <- XCMS4[order(XCMS4$MeAb), ]
  
  # Compute CorrPkGrp
  CorrPkGrp <- rep(NA_integer_, nrow(XCMS4))
  CorrPkGrp[1:2] <- c(1, 2)
  
  i <- 2
  while (i < nrow(XCMS4)) {
    j <- i + 1
    if ((XCMS4[[subcol]][j] == 1) | (XCMS4[[subcol]][j] != XCMS4[[subcol]][i])) {
      CorrPkGrp[j] <- CorrPkGrp[i] + 1
    } else {
      CorrPkGrp[j] <- CorrPkGrp[i]
    }
    i <- i + 1
  }
  
  XCMS4$CorrPkGrp <- CorrPkGrp
  
  return(XCMS4)
}



 
