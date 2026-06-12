#' @title searching for combination and in-source fragments within the FGs 
#' using CID spectral data. 
#'
#' @description
#' The `XCMSgrouping2_SQL()` function searches for in-source fragments within 
#' the FGs using CID spectral data via the nested AutomFindingISF()
#' function In addition, the function searches for features that represent  
#' combinations of other features present in the proper and in the neighboring 
#' FGs via the SearchCombinat() function. The function takes the xcms resulting 
#' object from the XCMSgrouping() function as input and creates an upgraded 
#' xcms object.
#' 
#' @param x `XcmsExperiment` object resulting from the XCMSgrouping() function.
#' 
#' @param adjPeak `numeric()` number of adjacent feature groups in the 
#' chromatogram to search for possible combinations of features, 1 by default.
#' 
#' @param mode `character(1)` is the mode of the experiment, either neg or pos.
#' 
#' @return An update `XcmsExperiment` object with new and updated columns in 
#'         featureDefinitions .
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
XCMSgrouping2_SQL <- function(x, adjPeak = 1, mode = NULL) {
  

  dat <- featureDefinitions(x)
  dat$feature_id <- rownames(dat)
  
  dat$name <- .create_name(dat$mzmed, dat$rtmed)
  

  s <- featureSpectra(x)
  
  spd    <- Spectra::spectraData(s)
  ms2.ft <- Spectra::peaksData(s)
  
  stopifnot("feature_id" %in% colnames(spd))
  

  if (is.null(mode)) {
    mode <- ifelse(spd$polarity[1] < 0, "neg", "pos")
  }
  
  # here we detect the isf by feature groups
  fg_list <- unique(dat$CON.new)
  
  XCMS_ISF <- do.call(
    rbind,
    lapply(fg_list, FindingISF_SQL,
           XCMS = dat,
           spd = spd,
           ms2.ft = ms2.ft)
  )
  
  XCMS4 <- ConvertXCMS4(XCMS_ISF)
  XCMS4 <- SearchCombinat_SQL(XCMS4, adjPeak, mode)

  fd <- featureDefinitions(x)
  
  new_cols <- setdiff(colnames(XCMS4), colnames(fd))
  
  for (col in new_cols) {
    fd[[col]] <- XCMS4[[col]][match(rownames(fd), XCMS4$feature_id)]
  }
  
  featureDefinitions(x) <- fd
  
  x
}


#'helper function to a create a name for the features.
#'
.create_name <- function(mzmed, rtmed) {
  paste0(
    "M", round(mzmed, 0),
    "T", round(rtmed, 0)
  )
}
