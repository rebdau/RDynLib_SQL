#' @title Grouping and annotating the features of an `XcmsExperiment` object. 
#'
#' @description
#' The `XCMSgrouping_SQL()` takes as input a preprocessed `XcmsExperiment` 
#' object, then extract the featureDefinitions for mzmed, rtmed, and name columns
#' also we call featureValues to get the individual intensities of the features 
#' in each sample columns from it. Hence we add an fc column for the sample 
#' groups. Then we call a set of functions to perform the grouping and
#' the annotation of the features, and add those information as additional
#' columns to featureDefinitions.
#' 
#' @param tR `numeric(1)` minimum break-up retention time difference: determines
#'  the minimum difference in retention time between two subsequent rows 
#'  (i.e., features) prior to considering the row as a candidate for splitting 
#'  the ‘featureDefinitions()’-derived spreadsheet (\link{BreakRowsXCMS_SQL}), 
#'  1 sec by default.
#'
#' @param y10 maximum retention time deviation for feature alignment across 
#'  chromatograms, 1 sec by default.
#' 
#' @param y20 minimum feature abundance-based Pearson correlation for feature  
#'  alignment across chromatograms, 0.8 by default 
#' 
#' @param y11 monoisotopic mass of the buffer adduct, 46 Da (formic acid) 
#'  by default.
#' 
#' @param x `XcmsExperiment` object with the results from xcms preprocessing
#' 
#' @return An update `XcmsExperiment` object with new columns in 
#'         featureDefinitions of the grouping and annotation information.
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
XCMSgrouping_SQL <- function(x, tR = 1, y10 = 1, y20 = 0.8, y11 = 46.0) {
  
  t0 <- Sys.time()
  message("[XCMSgrouping_SQL] Start: ", t0)
  
  ## Extract original feature information from the xcmsExpeiment object
  message("[1/9] Extracting featureDefinitions()")
  dat <- featureDefinitions(x)
  
  rtmed <- dat$rtmed
  mzmed <- dat$mzmed
  dat$name <- .create_name(mzmed, rtmed)
  
  message("The Number of features is: ", nrow(dat))
  
  ## Sample groups
  message("[2/9] Extracting sample groups (fc)")
  fc <- sampleData(x)$fc
  
  ## Feature intensities
  message("[3/9] Extracting featureValues()")
  ints <- featureValues(x)
  
  ## in case there is just one sample
  if (is.null(dim(ints))) {
    ints <- matrix(ints, ncol = 1,
                   dimnames = list(names(ints), "sample1"))
  }
  
  ## Align intensities using rownames
  message("[4/9] Aligning intensities to featureDefinitions")
  fd_rn <- rownames(featureDefinitions(x))
  fv_rn <- rownames(ints)
  ints <- ints[match(fd_rn, fv_rn), , drop = FALSE]
  
  ## Break rows
  message("[5/9] Computing breakpoints (BreakRowsXCMS_SQL)")
  stp <- BreakRowsXCMS_SQL(dat, tR, rtmed)
  
  CBR <- CheckBreakRows(stp)
  message("      Breakpoints found: ", length(stp) - 1)
  
  if (CBR == 0) {
    message("      No valid breakpoints → returning original object")
    return(x)
  }
  
  ## Processing info
  message("[6/9] Loading processing info (LoadProc_SQL)")
  procfc <- LoadProc_SQL(x)
  fc <- procfc[[2]]
  
  ## Main grouping step (MOST EXPENSIVE)
  message("[7/9] Starting AutomFeatureGroup_SQL()")
  t_group <- Sys.time()
  
  XCMS <- AutomFeatureGroup_SQL(dat, rtmed, ints, stp, y10, y20, fc)
  
  message("      AutomFeatureGroup_SQL finished in ",
          round(difftime(Sys.time(), t_group, units = "secs"), 2), " sec")
  
  ## Generalize CON
  message("[8/9] Generalizing consensus groups")
  CON.new <- GeneralizeConXCMS(XCMS)
  XCMS$CON.new <- CON.new
  
  ## Sort by CON.new
  XCMS <- XCMS[order(XCMS$CON.new), ]
  
  ## Secondary grouping
  message("[9/9] Secondary RT-based subgrouping")
  XCMS <- FeatureSubgroupXCMS_SQL(XCMS, y10)
  
  ## Generalize subgroups
  message("      Generalizing subgroups")
  XCMS <- GeneralizeSubgroupXCMS(XCMS)
  
  ## Feature annotation
  message("      Feature annotation")
  XCMS <- FeatureAnnotXCMS_SQL(XCMS, y11, fc, ints)
  
  ## Add new columns back to featureDefinitions
  message("      Writing results back to featureDefinitions")
  fd <- featureDefinitions(x)
  new_cols <- setdiff(colnames(XCMS), colnames(fd))
  for (col in new_cols) {
    fd[[col]] <- XCMS[[col]]
  }
  featureDefinitions(x) <- fd
  
  message("[XCMSgrouping_SQL] Finished in ",
          round(difftime(Sys.time(), t0, units = "secs"), 2), " sec")
  
  return(x)
}



#'helper function to a create a name for the features.
#'
.create_name <- function(mzmed, rtmed) {
  paste0(
    "M", round(mzmed, 0),
    "T", round(rtmed, 0)
  )
}
