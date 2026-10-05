#' @title Remove retention-time outliers
#'
#' @description
#' Helper function that removes alignment pairs with retention-time differences
#' outside a specified range around the mean. The lower and upper limits are
#' calculated as the mean retention-time difference plus or minus a specified
#' number of standard deviations.
#'
#' @param LCal A matrix or data frame containing LC alignment results. The
#'   sixth column contains the retention-time differences used for outlier
#'   detection.
#' @param rng Numeric. Number of standard deviations from the mean used to
#'   define the lower and upper limits for retaining alignment pairs.
#'
#' @return The filtered `LCal` object with retention-time outliers removed.
#'
#' @noRd
RemoveOutliers <- function(LCal, rng) {
	m.trdiff<-mean(as.numeric(LCal[,6]))
	sd.trdiff<-sd(as.numeric(LCal[,6]))
	diff.low=m.trdiff-rng*sd.trdiff
	diff.high=m.trdiff+rng*sd.trdiff
	i=1
	while(i<=dim(LCal)[1]){
	 if((as.numeric(LCal[i,6])<diff.low)|(as.numeric(LCal[i,6])>diff.high)){
	  LCal<-LCal[-i,]
	  next
	 }
	 i=i+1
	}
	return(LCal)
}