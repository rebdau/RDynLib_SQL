#' @title Select the best compound match
#'
#' @description
#' Helper function that selects the best candidate compound match and stores
#' its compound identifier in an association table. If only one candidate is
#' available, it is selected directly. If multiple candidates are available,
#' the expected retention time is calculated from the regression parameters,
#' and the candidate with the smallest absolute difference from the predicted
#' retention time is selected.
#'
#' @param pres A matrix or data frame containing candidate compound matches.
#'   The first column contains retention times and the third column contains
#'   compound identifiers.
#' @param Assoc A matrix or data frame containing compound associations to be
#'   updated.
#' @param time.ftng Numeric. Retention time of the reference compound used to
#'   calculate its predicted retention time.
#' @param COMPID.ftng Integer. Row index in `Assoc` corresponding to the
#'   reference compound for which the selected match is stored.
#' @param rg A numeric vector containing the regression parameters. The first
#'   two elements represent the intercept and slope used to predict retention
#'   time.
#'
#' @return The updated `Assoc` object with the compound identifier of the
#'   selected candidate stored in the second column.
#'
#' @author Ahlam Mentag
#'
#' @noRd
Sort_comp_matches <- function(pres, Assoc, time.ftng, COMPID.ftng, rg) {
	if(dim(pres)[1]==1){
	 Assoc[COMPID.ftng,2]<-pres[1,3]
		# COMPID.ftng refers to COMPID entry, and is used here to
		# select the row (!) in the Assoc dataframe rather than the 
		# row corresponding with the COMPID itself.   
	}else{
	 pred_time=rg[1]+rg[2]*time.ftng
	 timediff<-abs(pres[,1]-pred_time)
	 pres.corrected<-cbind(pres,timediff)
	 o<-order(pres.corrected[,5])
	 pres.corr<-pres.corrected[o,]
	 Assoc[COMPID.ftng,2]<-pres.corr[1,3]
	}
	return(Assoc)
}
