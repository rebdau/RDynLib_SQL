#' @title Generalize XCMS conversion identifiers
#'
#' @description
#' Helper function that converts the conversion identifiers stored in the
#' last column of an XCMS-derived object into a generalized sequential
#' numbering scheme. Repeated conversion identifiers are assigned to the same
#' group, while increasing identifiers define new groups.
#'
#' @param XCMS A matrix or data frame containing XCMS-derived data. The last
#'   column is expected to contain the conversion identifiers to be
#'   generalized.
#'
#' @return A numeric vector containing the generalized conversion identifiers
#'   for all rows of `XCMS`.
#'
#'
#' @noRd
GeneralizeConXCMS <- function(XCMS) {
	CON<-XCMS[,dim(XCMS)[2]]
	tim<-length(CON)-1
	CON.new<-c(1,rep(NA,tim))
	a=1
	i=2
	while(i<=length(CON)){
	 if(CON[i]==1){
	  CON.new[i]<-1
	 }else{
	  if(CON[i]==CON[i-1]){
	   CON.new[i]<-a
	  }else{
	   if(CON[i]>CON[i-1]){
	    a=a+1
	    CON.new[i]<-a
	   }
	  }
	 }
	 i=i+1
	}
	return(CON.new)
}
