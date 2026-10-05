#' @title Prepare compound data for CSPP name transfer
#'
#' @description
#' Helper function that extracts compound identifiers, compound names, and
#' conversion information from a DynLib compound table and formats them for
#' CSPP-based name transfer. Missing compound names and conversion values are
#' replaced with `"NULL"`.
#'
#' @param inp.x A data frame containing DynLib compound information. The
#'   function extracts the compound identifier, compound name, and conversion
#'   information from columns 1, 3, and 11, respectively.
#'
#' @return A data frame with three columns: `COMPID`, `COMPNAME`, and
#'   `CONVERSION`, with missing compound names and conversion values replaced
#'   by `"NULL"`.
#'
#' @author Ahlam Mentag
#'
#' @noRd
File_CSPPname <- function(inp.x) {
	msdet<-data.frame(inp.x[,1],inp.x[,3],inp.x[,11])
	msdet[,2]<-as.character(msdet[,2])
	msdet[,3]<-as.character(msdet[,3])
	colnames(msdet)<-c("COMPID","COMPNAME","CONVERSION")
	i=1
	while (i<=dim(msdet)[1]) {
	 if (is.na(msdet[i,2])) msdet[i,2]<-"NULL"
	 i=i+1
	}
	i=1
	while (i<=dim(msdet)[1]) {
	 if (is.na(msdet[i,3])) msdet[i,3]<-"NULL"
	 i=i+1
	}
	return(msdet)
}
