#' @title Check XCMS subset sizes
#'
#' @description
#' Helper function that checks the number of rows in consecutive XCMS subsets.
#' If a subset contains more than 3000 rows, the user is warned and prompted
#' to choose whether to continue the analysis.
#'
#' @param stp A numeric vector defining the row boundaries of the XCMS subsets.
#'
#' @return If at least one subset contains more than 3000 rows, returns the
#'   user's response (`"1"` or `"0"`) indicating whether to continue.
#'   Otherwise, the function prints a message and returns `NULL`.
#'
#' @author Ahlam Mentag
#'
#' @noRd
CheckBreakRows<-function(stp){
	stpm<-outer(stp,stp,"-")
	stp.l<-as.integer(rep(NA,length(stp)))
	i=1
	repeat{
	 j=i+1
	 stp.l[i]<-stpm[j,i]
	 if(j==length(stp))break
	 i=i+1
	}
	if(max(stp.l,na.rm=T)>3000){
	 print("one XCMS subset counts more than 3000 rows, consider using
	  a smaller tR value.")
	 CBR<-readline("Enter 1 or 0 whether you want to continue or not:")
	 return(CBR)
	}else{
	 print("all XCMS subsets contain less than 3000 rows")
	}
}
