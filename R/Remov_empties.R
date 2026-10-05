#' @title Remove unmatched compounds based on m/z
#'
#' @description
#' Helper function that filters compounds from one experiment by retaining
#' only entries for which a corresponding m/z value is found in a second
#' experiment within a specified mass tolerance. Both datasets are first
#' ordered by m/z to enable sequential matching. If fewer than 30 percent
#' of the expected entries are retained, a warning message is printed.
#'
#' @param err Numeric. Mass tolerance used to define the m/z matching window.
#' @param syn.exp A matrix or data frame containing compounds from the
#'   comparison experiment. The second column contains the m/z values used
#'   for matching.
#' @param ft.sh A matrix or data frame containing compounds to be filtered.
#'   The second column contains the m/z values used for matching.
#' @param reg Numeric. Reference number of entries used to evaluate whether
#'   at least 30 percent of the entries are retained after filtering.
#'
#' @return The filtered `ft.sh` object containing only entries with a
#'   corresponding m/z value in `syn.exp` within the specified mass
#'   tolerance. The returned object is ordered by its first column.
#'
#' @noRd
Remov_empties <- function(err, syn.exp, ft.sh, reg) {
	ft.sh<-ft.sh[order(ft.sh[,2]),] #order on m/z value
	syn.o<-syn.exp[order(syn.exp[,2]),] #order on m/z value
	i=1
	j=1
	repeat{
	 lb=ft.sh[i,2]-err #lower threshold QTOF m/z value
	 ub=ft.sh[i,2]+err #higher threshold QTOF m/z value
	 while (j<=length(syn.o[,2])) {
	  if (syn.o[j,2]>lb&syn.o[j,2]<ub){
		#if a QTOF row is found of which the m/z is in the selected
		#m/z window, don't do anything, just select the next FT m/z
		#value. Yet, if the QTOF row is not the first row, lower the
		#QTOF starting row so that the current QTOF-based m/z value
		#will be refound if the next FT-based m/z value is an isomer
		#of the current FT-based m/z value
	   if (j>1) j=j-1
	   i=i+1
	   break
	  }
	  if (syn.o[j,2]>=ub){
		#if QTOF row has higher m/z than FT selected m/z, remove
		#FT selected m/z row, break searching the QTOF DynLib
		#and go to the next FT m/z value (for which i remains 1)
	   ft.sh<-data.frame(ft.sh[-i,], check.names = FALSE) 
	   rownames(ft.sh)=1:(dim(ft.sh)[1])
	   if (j>1) j=j-1
	   break
	  }
		#if QTOF row m/z value is lower than the lower m/z threshold,
		#go to next QTOF row
	  j=j+1
	 }
	 if (i>dim(ft.sh)[1]) break
	 if (j>dim(syn.o)[1]) {
	  ft.sh<-ft.sh[1:i-1,]
	  break
	 }
	}
	if(dim(ft.sh)[1]<0.3*reg){
	 print("Less than 30% of the entries are retained.")
	}
	ft.sh<-ft.sh[order(ft.sh[,1]),]
	return(ft.sh)
}

