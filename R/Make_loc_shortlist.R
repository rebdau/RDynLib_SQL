#' @title Create a local shortlist around a compound
#'
#' @description
#' Helper function that creates a local shortlist of compounds surrounding a
#' selected compound in retention-time order. The compounds are first ordered
#' by retention time, and a specified number of neighboring entries before and
#' after the selected compound are retained. The selected compound itself is
#' excluded from the resulting shortlist.
#'
#' @param COMPID Integer. Compound identifier of the compound around which the
#'   local shortlist is created.
#' @param reg Numeric. Size of the retention-time neighborhood used to select
#'   surrounding compounds. Half of `reg` entries are selected before and half
#'   after the compound of interest.
#' @param ft.exp A matrix or data frame containing compound information. The
#'   first column is expected to contain retention times and the third column
#'   compound identifiers.
#'
#' @return A matrix or data frame containing the neighboring compounds around
#'   `COMPID`, ordered by retention time and excluding the compound of
#'   interest.
#'
#' @noRd
Make_loc_shortlist <- function(COMPID, reg, ft.exp) {
	#order on retention time
	ft.o<-ft.exp[order(ft.exp[,1]),] #reg is equal to nr of rows in ft.o or 1 less
	i<-which(ft.o[,3]==COMPID) #middle of chromatogram
	regdiv<-reg/2
	i.min<-i-regdiv
	if(i.min<0){
	 stop("Sorry,your compound is eluting very early. You have to decrease 'reg'")
	}
	i.max<-i+regdiv
	if(i.max>dim(ft.o)[1]){
	 stop("Sorry,your compound is eluting very late. You have to decrease 'reg'")
	}
	rng<-c(i.min:(i-1),(i+1):i.max) 
	ft.sh<-ft.o[rng,]	#all entries except peak in middle of chromatogram
	return(ft.sh)
}

