#' @title Filter targeted MS2 comparison results
#'
#' @description
#' Helper function that evaluates the result of a targeted MS2 spectral
#' comparison performed by `targMS2comp()`. A match is retained when either
#' the forward or reverse ion criterion is satisfied, both spectra contain
#' more than one ion, and the spectral similarity score (`dot_ions`) is
#' greater than 0.9.
#'
#' @param dbkey1 Identifier of the first database entry passed to
#'   `targMS2comp()`.
#' @param dbkey2 Identifier of the second database entry passed to
#'   `targMS2comp()`.
#' @param subDB Database or database subset passed to `targMS2comp()`.
#' @param AnalMS MS analysis information passed to `targMS2comp()`.
#'
#' @return A numeric spectral similarity score (`dot_ions`) when the matching
#'   criteria are satisfied and the score is greater than `0.9`. Returns `0`
#'   when the matching criteria are satisfied but `dot_ions` is `NaN`, and
#'   `NA` when the comparison does not satisfy the required criteria.
#'
#' @noRd
out.targMS2comp <- function(dbkey1, dbkey2, subDB, AnalMS) {
	out<-targMS2comp(dbkey1,dbkey2,subDB,AnalMS)
	outcm<-((out$forward_ions==1)|(out$reverse_ions==1))&(out$prec_nr_ions>1)&(out$cmp_nr_ions>1)
	if (is.na(outcm)){
	 return(NA)
	}else{
	 if (outcm==TRUE){
	  if (is.nan(out$dot_ions)){
	   return(0)
	  }else{
	   if (out$dot_ions>0.9){
	    return(out$dot_ions)
	   }else{
	    return(NA)
	   }
	  }
	 }else{
	  return(NA)
	 }
	}
}
