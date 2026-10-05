#' @title Match product ions
#'
#' @description
#' Helper function that compares observed product-ion m/z values with a
#' reference list of product ions within a specified mass tolerance. For each
#' match, the observed m/z value and the corresponding reference annotation
#' information are printed.
#'
#' @param prod_ion.num A numeric vector containing the observed product-ion
#'   m/z values to be matched.
#' @param prdion A matrix or data frame containing reference product-ion
#'   information. The first column contains the reference m/z values, while
#'   the second and third columns contain the corresponding annotation
#'   information.
#' @param err Numeric. Mass tolerance used to define the matching window
#'   around each reference product-ion m/z value.
#'
#' @return This function is called for its side effect of printing matching
#'   product ions and their associated annotation information. It does not
#'   explicitly return a value.
#'
#' @noRd
ProdIonMatch <- function(prod_ion.num, prdion, err) {
	i=1
	repeat{
	 lowmz<-prdion[i,1]-err
	 highmz<-prdion[i,1]+err
	 j=1
	 repeat{
	  if((prod_ion.num[j]>=lowmz)&(prod_ion.num[j]<=highmz)){
	   infostr<-paste("m/z",prod_ion.num[j],"--",prdion[i,2],
				"--",prdion[i,3],sep=" ")
	   print(infostr)
	  }
	  if(j==length(prod_ion.num))break
	  j=j+1
	 }
	 if(i==dim(prdion)[1])break
	 i=i+1
	}
	writeLines("\n")
}

