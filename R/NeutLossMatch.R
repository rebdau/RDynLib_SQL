#' @title Match neutral losses to product ions
#'
#' @description
#' Helper function that compares calculated neutral-loss values with a
#' reference list of neutral losses within a specified mass tolerance.
#' For each match, the corresponding neutral-loss mass and associated
#' annotation information are printed.
#'
#' @param neutprod.num A numeric vector containing the calculated neutral-loss
#'   values to be matched.
#' @param neutloss A matrix or data frame containing reference neutral-loss
#'   information. The first column contains the neutral-loss masses, while the
#'   second and third columns contain the corresponding annotation information.
#' @param err Numeric. Mass tolerance used to define the matching window around
#'   each reference neutral-loss mass.
#'
#' @return This function is called for its side effect of printing matching
#'   neutral losses and their associated annotation information. It does not
#'   explicitly return a value.
#'
#' @author Ahlam Mentag
#'
#' @noRd
NeutLossMatch <- function(neutprod.num, neutloss, err) {
	i=1
	repeat{
	 lowmass<-neutloss[i,1]-err
	 highmass<-neutloss[i,1]+err
	 j=1
	 repeat{
	  if((neutprod.num[j]>=lowmass)&(neutprod.num[j]<=highmass)){
	   infostr<-paste(neutprod.num[j],"Da --",neutloss[i,2],
				"--",neutloss[i,3],sep=" ")
	   print(infostr)
	  }
	  if(j==length(neutprod.num))break
	  j=j+1
	 }
	 if(i==dim(neutloss)[1])break
	 i=i+1
	}
	writeLines("\n")
}

