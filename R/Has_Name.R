#' @title Identify compounds with valid names
#'
#' @description
#' Helper function that checks whether compound names are considered valid
#' based on a set of values indicating missing or unavailable annotations.
#' Names matching one of the specified values, or beginning with a specified
#' exclusion character, are assigned `0`. Otherwise, the corresponding value
#' from `sel` is retained.
#'
#' @param name.cha A character vector containing compound names to evaluate.
#' @param No.name A character vector defining values or prefixes that indicate
#'   a missing or unavailable compound name, such as `"NULL"` or `"!"`.
#' @param sel A vector containing the values to retain for compounds with
#'   valid names.
#'
#' @return An integer vector in which entries without a valid compound name
#'   are assigned `0`, while entries with a valid name contain the
#'   corresponding value from `sel`.
#'
#' @noRd
Has_Name <- function(name.cha, No.name, sel) {
	HasName<-as.integer()
	j=1
	while (j<=length(name.cha)) {
	 Ent.name<-as.character(c(substr(name.cha[j],1,1),name.cha[j]))
	 if (length(which(Ent.name%in%No.name))!=0){ # No.name is NULL or !
	 # Ent.name contains first character of the name to check for "!"
	 # and full name to check for "NULL"
	  HasName<-append(HasName,0)
	 }else{
	  HasName<-append(HasName,sel[j])
	 }
	 j=j+1
	}
	return(HasName)
}
