#' @title Transfer compound names between CSPP-associated features
#'
#' @description
#' Helper function that propagates compound annotations between substrate and
#' product features linked by CSPP relationships. When the substrate is
#' annotated, its name is transferred to an unannotated product feature using
#' a `CSPP` prefix. Conversely, when the product is annotated, its name is
#' transferred to an unannotated substrate feature using a `CSPPr` prefix.
#'
#' @param msdet A data frame containing compound information, including
#'   compound names and CSPP relationships between substrate and product
#'   features.
#'
#' @return The updated `msdet` data frame with compound names propagated
#'   between CSPP-associated substrate and product features where applicable.
#'
#' @noRd
CSPPnameTransfer<-function(msdet){
	i=1
	# When CSPP 'substrate' m/z feature is known
	repeat{
	 No.name=as.character(c("NULL","!","CSPP")) # enter new name for the product m/z feature
	 # if one of these strings is present in the COMPNAME of the product m/z feature 
	 No.name2=as.character(c("NULL","!")) # enter new name for the product m/z feature 
	 # if none of these strings are present in the COMPNAME of the substrate m/z feature
	 if (msdet[i,3]!="NULL") {
	  sub_nam<-msdet[i,2]
	  Ent.namsub<-as.character(c(substr(sub_nam,1,1),sub_nam)) # character vector 
	  # containing the first character of the COMPNAME of the substrate m/z feature and 
	  # the full COMPNAME
	  cspp<-strsplit(msdet[i,3]," ")
	  prod_compid<-as.integer(cspp[[1]][3])
	  prod_nam<-msdet[prod_compid,2]
	  Ent.namprod<-as.character(c(substr(prod_nam,1,1),substr(prod_nam,1,4),prod_nam)) 
	  # character vector containing the first character of the COMPNAME of the product 
	  # m/z feature, the first four characters of its COMPNAME and the COMPNAME itself
	  if (length(which(Ent.namprod%in%No.name))!=0) {
	   if (length(which(Ent.namsub%in%No.name2))==0) {
	    new_nam<-paste("CSPP",cspp[[1]][1],cspp[[1]][2],sub_nam,sep=" ")
	    msdet[prod_compid,2]<-new_nam
	   }
  	  }
	 }
	 if (i==dim(msdet)[1]) break
	 i=i+1
	}
	i=1
	# When CSPP 'product' m/z feature is known
	repeat{
	 No.name=as.character(c("NULL","!"))
	 No.name2=as.character(c("NULL","!","CSPPr"))
	 if (msdet[i,3]!="NULL") {
	  sub_nam<-msdet[i,2]
	  Ent.namsub<-as.character(c(substr(sub_nam,1,1),substr(sub_nam,1,5),sub_nam))
	  cspp<-strsplit(msdet[i,3]," ")
	  prod_compid<-as.integer(cspp[[1]][3])
	  prod_nam<-msdet[prod_compid,2]
	  Ent.namprod<-as.character(c(substr(prod_nam,1,1),prod_nam))
	  if (length(which(Ent.namprod%in%No.name))==0) {
	   if (length(which(Ent.namsub%in%No.name2))!=0) {
	    new_nam<-paste("CSPPr",cspp[[1]][1],cspp[[1]][2],prod_nam,sep=" ")
	    msdet[i,2]<-new_nam
	   }
  	  }
	 }
	 if (i==dim(msdet)[1]) break
	 i=i+1
	}
	return(msdet)
}
