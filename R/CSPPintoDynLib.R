#' @title Integrate CSPP results into a DynLib compound table
#'
#' @description
#' Helper function that filters CSPP matches based on spectral similarity
#' criteria and integrates the best retained conversion into the corresponding
#' compound entry of a DynLib database table.
#'
#' @param inp.x A data frame containing the DynLib compound information to be
#'   updated with CSPP conversion annotations.
#' @param comp_add A data frame containing CSPP results for the compounds and
#'   the different biotransformations or conversions evaluated.
#' @param thr1 Numeric. Minimum spectral similarity score required to retain
#'   a CSPP match.
#' @param thr2 Numeric. Minimum combined score, calculated from the number of
#'   matched ions, the proportion of shared ions, and the spectral similarity
#'   score, required to retain a CSPP match.
#'
#' @return The updated `inp.x` data frame, with the best retained CSPP
#'   conversion annotation added to the corresponding compound entries.
#'
#'
#' @noRd
CSPPintoDynLib<-function(inp.x,comp_add,thr1,thr2){
	biotr<-colnames(comp_add)
	i=1
	repeat{
	 compid_sub<-comp_add[i,]
	 cspp_nr<-c()
	 cspp_shared<-c()
	 cspp_dot<-c()
	 cspp_prod<-c()
	 cspp_conv<-c()
	 j=3
	 while (j<=length(compid_sub)){
	  # For a particular row (i.e. m/z feature) in compound_add.txt, go through all
	  # conversions, and check for which of them a CSPP was found, add them to cspp_res 
	  if(1%in%grep("!!",compid_sub[j])){  
	   cspp_inf<-strsplit(as.character(compid_sub[j]),"!")
	   cspp_conv<-append(cspp_conv,biotr[j])
	   cspp_nr<-append(cspp_nr,as.integer(cspp_inf[[1]][3]))
	   cspp_shared<-append(cspp_shared,as.numeric(cspp_inf[[1]][4]))
	   cspp_dot<-append(cspp_dot,as.numeric(cspp_inf[[1]][5]))
	   cspp_prod<-append(cspp_prod,as.integer(cspp_inf[[1]][7]))
	  }
	  j=j+1
	 }
	 cspp_res<-data.frame(cspp_conv,cspp_prod,cspp_nr,cspp_shared,cspp_dot)
	 if(nrow(cspp_res)!=0){
	  # retain from the cspp_res only those CSPPs for which the CID spectral
	  # similarity match surpasses a certain threshold and order the resulting
	  # cspp_ka output. 
	  cspp_ka<-cspp_res[(cspp_res[,5]>thr1)&(cspp_res[,3]*cspp_res[,4]*cspp_res[,5]>thr2),]
	  if(nrow(cspp_ka)!=0){
	   cspp_ka<-cspp_ka[order(cspp_ka[,5],decreasing=T),]
	   # the cspp with the highest CID spectral similarity (first entry in cspp_ka)
	   # is used to enter the COMPID of the 'product' m/z feature into the 'substrate'
	   # m/z feature cell of the CONVERSION column of compound.csv 
	   if(!is.na(cspp_ka[1,5])){
	    conv_name=paste(as.character(cspp_ka[1,1]),"KA",as.character(cspp_ka[1,2]),collapse="$")
		inp.x[i,11]<-conv_name
	   }
	  }
	 }
	 if (i>dim(comp_add)[1]) break
	 i=i+1
	}
	return(inp.x)
}