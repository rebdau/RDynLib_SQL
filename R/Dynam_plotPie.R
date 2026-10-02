#' @title Add selected feature matches to a retention-time alignment plot
#'
#' @description
#' Helper function that generates a retention-time alignment plot and
#' highlights selected feature matches. Dashed lines indicate the retention
#' times of candidate FT and QTOF features, while the selected feature pair
#' is highlighted on the plot.
#'
#' @param LCal A data frame or matrix containing the retention-time alignment
#'   data used to generate the plot.
#' @param rg Numeric value passed to `PlotPie_LCalign()` to define the plotting
#'   range.
#' @param z Integer index identifying the selected feature in `pres1`.
#' @param v Integer index identifying the selected feature in `pres2`.
#' @param pres1 A data frame or matrix containing candidate FT features, with
#'   retention times in the first column.
#' @param pres2 A data frame or matrix containing candidate QTOF features, with
#'   retention times in the first column.
#'
#' @return This function is called for its side effect of generating and
#'   annotating a retention-time alignment plot. It does not explicitly
#'   return a value.
#'
#' @author Ahlam Mentag
#'
#' @noRd
Dynam_plotPie<-function(LCal,rg,z,v,pres1,pres2){
	PlotPie_LCalign(LCal,rg)
	iso.ft<-round(pres1[,1],digits=2)
	iso.syn<-round(pres2[,1],digits=2)
	for (i in 1:length(iso.ft)){
	 abline(h=iso.ft[i],lty=2,col=4)
	}
	for(i in 1:length(iso.syn)){
	 abline(v=iso.syn[i],lty=2,col=4)
	}
	points(pres2[v,1],pres1[z,1],col=2,cex=1.2,pch=21)
}