#' @title Plot FTMS-Synapt retention-time alignment
#'
#' @description
#' Helper function that visualizes the retention-time alignment between FTMS
#' and Synapt data using a previously fitted piecewise regression model.
#' FTMS and Synapt retention times are extracted from the alignment table and
#' displayed as a scatter plot. The regression curve defined by `rg` is added
#' to the plot, together with the corresponding regression equation.
#'
#' @param LCal A matrix or data frame containing the LC alignment results.
#'   The second column contains FTMS retention times and the fourth column
#'   contains Synapt retention times.
#' @param rg A numeric vector containing the parameters of the piecewise
#'   retention-time regression model.
#'
#' @return This function is called for its side effect of generating a
#'   retention-time alignment plot. It does not explicitly return a value.
#'
#'
#' @noRd
PlotPie_LCalign <- function(LCal, rg) {
	LCal.FT<-as.numeric(LCal[,2]) # FT retention times
	LCal.Syn<-as.numeric(LCal[,4]) # Synapt retention times
	plot(LCal.FT~LCal.Syn,xlab="Synapt",ylab="FTMS",mgp=c(2,0.5,0),col=5,cex.axis=0.8)
	min.FT=max(LCal.FT)/1000
	x.p=seq(min.FT,max(LCal.FT),min.FT)
	t1.p=as.numeric(rep(0,1000))
	t2.p=as.numeric(rep(0,1000))
	for (i in 1:1000){
	 if(rg[3]!=0){
	  t1.p[i]=ifelse(x.p[i]>rg[3],x.p[i]-rg[3],0)
	 }
	 if(rg[5]!=0){
	  t2.p[i]=ifelse(x.p[i]>rg[5],x.p[i]-rg[5],0)
	 }
	}
	y.p=rg[1]+rg[2]*x.p+rg[4]*t1.p+rg[6]*t2.p
	lines(y.p,x.p) 
	pst<-paste("IWLS: Syn=",round(rg[1],digits=2),"+",round(rg[2],digits=2),"FT+",
		round(rg[4],digits=2),"FT[time>",round(rg[3],digits=0),"]+",
		round(rg[6],digits=2),"FT[time>",round(rg[5],digits=0),"]")
	title(main="FT vs Synapt retention times (min)",sub=pst,cex.main=1,cex.sub=0.8)
}
