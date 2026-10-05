#' @title Fit a piecewise retention-time regression model
#'
#' @description
#' Helper function that models the relationship between retention times from
#' two experiments using linear and piecewise linear regression. Models with
#' zero, one, or two knot points are evaluated, and the model with the highest
#' adjusted R-squared is selected. The selected model is then refitted using
#' robust linear regression with `MASS::rlm()`.
#'
#' @param LCal.FT A numeric vector containing retention times from the first
#'   experiment.
#' @param LCal.Syn A numeric vector containing the corresponding retention
#'   times from the second experiment.
#' @param startpoint Numeric. Starting retention-time value used when
#'   evaluating candidate positions for the first knot of the piecewise
#'   regression model.
#'
#' @return A numeric vector containing six parameters of the selected
#'   retention-time regression model: the intercept, initial slope, first
#'   knot position, slope change after the first knot, second knot position,
#'   and slope change after the second knot.
#'
#'
#' @noRd
Regression_LCalign <- function(LCal.FT, LCal.Syn, startpoint) {
	knottime1<-c()
	knottime2<-c()
	VarExp<-c()
	y=LCal.Syn
	x=LCal.FT
	out<-lm(y~x)
	knottime1<-append(knottime1,0)
	knottime2<-append(knottime2,0)
	VarExp<-append(VarExp,summary(out)$adj.r.squared)
	i=startpoint
	while(i<as.integer(max(x))){
	 t1<-ifelse(x>i,x-i,0)
	 out<-lm(y~x+t1)
	 knottime1<-append(knottime1,i)
	 knottime2<-append(knottime2,0)
	 VarExp<-append(VarExp,summary(out)$adj.r.squared)
	 i=i+1
	}
	i=2
	while(i<as.integer(max(x))){
	 t1<-ifelse(x>i,x-i,0)
	 j=i+1
	 while(j<as.integer(max(x))){
	  t2<-ifelse(x>j,x-j,0)
	  out<-lm(y~x+t1+t2)
	  knottime1<-append(knottime1,i)
	  knottime2<-append(knottime2,j)
	  VarExp<-append(VarExp,summary(out)$adj.r.squared)
	  j=j+1
	 }
	 i=i+1
	}
	best.model<-data.frame(knottime1,knottime2,VarExp)
	best.model<-best.model[order(best.model[,3],decreasing=T),]
	t1<-ifelse(x>best.model[1,1],x-best.model[1,1],0)
	t2<-ifelse(x>best.model[1,2],x-best.model[1,2],0)
	library(MASS)
	if((best.model[1,1]==0)&(best.model[1,2]==0)){
	 exp.t1=0
	 exp.t2=0
	 out<-rlm(y~x)
	 exp.intercept<-summary(out)$coef[1,1]
	 exp.slope<-summary(out)$coef[2,1]
	}else{
	 if(best.model[1,2]==0){
	  exp.t2=0
	  out<-rlm(y~x+t1)
	  exp.intercept<-summary(out)$coef[1,1]
	  exp.slope<-summary(out)$coef[2,1]
	  exp.t1<-summary(out)$coef[3,1]
	 }else{
	  out<-rlm(y~x+t1+t2)
	  exp.intercept<-summary(out)$coef[1,1]
	  exp.slope<-summary(out)$coef[2,1]
	  exp.t1<-summary(out)$coef[3,1]
	  exp.t2<-summary(out)$coef[4,1]
	 }
	}
	rg<-c(exp.intercept,exp.slope,best.model[1,1],exp.t1,best.model[1,2],exp.t2)
	return(rg)
}
