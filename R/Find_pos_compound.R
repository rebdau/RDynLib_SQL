#' @title Find matching compounds based on mass and retention time
#'
#' @description
#' Helper function that searches for candidate compounds matching a reference
#' compound based on mass and retention-time criteria. The expected target
#' mass is calculated by adding 2.01456 Da to the reference mass, while the
#' expected retention time is calculated using a linear regression model.
#' Candidates falling within the specified mass and retention-time tolerances
#' are retained.
#'
#' @param mass.ftng Numeric. Measured mass of the reference compound.
#' @param time.ftng Numeric. Retention time of the reference compound.
#' @param ftps.o A matrix or data frame containing candidate target compounds,
#'   ordered by measured mass. The first and second columns contain retention
#'   time and measured mass, respectively.
#' @param err Numeric. Mass tolerance used to define the matching mass window.
#' @param lc.err Numeric. Retention-time tolerance used to define the matching
#'   retention-time window.
#' @param rg Numeric vector containing the intercept and slope of the
#'   retention-time regression model, respectively.
#'
#' @return A matrix containing the candidate compounds that satisfy both the
#'   mass and retention-time matching criteria. If no candidates are found,
#'   an empty matrix is returned.
#'
#'
#' @noRd
Find_pos_compound <- function(mass.ftng, time.ftng, ftps.o,
                              err, lc.err, rg) {
	mass.sel<-mass.ftng+2.01456
	mass.lb=mass.sel-err
	mass.ub=mass.sel+err
	pred_time=rg[1]+rg[2]*time.ftng
	time.lb=pred_time-lc.err
	time.ub=pred_time+lc.err
	pres<-array(dim=c(0,4))
	j=1
	while(j<=dim(ftps.o)[1]){
	 jak<-as.integer(rep(0,4))
	 jak[1]<-ifelse(ftps.o[j,2]>mass.lb,1,0)
	 jak[2]<-ifelse(ftps.o[j,2]<mass.ub,1,0)
	 jak[3]<-ifelse(ftps.o[j,1]>time.lb,1,0)
	 jak[4]<-ifelse(ftps.o[j,1]<time.ub,1,0)
	 if(all(jak==c(1,1,1,1))){
	  int<-ftps.o[j,1:4]
	  pres<-rbind(pres,int)
	 }
	 if(ftps.o[j,2]>=mass.ub)break
	 j=j+1
	}
	return(pres)
}
