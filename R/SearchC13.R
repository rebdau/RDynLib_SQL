#' @title Annotate C12 and C13 isotope features
#'
#' @description
#' Helper function that evaluates candidate C12-C13 feature pairs and assigns
#' isotope annotations based on their relative abundances. For each candidate
#' pair, the feature with the lower abundance is annotated as `"C13"` and the
#' other as `"C12"`. A candidate C12 feature is not reassigned as `"C12"` if
#' it has already been annotated as `"C13"` in a previous comparison.
#'
#' @param mtr A data frame containing candidate isotope matches. Entries with
#'   `z3 == 1.0` are evaluated, while the `cw` and `rw` columns contain row
#'   indices referring to the candidate C12 and C13 features in `dt.s`.
#' @param dt.s A matrix or data frame containing feature information and
#'   abundance values used to evaluate candidate isotope pairs.
#' @param dt A matrix or data frame containing feature information. The first
#'   column contains the feature identifiers used to locate entries in
#'   `anno`.
#' @param s1 Integer. Number of candidate C12-C13 pairs to process.
#' @param anno A vector containing the current feature annotations.
#' @param findmn Integer or column index identifying the column in `dt.s`
#'   containing the abundance values used to compare candidate C12 and C13
#'   features.
#'
#' @return The updated `anno` vector with candidate isotope features annotated
#'   as `"C12"` or `"C13"`.
#'
#' @noRd
SearchC13 <- function(mtr, dt.s, dt, s1, anno, findmn) {
	t=0
	anno.t<-as.character(c())
	repeat{
		t=t+1
		nms<-dt.s[mtr[mtr$z3==1.0,]$cw[t],1]		
			# 'substrate' name, i.e. peak name of the candidate C12 (A)
		mngs<-dt.s[mtr[mtr$z3==1.0,]$cw[t],findmn]		
			# amount of candidate C12 peak				(A)
		nmp<-dt.s[mtr[mtr$z3==1.0,]$rw[t],1]			
			# 'product' name, i.e. peak name of candidate C13	(A)
		mngp<-dt.s[mtr[mtr$z3==1.0,]$rw[t],findmn]		
			# amount of candidate C13 peak				(A)
		if (mngs>mngp) {						
			# candidate C13 is only annotated as 'C13' if its amounts 
			#are lower than that of the candidate C12, otherwise it is 
			#annotated as 'C12'
			anno[which(dt[,1]==nmp)]="C13"
		}else{
			anno[which(dt[,1]==nmp)]="C12"
		}	
		anno.t<-anno[which(dt[,1]==nms)]
			# candidate C12 wordt geannoteerd als 'C12' enkel als er 
			#al geen annotatie als 'C13' eerder gebeurd was 
		if (anno.t=="C13"){
			if (t==s1) break
			next
		}else{
			anno[which(dt[,1]==nms)]="C12"
		}
		if (t==s1) break
	}
	return(anno)
}
