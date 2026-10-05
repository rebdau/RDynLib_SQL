#' @title Annotate adduct-related features
#'
#' @description
#' Helper function that identifies features associated with a selected adduct
#' group and updates their annotations to `"adduct"`. Feature identifiers are
#' retrieved from `dt.s` using row indices stored in `mtr`, and the
#' corresponding entries in `anno.1` are updated.
#'
#' @param mtr A data frame containing adduct matching information. The `z3`
#'   column is used to select the adduct group and the `rw` column contains
#'   row indices referring to `dt.s`.
#' @param buf Value used to select the adduct group from the `z3` column of
#'   `mtr`.
#' @param dt.s A matrix or data frame containing feature information. The
#'   first column contains the feature identifiers used for annotation.
#' @param dt A matrix or data frame containing feature information. The first
#'   column contains the feature identifiers used to locate entries in
#'   `anno.1`.
#' @param s2 Integer. Number of matched entries to process.
#' @param anno.1 A vector containing the current feature annotations.
#'
#' @return The updated `anno.1` vector, with matched features annotated as
#'   `"adduct"`.
#'
#' @noRd
SearchAdduct <- function(mtr, buf, dt.s, dt, s2, anno.1) {
	t=0	
	repeat{
		t=t+1
		nmp<-dt.s[mtr[mtr$z3==buf,]$rw[t],1]	# (A)
		anno.1[which(dt[,1]==nmp)]="adduct"
		if (t==s2) break
	}
	return(anno.1)
}
