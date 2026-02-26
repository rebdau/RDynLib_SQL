#' @title Rank GNPS-like MS/MS relationships and update a GNPS annotation table.
#'
#' @description
#'  This function ranks MS/MS relationships produced by GNPS-style spectral
#'  comparisons, based (i) on the average of the numbers
#'  of common product ions and common neutral losses and (ii) on the average 
#'  of the dot products obtained for the common product ions and for the common.
#'
#' @param gnps.df Data frame containing GNPS-style MS/MS comparison results.
#'   Expected to include columns for compound IDs, ion counts, dot products,
#'   and forward/reverse similarity metrics.
#'   
#' @param gnps_add Data frame representing the GNPS annotation table to be
#'   updated. Must contain a \code{compound_id} column and sufficient columns
#'   to store ranked GNPS annotations.
#'   
#' @param inp.x Data frame of compounds containing at least \code{compound_id}
#'   and \code{mass_measured}.
#'
#' @return
#' A data frame corresponding to the updated GNPS annotation table.
#'
rank.GNPS_SQL <- function(gnps.df, gnps_add, inp.x) {
  
  # Scoring
  gnps.df$scor <- gnps.df$COMMON_IONS * gnps.df$DOT_IONS +
    gnps.df$COMMON_LOSS * gnps.df$DOT_LOSS
  
  gnps.df <- gnps.df[order(gnps.df$scor, decreasing = TRUE), ]
  
  gnps.df$ave_cnt <- rowMeans(cbind(gnps.df$FORW_IONS,
                                    gnps.df$REV_IONS), na.rm = TRUE)
  
  gnps.df$ave_dot <- rowMeans(cbind(gnps.df$DOT_IONS,
                                    gnps.df$DOT_LOSS), na.rm = TRUE)
  
  gnps.df$massdiff <- abs(
    inp.x$mass_measured[match(gnps.df$COMPID.prod, inp.x$compound_id)] -
      inp.x$mass_measured[match(gnps.df$COMPID.sub, inp.x$compound_id)]
  )
  
  for (k in seq_len(min(5, nrow(gnps.df)))) {
    
    compid_sub  <- gnps.df$COMPID.sub[k]
    compid_prod <- gnps.df$COMPID.prod[k]
    
    gnps_name <- paste(
      "GNPS","!!",
      gnps.df$COMMON_LOSS[k],"!",
      gnps.df$FORW_IONS[k],"!",
      gnps.df$REV_IONS[k],"!",
      compid_prod,"!!",
      round(gnps.df$massdiff[k],3),
      sep=""
    )
    
    row_idx <- match(compid_sub, gnps_add$compound_id)
    
    if (!is.na(row_idx))
      gnps_add[row_idx, k + 1] <- gnps_name
  }
  
  return(gnps_add)
}