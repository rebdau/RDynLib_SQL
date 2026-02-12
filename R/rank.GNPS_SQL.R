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
#' @param Prod.dat Data frame of compound metadata.
#'   
#' @param inp.x Data frame of compounds containing at least \code{compound_id}
#'   and \code{mass_measured}.
#'
#' @return
#' A data frame corresponding to the updated GNPS annotation table.
#'
rank.GNPS_SQL <- function(gnps.df, gnps_add, Prod.dat, inp.x) {
  
  scor <- gnps.df[,7] * gnps.df[,8] +
    gnps.df[,9] * gnps.df[,10]
  
  gnps.df <- data.frame(gnps.df, scor)
  gnps.df <- gnps.df[order(gnps.df$scor, decreasing = TRUE), ]
  
  ave.cnt <- (gnps.df[,12] + gnps.df[,14]) / 2
  ave.dot <- (gnps.df[,8]  + gnps.df[,10]) / 2
  
  massdiff <- abs(
    inp.x$mass_measured[match(gnps.df[,4], inp.x$compound_id)] -
      inp.x$mass_measured[match(gnps.df[,1], inp.x$compound_id)]
  )
  
  gnps.df <- data.frame(gnps.df, ave.cnt, ave.dot, massdiff)
  
  if (nrow(gnps.df) > 0) {
    
    k <- 1
    while (k <= nrow(gnps.df) && k <= 5) {
      
      compid.sub  <- gnps.df[k, 1]
      compid.prod <- gnps.df[k, 4]
      
      gnps.name <- paste(
        "GNPS!!",
        gnps.df[k, 6], "!",
        gnps.df[k, 16], "!",
        gnps.df[k, 17], "!",
        compid.prod,
        "!!",
        round(gnps.df[k, 18], 3),
        sep = ""
      )
      
      ## SQL row selection
      row_idx <- which(gnps_add$compound_id == compid.sub)
      
      if (length(row_idx) == 1) {
        gnps_add[row_idx, k + 1] <- gnps.name
      }
      
      k <- k + 1
    }
  }
  
  gnps_add
}
