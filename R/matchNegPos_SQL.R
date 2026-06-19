#' @title Match  MS2 spectra between the same instruments type 
#'        using SQL databases. 
#'
#' @description
#' 
#' `matchNegPos_SQL()` Removes false positive negFTm/z,posFTm/z peak pairs by 
#                      checking whether a certain percentage of the neg MS2  
#                      spectral peaks can be traced in the pos MS2 spectrum. 
#                      By default, 20% of the neg MS2 peaks are traced in the
#                      pos MS2 spectrum.
#'
#' @param LCal 
#' `data.frame` returned by `Aligning_General_SQL()` function.
#' 
#' @param FTn_con `DBIConnection`  
#'   A `DBI connection` object to the reference SQLite database.
#'
#' @param FTp_con `DBIConnection`  
#'   A `DBI connection` object to the SQLite database to align with.
#'
#' @param polarity_ftn 
#'        `Integer (0/1)` Polarity of negative FTMS or QTOF experiment.
#'        
#' @param polarity_ftp 
#'        `Integer (0/1)` Polarity of positive FTMS or QTOF experiment.
#' 
#' @param minPeaks `Numeric(1)` Minimum matching ratio (default 0.2).
#' 
#'
#' @return Filtered `LCal` with MS2-supported neg–pos pairs.
#'
#' @import DBI RSQLite
#' @export
matchNegPos_SQL <- function(
    LCal,
    FTn_con,
    FTp_con,
    polarity_ftn = 0,
    polarity_ftp = 1,
    minIon = 0.2,
    tol = 0.01
) {
  
  
  
  ms2_neg_df <- dbGetQuery(FTn_con, sprintf("
    SELECT s.compound_id, p.mz
    FROM msms_spectrum s
    JOIN msms_spectrum_peak p USING(spectrum_id)
    WHERE s.ms_level = 2
      AND s.polarity = %d
  ", polarity_ftn))
  
  
  
  ms2_pos_df <- dbGetQuery(FTp_con, sprintf("
    SELECT s.compound_id, p.mz
    FROM msms_spectrum s
    JOIN msms_spectrum_peak p USING(spectrum_id)
    WHERE s.ms_level = 2
      AND s.polarity = %d
  ", polarity_ftp))
  
  
  
  ms2_neg <- split(ms2_neg_df$mz, ms2_neg_df$compound_id)
  ms2_pos <- split(ms2_pos_df$mz, ms2_pos_df$compound_id)
  
  
  
  if (is.null(LCal) || nrow(LCal) == 0) return(LCal)
  
  i <- 1
  
  
  
  while (i <= nrow(LCal)) {
    
    neg_id <- LCal[i, 1]
    pos_id <- LCal[i, 7]
    
    neg_peaks <- ms2_neg[[as.character(neg_id)]]
    pos_peaks <- ms2_pos[[as.character(pos_id)]]
    
    # remove empty spectra
    if (is.null(neg_peaks) || is.null(pos_peaks)) {
      LCal <- LCal[-i, ]
      next
    }
    
    
     pos_peaks <- pos_peaks - 2
    
    
    
    same_ion <- sum(sapply(neg_peaks, function(x) {
      any(abs(pos_peaks - x) <= tol)
    }))
    
    ratio <- same_ion / length(neg_peaks)
    
    
    
    if (ratio < minIon) {
      LCal <- LCal[-i, ]
      next
    }
    
    i <- i + 1
  }
  
  unique(LCal)
}