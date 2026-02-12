#' @title Convert and rank GNPS-like MS/MS relationships and store results
#'  in SQLite
#'
#' @description
#' 
#' This function iterates the GNPS-like feature pair generation procedure over
#' all compound_id belonging to an experiment. the features (cfr. compound_id)   
#' are first sorted with increasing retention times. 
#' For each feature, CID spectral similarities are computed with all features 
#' that elute later. Only those feature pairs for which the CID spectral
#' similarity passes certain thresholds (thr1, thr2 and thr3) are retained. 
#' Subsequently, for a particular feature, all generated feature pairs are
#' ranked, via the rank.GNPS() function.
#'
#' @param sql_path `Character(1)` string giving the path to the SQLite database.
#' 
#' @param expid Integer. Experiment ID used to select compounds from the database.
#' 
#' @param peakwidth `Numeric(1)` half of the retention time window that is 
#' defined for the CSPP ‘substrate’ feature allowing
#' a minimum retention time difference with the CSPP ‘product’ feature,
#' by default set on \code{0.2}. 
#' 
#' @param mzerr `Numeric(1)` mass tolerance (in Da) used to match precursor and
#'   product ions. Default it is configured for ftms data
#'   \code{0.01}, for qtof data the user could set it to \code{0.015}.
#' 
#' @param min minimum absolute intensity of product ion to be withheld for CID 
#' spectral matching, by default \code{5} for QTOF, but the user can set it to
#'  \code{100} for FTMS spectra.
#'  
#' @param adduct mass of buffer, by default \code{46.0055 Da} representing 
#'  formic acid. 
#'  
#' @param thr1 average of the number of common ions and the number of common 
#' losses, set on \code{3} by default.
#' 
#' @param thr2 ratio of thr1 to the number of product ions in the CID spectrum 
#' of the feature in the feature pair having the lowest number of product ions,
#' by default set on \code{0.1}.
#' 
#' @param thr3 average of the dot products obtained for the common product ions 
#' and the common neutral losses, by default set on \code{0.4}.
#' 
#' @return
#' Invisibly returns a data frame corresponding to the updated GNPS SQL table
#' for the processed experiment.
#'
#' @author Ahlam Mentag
#' 
#' @export
conv.GNPS_SQL <- function(sql_path, expid, peakwidth = NULL, mzerr = 0.015,
                          min = 5, adduct = 46.0055, thr1 = 3, thr2 = 0.1,
                          thr3 = 0.4) {
  
  if (is.null(peakwidth)) peakwidth <- 0.2
  
  con <- dbConnect(RSQLite::SQLite(), sql_path)
  on.exit(dbDisconnect(con), add = TRUE)
  
  ## Load compounds
  inp.x <- dbGetQuery(
    con,
    sprintf(
      "SELECT compound_id,
              mass_measured,
              retention_time
       FROM ms_compound
       WHERE expid = %d
       ORDER BY retention_time",
      expid
    )
  )
  
  if (nrow(inp.x) == 0)
    stop("No compounds found for expid = ", expid)
  
  ## Load GNPS SQL table
  gnps_sql <- dbGetQuery(con, "SELECT * FROM gnps_add")
  
  i <- 1
  while (i <= nrow(inp.x)) {
    
    gnps.df <- data.frame(
      COMPID.sub   = integer(),
      MZ.sub       = numeric(),
      IONS.sub     = integer(),
      COMPID.prod  = integer(),
      MZ.prod      = numeric(),
      IONS.prod    = integer(),
      COMMON_IONS  = integer(),
      DOT_IONS     = numeric(),
      COMMON_LOSS  = integer(),
      DOT_LOSS     = numeric(),
      FORW_IONS    = numeric(),
      REV_IONS     = numeric(),
      FORW_LOSS    = numeric(),
      REV_LOSS     = numeric(),
      stringsAsFactors = FALSE
    )
    
    sub.mz <- inp.x$mass_measured[i]
    sub.rt <- inp.x$retention_time[i] + peakwidth
    
    sub.lowmz  <- sub.mz - mzerr
    sub.highmz <- sub.mz + 1.0034 + mzerr
    
    forb.lowmz  <- sub.mz + adduct - mzerr
    forb.highmz <- sub.mz + adduct + 1.0034 + mzerr
    
    j <- i + 1
    while (j <= nrow(inp.x)) {
      
      if (inp.x$retention_time[j] > sub.rt) {
        
        mzj <- inp.x$mass_measured[j]
        
        if ((mzj < sub.lowmz || mzj > sub.highmz) &&
            (mzj < forb.lowmz || mzj > forb.highmz)) {
          
          out <- targMS2comp_SQL(
            inp.x$compound_id[i],
            inp.x$compound_id[j],
            con,
            IntThres = min
          )
          
          thresh1 <- (out[1, 7] + out[1, 9]) / 2
          thresh2 <- thresh1 / min(out[1, c(3, 6)])
          thresh3 <- (out[1, 8] + out[1, 10]) / 2
          
          if (!is.na(thresh1) &&
              !is.na(thresh2) &&
              !is.na(thresh3) &&
              thresh1 > thr1 &&
              thresh2 > thr2 &&
              thresh3 > thr3) {
            
            gnps.df <- rbind(gnps.df, out)
          }
        }
      }
      j <- j + 1
    }
    
    if (nrow(gnps.df) > 0) {
      gnps_sql <- rank.GNPS_SQL(
        gnps.df,
        gnps_sql,
        inp.x
      )
    }
    
    i <- i + 1
  }
  
  ## Align compound_id
  gnps_sql$compound_id <- inp.x$compound_id[seq_len(nrow(gnps_sql))]
  
  ## Write back to SQL
  cols_to_update <- setdiff(colnames(gnps_sql), "compound_id")
  
  dbBegin(con)
  
  for (i in seq_len(nrow(gnps_sql))) {
    cid <- gnps_sql$compound_id[i]
    for (col in cols_to_update) {
      dbExecute(
        con,
        sprintf(
          "UPDATE gnps_add
           SET %s = ?
           WHERE compound_id = ?",
          col
        ),
        params = list(gnps_sql[i, col], cid)
      )
    }
  }
  
  dbCommit(con)
  
  invisible(gnps_sql)
}
