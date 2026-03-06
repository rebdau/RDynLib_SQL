#' @title plot MSnspectra 
#'
#' @import DBI
#' @import RSQLite
#' 
#' @author Ahlam Mentag
#' 
#' @export
MSnspecplot_SQL <- function(sql_path,
                            dbkey,
                            nr_col = 35,
                            nr_col2 = 6,
                            lc.err = 0.02,
                            mz.err = 0.001,
                            MS1 = NULL,
                            prcx = NULL) {
  
  if (is.null(prcx)) prcx <- 0.6
  
  con <- dbConnect(SQLite(), sql_path)
  on.exit(dbDisconnect(con), add = TRUE)
  
  # CSPP + GNPS
  cspp.res <- cspp.display_SQL(sql_path, dbkey, nr_col)
  cat("\nAssociated CSPPs:\n")
  print(cspp.res)
  
  gnps.res <- gnps.display_SQL(sql_path, dbkey, nr_col2)
  cat("\nAssociated GNPS conversions:\n")
  print(gnps.res)
  
  # MS1 
  if (!is.null(MS1)) {
    
    # count MS3 spectra
    ms3_count <- dbGetQuery(con, sprintf(
      "SELECT COUNT(*) AS n
       FROM msms_spectrum
       WHERE compound_id = %d
         AND ms_level = 3",
      dbkey
    ))$n
    
    msnle <- ifelse(ms3_count > 0, ms3_count + 2, 2)
    
    row.nr <- floor(sqrt(msnle))
    col.nr <- ceiling(msnle / row.nr)
    
    oldpar <- par(no.readonly = TRUE)
    on.exit(par(oldpar), add = TRUE)
    
    par(mfrow = c(row.nr, col.nr))
    
    single <- "single"
    
    MS1peaks <- MS1specplot_SQL(sql_path,
                                dbkey,
                                lc.err = lc.err,
                                mz.err = mz.err,
                                prcx = prcx)
    
    MSnplot_SQL(sql_path,
                dbkey,
                prcx = prcx,
                single = single)
    
    cat("\nMS1 peaks:\n")
    print(MS1peaks)
    
    return(list(MS1peaks = MS1peaks,
                cspp = cspp.res,
                gnps = gnps.res))
  }
  
  
  ms3_count <- dbGetQuery(con, sprintf(
    "SELECT COUNT(*) AS n
     FROM msms_spectrum
     WHERE compound_id = %d
       AND ms_level = 3",
    dbkey
  ))$n
  
  msnle <- 1 + ms3_count
  
  row.nr <- floor(sqrt(msnle))
  col.nr <- ceiling(msnle / row.nr)
  
  oldpar <- par(no.readonly = TRUE)
  on.exit(par(oldpar), add = TRUE)
  
  par(mfrow = c(row.nr, col.nr))
  
  single <- "single"
  
  # Plot MS2
  MSnplot_SQL(sql_path,
              dbkey,
              prcx = prcx,
              single = single)
  
  # Plot ALL MS3 
  if (ms3_count > 0) {
    for (i in seq_len(ms3_count)) {
      MS3specplot_SQL(sql_path = sql_path,
                      dbkey = dbkey,
                      prcx = prcx,
                      wh = i,
                      single = single)
    }
  }
  
  return(list(cspp = cspp.res,
              gnps = gnps.res))
}
