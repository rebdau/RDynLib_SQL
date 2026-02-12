MSnspecplot_SQL <- function(sql_path,
                            dbkey,
                            nr_col = 35,
                            nr_col2 = 6,
                            lc.err = 0.02,
                            mz.err = 0.001,
                            MS1 = NULL,
                            prcx = 0.6) {
  
  library(DBI)
  library(RSQLite)
  
  con <- dbConnect(SQLite(), sql_path)
  on.exit(dbDisconnect(con), add = TRUE)
  
  
  # CSPP relationships
  cspp.res <- cspp.display_SQL(
    sql_path = sql_path,
    dbkey = dbkey,
    nr_col = nr_col
  )
  
  cat("\nAssociated CSPPs:\n")
  print(cspp.res)
  
  # GNPS relationships
  gnps.res <- gnps.display_SQL(
    sql_path = sql_path,
    dbkey = dbkey,
    nr_col = nr_col2
  )
  
  cat("\nAssociated GNPS conversions:\n")
  print(gnps.res)
  
  #  MS1 + MSn plotting
  if (!is.null(MS1)) {
    
    oldpar <- par(no.readonly = TRUE)
    on.exit(par(oldpar), add = TRUE)
    
    par(mfrow = c(1, 2))
    
    MS1peaks <- MS1specplot_SQL(
      sql_path = sql_path,
      dbkey = dbkey,
      lc.err = lc.err,
      mz.err = mz.err,
      prcx = prcx
    )
    
    MSnplot_SQL(
      sql_path = sql_path,
      dbkey = dbkey,
      prcx = prcx
    )
    
    cat("\nMS1 peaks:\n")
    print(MS1peaks)
    
    return(list(MS1peaks = MS1peaks,
                cspp = cspp.res,
                gnps = gnps.res))
    
  } else {
    
    MSnplot_SQL(
      sql_path = sql_path,
      dbkey = dbkey,
      prcx = prcx
    )
    
    return(list(cspp = cspp.res,
                gnps = gnps.res))
  }
}
