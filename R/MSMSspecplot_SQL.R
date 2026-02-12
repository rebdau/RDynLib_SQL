MSMSspecplot_SQL <- function(sql_path,
                             dbkey,
                             err = NULL,
                             minum = NULL,
                             nr_col = NULL,
                             nr_col2 = NULL) {
  
  library(DBI)
  library(RSQLite)
  
  # Defaults
  if (is.null(err)) err <- 0.015
  if (is.null(minum)) minum <- 2
  if (is.null(nr_col)) nr_col <- 35
  if (is.null(nr_col2)) nr_col2 <- 6
  
  # Connect to SQLite
  con <- dbConnect(SQLite(), sql_path)
  on.exit(dbDisconnect(con), add = TRUE)
  
  # Check compound exists
  compound <- dbGetQuery(con, sprintf(
    "SELECT compound_id, name, retention_time, expid
     FROM ms_compound
     WHERE compound_id = %d",
    dbkey
  ))
  
  if (nrow(compound) == 0) {
    stop("Compound not found in ms_compound")
  }
  
  cat("\nCompound:", compound$name,
      "\nRetention time:", compound$retention_time,
      "\nExperiment ID:", compound$expid, "\n")
  
  # CSPP display
  cspp.res <- cspp.display_SQL(sql_path = sql_path,
                               dbkey = dbkey,
                               nr_col = nr_col)
  
  cat("\nAssociated CSPPs:\n")
  print(cspp.res)
  
  #  GNPS display 
  gnps.res <- gnps.display_SQL(sql_path = sql_path,
                               dbkey = dbkey,
                               nr_col2 = nr_col2)
  
  cat("\nAssociated GNPS conversions:\n")
  print(gnps.res)
  
  #Plot MS/MS spectrum 
  oldpar <- par(no.readonly = TRUE)
  
  MSMSplot_SQL(sql_path = sql_path,
               dbkey = dbkey,
               err = err,
               minum = minum,
               oldpar = oldpar,
               nl = "nl")
  
  par(oldpar)
  
  return(list(cspp.res, gnps.res))
}
