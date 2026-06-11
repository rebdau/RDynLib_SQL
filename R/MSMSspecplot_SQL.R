#' @title plot MS2 spectra 
#'
#' @import DBI
#' @import RSQLite
#' 
#' @author Ahlam Mentag
#' 
#' @export

MSMSspecplot_SQL <- function(sql_path,
                             spectrum_id,
                             prdion,
                             neutloss,
                             err = NULL,
                             minum = NULL,
                             nr_col = NULL,
                             nr_col2 = NULL) {
  
  # Defaults
  if (is.null(err)) err <- 0.015
  if (is.null(minum)) minum <- 2
  if (is.null(nr_col)) nr_col <- 35
  if (is.null(nr_col2)) nr_col2 <- 6
  
  # Connect to SQLite
  con <- DBI::dbConnect(RSQLite::SQLite(), sql_path)
  on.exit(DBI::dbDisconnect(con), add = TRUE)
  
  # Get compound_id corresponding to spectrum_id
  dbKey <- DBI::dbGetQuery(
    con,
    sprintf(
      "SELECT compound_id
       FROM msms_spectrum
       WHERE spectrum_id = %d",
      spectrum_id
    )
  )
  
  if (nrow(dbKey) == 0)
    stop("Spectrum not found in msms_spectrum")
  
  compound_id <- dbKey$compound_id[1]
  
  # Retrieve compound information
  compound <- DBI::dbGetQuery(
    con,
    sprintf(
      "SELECT compound_id, name, retention_time, expid
       FROM ms_compound
       WHERE compound_id = %s",
      compound_id
    )
  )
  
  if (nrow(compound) == 0)
    stop(
      sprintf(
        "Compound %s not found in ms_compound",
        compound_id
      )
    )
  
  cat(
    "\nCompound:", compound$name,
    "\nRetention time:", compound$retention_time,
    "\nExperiment ID:", compound$expid,
    "\n"
  )
  
  # Use compound_id as dbkey
  dbkey <- compound_id
  
  # CSPP display
  cspp.res <- cspp.display_SQL(
    sql_path = sql_path,
    dbkey = dbkey,
    nr_col = nr_col
  )
  
  cat("\nAssociated CSPPs:\n")
  print(cspp.res)
  
  # GNPS display
  gnps.res <- gnps.display_SQL(
    sql_path = sql_path,
    dbkey = dbkey,
    nr_col2 = nr_col2
  )
  
  cat("\nAssociated GNPS conversions:\n")
  print(gnps.res)
  
  # Plot MS/MS spectrum
  oldpar <- par(no.readonly = TRUE)
  on.exit(par(oldpar), add = TRUE)
  
  MSMSplot_SQL(
    sql_path = sql_path,
    dbkey = dbkey,
    spectrum_id = spectrum_id,
    prdion = prdion,
    neutloss = neutloss,
    err = err,
    minum = minum,
    oldpar = oldpar
  )
  
  invisible(
    list(
      compound = compound,
      cspp = cspp.res,
      gnps = gnps.res
    )
  )
}