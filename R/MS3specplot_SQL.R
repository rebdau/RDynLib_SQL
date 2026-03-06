#' @title plot MS3 spectra 
#'
#' @import DBI
#' @import RSQLite
#' 
#' @author Ahlam Mentag
#' 
#' @export
MS3specplot_SQL <- function(sql_path,
                            dbkey,
                            prcx = NULL,
                            wh = NULL,
                            single = NULL) {


  if (is.null(prcx)) prcx <- 0.7
  if (is.null(wh))   wh   <- 1
  
  # Connect to SQLite
  con <- dbConnect(SQLite(), sql_path)
  on.exit(dbDisconnect(con), add = TRUE)
  
  # Check compound exists
  compound <- dbGetQuery(con, sprintf(
    "SELECT compound_id, name
     FROM ms_compound
     WHERE compound_id = %d",
    dbkey
  ))
  
  if (nrow(compound) == 0) {
    stop("Compound not found in ms_compound")
  }
  
  cat("\nCompound:", compound$name, "\n")
  
  # Layout handling
  if (is.null(single)) {
    oldpar <- par(no.readonly = TRUE)
    on.exit(par(oldpar), add = TRUE)
    par(mfrow = c(1, 1))
  }
  
  # Call MS3 plot
  MS3plot_SQL(sql_path = sql_path,
              dbkey = dbkey,
              prcx = prcx,
              wh = wh)
  
  invisible(NULL)
}
