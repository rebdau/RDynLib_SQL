MSnplot_SQL <- function(sql_path, dbkey, prcx = 0.6) {
  
  library(DBI)
  library(RSQLite)
  
  con <- dbConnect(SQLite(), sql_path)
  on.exit(dbDisconnect(con), add = TRUE)
  
  ms3count <- dbGetQuery(con, sprintf(
    "SELECT COUNT(*) as n
     FROM msms_spectrum
     WHERE compound_id = %d
       AND ms_level = 3",
    dbkey
  ))$n
  
  totalplots <- 1 + ms3count
  
  row.nr <- floor(sqrt(totalplots))
  col.nr <- ceiling(totalplots / row.nr)
  
  oldpar <- par(no.readonly = TRUE)
  on.exit(par(oldpar), add = TRUE)
  
  par(mfrow = c(row.nr, col.nr))
  
  MS2plot_SQL(sql_path, dbkey, prcx)
  
  if (ms3count > 0) {
    for (i in 1:ms3count) {
      MS3plot_SQL(sql_path, dbkey, prcx, wh = i)
    }
  }
}
