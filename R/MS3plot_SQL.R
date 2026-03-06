#' @title plot MS3 spectra 
#'
#' @import DBI
#' @import RSQLite
#' 
#' @author Ahlam Mentag
#' 
#' @export
MS3plot_SQL <- function(sql_path, dbkey, prcx = 0.7, wh = 1) {
  
  con <- dbConnect(SQLite(), sql_path)
  on.exit(dbDisconnect(con), add = TRUE)
  
  ms3 <- dbGetQuery(con, sprintf(
    "SELECT spectrum_id, precursor_mz
     FROM msms_spectrum
     WHERE compound_id = %d
       AND ms_level = 3
     ORDER BY spectrum_id",
    dbkey
  ))
  
  if (nrow(ms3) < wh) return()
  
  spectrum_id <- ms3$spectrum_id[wh]
  precursor   <- ms3$precursor_mz[wh]
  
  peaks <- dbGetQuery(con, sprintf(
    "SELECT mz, intensity
   FROM msms_spectrum_peak
   WHERE spectrum_id = %d
   ORDER BY mz",
    spectrum_id
  ))
  
  if (nrow(peaks) == 0 || any(is.na(peaks$mz))) return()
  
  prod_ion <- peaks$mz
  intens   <- peaks$intensity
  
  
  par(cex = prcx)
  
  plot(prod_ion,
       intens,
       type = "h",
       xlab = "m/z",
       ylab = "ion intensity",
       main = paste("MS3: precursor m/z", round(precursor, 2)),
       xlim = c(min(prod_ion) * 0.9, max(prod_ion) * 1.1),
       ylim = c(0, max(intens) * 1.1))
  
  text(prod_ion, intens, round(prod_ion, 2), pos = 3)
}
