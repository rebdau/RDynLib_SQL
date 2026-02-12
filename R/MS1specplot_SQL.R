MS1specplot_SQL <- function(sql_path,
                            dbkey,
                            lc.err = 0.02,
                            mz.err = 0.001,
                            prcx = 0.7) {
  
  library(DBI)
  library(RSQLite)
  
  con <- dbConnect(SQLite(), sql_path)
  on.exit(dbDisconnect(con), add = TRUE)
  
  compound <- dbGetQuery(con, sprintf(
    "SELECT retention_time, mass_measured
     FROM ms_compound
     WHERE compound_id = %d",
    dbkey
  ))
  
  if (nrow(compound) == 0) return()
  
  rt0 <- compound$retention_time
  mz0 <- compound$mass_measured
  
  # Get MS1 spectrum near RT
  ms1_spec <- dbGetQuery(con, sprintf(
    "SELECT spectrum_id
     FROM msms_spectrum
     WHERE ms_level = 1
       AND rtime BETWEEN %f AND %f
     LIMIT 1",
    rt0 - lc.err, rt0 + lc.err
  ))
  
  if (nrow(ms1_spec) == 0) return()
  
  spectrum_id <- ms1_spec$spectrum_id[1]
  
  peaks <- dbGetQuery(con, sprintf(
    "SELECT Mz, Intensity
     FROM msms_spectrum_peak
     WHERE Spectrum_id = %d
     ORDER BY Mz",
    spectrum_id
  ))
  
  if (nrow(peaks) == 0) return()
  
  relabun <- round(100 * peaks$Intensity / max(peaks$Intensity), 1)
  mzdiff  <- round(peaks$Mz - mz0, 4)
  
  par(cex = prcx)
  
  plot(peaks$Mz,
       relabun,
       type = "h",
       xlab = "m/z",
       ylab = "relative intensity",
       main = paste("MS1: m/z", round(mz0, 4)),
       xlim = c(min(peaks$Mz) * 0.9, max(peaks$Mz) * 1.1),
       ylim = c(0, max(relabun) * 1.1))
  
  text(peaks$Mz, relabun, round(peaks$Mz, 3), pos = 3)
  text(peaks$Mz, relabun, round(mzdiff, 3), pos = 3,
       offset = 1.5, col = 2)
  
  return(data.frame(
    mz = peaks$Mz,
    mz_diff = mzdiff,
    rel_int = relabun
  ))
}
