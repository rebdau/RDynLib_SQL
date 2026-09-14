available_spectrum_types <- function(con) {
  
  dbGetQuery(
    con,
    "SELECT DISTINCT spectrum_type
     FROM msms_spectrum
     ORDER BY spectrum_type"
  )
}