net.nodes_SQL <- function(sql_path, exp.id, min = 5) {
  
  library(DBI)
  library(RSQLite)
  
  con <- dbConnect(SQLite(), sql_path)
  on.exit(dbDisconnect(con), add = TRUE)
  
  query <- sprintf("
    SELECT c.compound_id,
           COUNT(p.Peak_id) AS nr_peaks
    FROM ms_compound c
    JOIN msms_spectrum s
         ON s.compound_id = c.compound_id
        AND s.ms_level = 2
    JOIN msms_spectrum_peak p
         ON p.Spectrum_id = s.spectrum_id
        AND p.Intensity >= %f
    WHERE c.expid = %d
    GROUP BY c.compound_id
    ORDER BY c.compound_id
  ", min, exp.id)
  
  nodes.net <- dbGetQuery(con, query)
  
  return(nodes.net)
}

