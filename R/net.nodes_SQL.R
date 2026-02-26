#' @title generating lists of network nodes and edges.
#'
#' @description
#' net.nodes_SQL() connects to a SQLite database and returns a table of 
#' compounds 'compound_id' for a given experiment 'exp.id', along with 
#' the number of MS2 peaks 'nr_peaks' per compound whose intensity is greater 
#' than or equal to a specified threshold 'min'.
#' 
#' @param sql_path 'character(1)' path to the sqlite database.
#' 
#' @param exp.id 'number(1)' experiment to use for network generation.
#' 
#' @param min 'numeric(1)' minimum intensity for the spectra.
#' 
#' @import DBI
#' @import RSQLite
#' 
#' @author Ahlam Mentag
#' 
#' @export
net.nodes_SQL <- function(sql_path, exp.id, min) {
  
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

