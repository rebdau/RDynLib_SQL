#' @title generating the nodes and edges for the experiment-based global network.
#'
#' @description 
#' The net.addit_SQL() function is responsible for generating lists of network  
#' nodes and edges by the net.nodes_SQL() and either the allCSPP_SQL() 
#' or allGNPS_SQL() functions, respectively. This generates the nodes and edges 
#' for the experiment-based global network. Only edges are retained that pass 
#' certain CID spectral similarity thresholds.
#' 
#' @param sql_path 'character(1)' path to the sqlite database.
#' 
#' @param exp.id 'number(1)' experiment to use for network generation.
#' 
#' @param nettype 'character(1)' the graph type the user want to display, 
#' it could be either 
#' 
#' @param nr_col number of columns in the compound_add and gnps_add tables.
#' 
#' @return list of edges and nodes.
#' 
#' @import DBI
#' @import RSQLite
#' @author Ahlam Mentag
#' 
#' @export
net.addit_SQL <- function(sql_path,
                          exp.id,
                          nettype,
                          nr_col = NULL,
                          min = NULL,
                          thr1 = 3,
                          thr2 = 0.1,
                          thr3 = 0.2) {
  
  con <- dbConnect(RSQLite::SQLite(), sql_path)
  on.exit(dbDisconnect(con), add = TRUE)
  
  compounds <- dbGetQuery(
    con,
    sprintf(
      "SELECT compound_id
       FROM ms_compound
       WHERE expid = %d
       ORDER BY compound_id",
      exp.id
    )
  )
  
  if (nrow(compounds) == 0)
    stop("No compounds found for expid = ", exp.id)
  
  if (is.null(min)) {
    min <- 5  
  }
  
  ## Nodes
  nodes.net <- net.nodes_SQL(
    sql_path = sql_path,
    exp.id = exp.id,
    min = min
  )
  
  ## If nr_col is NULL, detect automatically from database
  if (is.null(nr_col)) {
    
    if (nettype == "CSPP") {
      table_info <- dbGetQuery(con, "PRAGMA table_info(compound_add)")
    } else if (nettype == "GNPS") {
      table_info <- dbGetQuery(con, "PRAGMA table_info(gnps_add)")
    } else {
      stop("Type of net should be CSPP or GNPS.")
    }
    
    nr_col <- nrow(table_info)
  }
  
  ## Edges
  if (nettype == "CSPP") {
    
    edges.net <- allCSPP_SQL(
      sql_path = sql_path,
      exp.id = exp.id,
      nr_col = nr_col,
      thr1 = thr1,
      thr2 = thr2,
      thr3 = thr3
    )
    
  } else if (nettype == "GNPS") {
    
    edges.net <- allGNPS_SQL(
      sql_path = sql_path,
      exp.id = exp.id,
      nr_col = nr_col,
      thr1 = thr1,
      thr2 = thr2,
      thr3 = thr3
    )
  }
  
  ## Additional metrics
  common.nr <- edges.net[,4] * edges.net[,5]
  edges.net <- data.frame(edges.net, common.nr)
  
  common.nr.rel <- edges.net[,7] / max(edges.net[,7])
  edges.net <- data.frame(edges.net, common.nr.rel)
  
  dot.rel <- round(edges.net[,6] * 8, digits = 0)
  edges.net <- data.frame(edges.net, dot.rel)
  
  return(list(nodes.net, edges.net))
}