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
  
  compid.start <- min(compounds$compound_id)
  compid.end   <- max(compounds$compound_id)
  
  if (is.null(min)) {
    min <- 5   # default safe value (you can store analyzer in DB later)
  }
  
  ##Nodes 
  nodes.net <- net.nodes_SQL(
    sql_path = sql_path,
    exp.id = exp.id,
    min = min
  )
  
  ##Edges 
  if (nettype == "CSPP") {
    
    if (is.null(nr_col)) nr_col <- 35
    
    edges.net <- allCSPP_SQL(
      sql_path = sql_path,
      compid.start = compid.start,
      compid.end = compid.end,
      nr_col = nr_col,
      thr1 = thr1,
      thr2 = thr2,
      thr3 = thr3
    )
    
  } else if (nettype == "GNPS") {
    
    if (is.null(nr_col)) nr_col <- 6
    
    edges.net <- allGNPS_SQL(
      sql_path = sql_path,
      compid.start = compid.start,
      compid.end = compid.end,
      nr_col = nr_col,
      thr1 = thr1,
      thr2 = thr2,
      thr3 = thr3
    )
    
  } else {
    stop("Type of net should be CSPP or GNPS.")
  }
  
  ##Additional metrics
  common.nr <- edges.net[,4] * edges.net[,5]
  edges.net <- data.frame(edges.net, common.nr)
  
  common.nr.rel <- edges.net[,7] / max(edges.net[,7])
  edges.net <- data.frame(edges.net, common.nr.rel)
  
  dot.rel <- round(edges.net[,6] * 8, digits = 0)
  edges.net <- data.frame(edges.net, dot.rel)
  
  return(list(nodes.net, edges.net))
}
