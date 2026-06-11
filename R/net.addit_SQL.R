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
#' @param min 'numeric(1)' minimum intensity for the spectra.
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
  
  con <- DBI::dbConnect(RSQLite::SQLite(), sql_path)
  on.exit(DBI::dbDisconnect(con), add = TRUE)
  

  ## CHECK COMPOUNDS

  
  compounds <- DBI::dbGetQuery(
    con,
    sprintf(
      "SELECT compound_id
       FROM ms_compound
       WHERE expid = %d
       ORDER BY compound_id",
      exp.id
    )
  )
  
  if (nrow(compounds) == 0) {
    stop("No compounds found for expid = ", exp.id)
  }
  

  ## NODES

  
  nodes.net <- net.nodes_SQL(
    sql_path = sql_path,
    exp.id = exp.id,
    min = min
  )
  

  ## DETECT NUMBER OF COLUMNS

  
  if (is.null(nr_col)) {
    
    if (nettype == "CSPP") {
      
      table_info <- DBI::dbGetQuery(
        con,
        "PRAGMA table_info(compound_add)"
      )
      
    } else if (nettype == "GNPS") {
      
      table_info <- DBI::dbGetQuery(
        con,
        "PRAGMA table_info(gnps_add)"
      )
      
    } else {
      
      stop("Type of net should be CSPP or GNPS.")
    }
    
    nr_col <- nrow(table_info)
  }
  

  ## BUILD EDGES

  
  if (nettype == "CSPP") {
    
    edges.net <- allCSPP_SQL(
      sql_path = sql_path,
      exp.id = exp.id,
      nr_col = nr_col,
      thr1 = thr1,
      thr2 = thr2,
      thr3 = thr3
    )
    
  } else {
    
    edges.net <- allGNPS_SQL(
      sql_path = sql_path,
      exp.id = exp.id,
      nr_col = nr_col,
      thr1 = thr1,
      thr2 = thr2,
      thr3 = thr3
    )
  }
  
  edges.net <- as.data.frame(edges.net)
  

  ## HANDLE EMPTY NETWORKS

  
  if (nrow(edges.net) == 0 || ncol(edges.net) < 6) {
    
    message(nettype, " network is empty.")
    
    empty.df <- data.frame(
      compid.sub = character(),
      compid.prod = character(),
      common.nr = numeric(),
      common.nr.rel = numeric(),
      dot.rel = numeric(),
      ave.dot = numeric(),
      conv.type = character(),
      stringsAsFactors = FALSE
    )
    
    return(list(nodes.net, empty.df))
  }
  

  ## ADD METRICS

  
  common.nr <- edges.net[[4]] * edges.net[[5]]
  
  edges.net$common.nr <- common.nr
  
  max_common <- max(common.nr, na.rm = TRUE)
  
  if (is.na(max_common) || max_common == 0) {
    max_common <- 1
  }
  
  edges.net$common.nr.rel <- common.nr / max_common
  
  edges.net$dot.rel <- round(
    as.numeric(edges.net[[6]]) * 8,
    digits = 0
  )
  

  ## GNPS MASS DIFFERENCE LABELS

  
  if (nettype == "GNPS") {
    
    masses <- DBI::dbGetQuery(
      con,
      sprintf(
        "SELECT compound_id, mass_measured
         FROM ms_compound
         WHERE expid = %d",
        exp.id
      )
    )
    
    masses$compound_id <- as.character(masses$compound_id)
    
    mass_map <- setNames(
      masses$mass_measured,
      masses$compound_id
    )
    
    sub_mass <- mass_map[as.character(edges.net$compid.sub)]
    prod_mass <- mass_map[as.character(edges.net$compid.prod)]
    
    mass_diff <- round(abs(prod_mass - sub_mass), 2)
    
    edges.net$conv.type <- as.character(mass_diff)
    
  } else {
    
    ## KEEP ORIGINAL CSPP LABELS
    
    if (!("conv.type" %in% names(edges.net))) {
      edges.net$conv.type <- "CSPP"
    }
  }
  
  return(list(nodes.net, edges.net))
}