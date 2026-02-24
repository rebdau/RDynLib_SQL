#' @title Local networks for a selected COMPID 
#'
#' @description the Overall.net_SQL() function plots the CSPP or GNPS-like 
#' local networks, it iterates through the net.addit() and Local.net_SQL()
#' functions. The net.addit() function is responsible for generating lists 
#' of network nodes and edges by the net.nodes() and either the allCSPP() or 
#' allGNPS() functions, respectively. This generates the nodes and edges for 
#' the experiment-based global network. Only edges are retained that pass 
#' certain CID spectral similarity thresholds. The Local.net_SQL() function will 
#' select those nodes and edges belonging to the local network of the selected
#'COMPID and display the network graph.
#' 
#' @param sql_path 'character(1)' path to the sqlite database.
#' 
#' @param exp.id 'number(1)' experiment to use for network generation.
#' 
#' @param dbkey 'number(1)' compound_id, serving as the starting point of the 
#' molecular networking.Only nodes and edges connected to this COMPID within 
#' the specified number of sequential connections (`nr_of_seq`) are included 
#' in the displayed CSPP and GNPS local networks.
#' 
#' @param nr_of_seq 'number(1)' size of local network, i.e., number of 
#' subsequent edges starting from the node representing the selected COMPID.
#' 
#' @param thr1 'number(1)' minimum number of product ions (varying between 0 
#' and 1) of one CID spectrum that can be traced in another CID spectrum,
#' set on 1 by default.
#' 
#' @param thr2 'number(1)' dot product threshold for the common ions 
#' between two CID spectra, set on 0.9 by default.
#' 
#' @param thr3 'number(1)' average of the dot products obtained for the common 
#' product ions and the common neutral losses, by default set on 0.4.
#' 
#' @return returns a list of the cspp and gnps-like information and plots 
#'         the final network of the given compound_id in a given experiment.
#'         
#' @author Ahlam Mentag
#' 
#' @export
Overall.net_SQL <- function(sql_path, exp.id, dbkey, nr_of_seq = 2, thr1 = 1,
                            thr2 = 0.9, thr3 = 0.4) {
  
  oldpar <- par(no.readonly = TRUE)
  on.exit({
    try(par(oldpar), silent = TRUE)
  }, add = TRUE)
  
  par(mfrow = c(1,2), mar = c(1,1,1,1))
  
  ##CSPP
  net.lst1 <- net.addit_SQL(sql_path = sql_path,
                            exp.id = exp.id,
                            nettype = "CSPP",
                            thr1 = thr1,
                            thr2 = thr2,
                            thr3 = thr3)
  
  locNET1.list <- Local.net_SQL(net.lst1,
                                dbkey = dbkey,
                                nettype = "CSPP",
                                nr_of_seq = nr_of_seq)
  
  
  ##GNPS 
  net.lst2 <- net.addit_SQL(sql_path = sql_path,
                            exp.id = exp.id,
                            nettype = "GNPS",
                            thr1 = thr1,
                            thr2 = thr2,
                            thr3 = thr3)
  
  locNET2.list <- Local.net_SQL(net.lst2,
                                dbkey = dbkey,
                                nettype = "GNPS",
                                nr_of_seq = nr_of_seq)
  
  return(list(locNET1.list, locNET2.list))
}
