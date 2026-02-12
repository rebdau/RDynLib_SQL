Local.net_SQL <- function(net.lst, dbkey, nettype, nr_of_seq = 2) {
  
  if (is.null(nr_of_seq)) nr_of_seq <- 2
  
  nodes.net <- net.lst[[1]]  # nodes from net.addit_SQL
  edges.net <- net.lst[[2]]  # edges from net.addit_SQL
  
  # Check if dbkey exists in nodes
  if (!dbkey %in% nodes.net$compound_id) {
    message("compound_id not found in nodes: ", dbkey)
    return(list(nodes.netloc = NULL, edges.netloc = NULL))
  }
  
  # Build local network starting from dbkey
  sequ <- dbkey
  sequ.cont <- integer()
  counter <- 1
  
  repeat {
    # Find edges where the source matches sequ (compound_id)
    sequ.row <- which(edges.net[,1] %in% sequ)
    if (length(sequ.row) == 0) break  # no further edges
    
    sequ.cont <- c(sequ.cont, sequ.row)
    sequ <- edges.net[sequ.row, 2]  # next compounds in the sequence
    
    if (counter == nr_of_seq) break
    counter <- counter + 1
  }
  
  edges.netloc <- edges.net[sequ.cont, , drop = FALSE]
  
  if (nrow(edges.netloc) == 0) {
    message("No conversions found for compound_id: ", dbkey)
    return(list(nodes.netloc = NULL, edges.netloc = NULL))
  }
  
  # Filter nodes involved in this local network
  nodes.loc <- sort(unique(c(edges.netloc[,1], edges.netloc[,2])))
  nodes.netloc <- nodes.net[nodes.net$compound_id %in% nodes.loc, , drop = FALSE]
  
  # Plot network
  library(RColorBrewer)
  pal3 <- brewer.pal(8, "BuGn")
  library(igraph)
  
  net <- graph.data.frame(edges.netloc, nodes.netloc, directed = TRUE)
  
  if (nettype == "CSPP") {
    plot(net,
         edge.arrow.size = 0.6,
         vertex.label.cex = 0.7,
         vertex.size = V(net)$X2,
         edge.label = E(net)$conv.type,
         edge.label.cex = 0.7,
         edge.label.color = 1,
         edge.color = pal3[E(net)$dot.rel],
         edge.width = 40 * E(net)$common.nr.rel,
         main = "CSPP")
  } else if (nettype == "GNPS") {
    plot(net,
         edge.arrow.size = 0.6,
         vertex.label.cex = 0.7,
         vertex.size = V(net)$X2,
         edge.label = E(net)$mass.diff,
         edge.label.cex = 0.7,
         edge.label.color = 1,
         edge.color = pal3[E(net)$dot.rel],
         edge.width = 40 * E(net)$common.nr.rel,
         main = "GNPS-like")
  } else {
    message("Define type of net: CSPP or GNPS")
  }
  
  return(list(nodes.netloc = nodes.netloc, edges.netloc = edges.netloc))
}
