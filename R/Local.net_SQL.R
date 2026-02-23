Local.net_SQL <- function(net.lst, dbkey, nettype, nr_of_seq = 2) {
  
  if (is.null(nr_of_seq)) nr_of_seq <- 2
  
  nodes.net <- net.lst[[1]]
  edges.net <- net.lst[[2]]
  
  
  edges.net <- edges.net[
    complete.cases(edges.net[, c("compid.sub","compid.prod")]), ]
  
  valid.nodes <- nodes.net$compound_id
  edges.net <- edges.net[
    edges.net$compid.sub %in% valid.nodes &
      edges.net$compid.prod %in% valid.nodes, ]
  
  if (!dbkey %in% nodes.net$compound_id) {
    message("compound_id not found: ", dbkey)
    return(NULL)
  }
  
  
  sequ <- dbkey
  sequ.cont <- integer()
  counter <- 1
  
  repeat {
    sequ.row <- which(edges.net$compid.sub %in% sequ)
    if (length(sequ.row) == 0) break
    
    sequ.cont <- c(sequ.cont, sequ.row)
    sequ <- edges.net$compid.prod[sequ.row]
    
    if (counter == nr_of_seq) break
    counter <- counter + 1
  }
  
  edges.netloc <- edges.net[sequ.cont, , drop = FALSE]
  
  if (nrow(edges.netloc) == 0) {
    message("no conversions")
    return(NULL)
  }
  
  nodes.ids <- unique(c(edges.netloc$compid.sub,
                        edges.netloc$compid.prod))
  
  nodes.netloc <- nodes.net[
    nodes.net$compound_id %in% nodes.ids, , drop = FALSE]
  
  
  # dot.rel
  edges.netloc$dot.rel <- as.numeric(edges.netloc$dot.rel)
  edges.netloc$dot.rel[is.na(edges.netloc$dot.rel)] <- 1
  
  # common.nr.rel
  edges.netloc$common.nr.rel <- as.numeric(edges.netloc$common.nr.rel)
  edges.netloc$common.nr.rel[is.na(edges.netloc$common.nr.rel)] <- 0.01
  
  # mass.diff
  if ("mass.diff" %in% colnames(edges.netloc)) {
    edges.netloc$mass.diff <- as.numeric(edges.netloc$mass.diff)
    edges.netloc$mass.diff[is.na(edges.netloc$mass.diff)] <- 0
    edges.netloc$mass.diff <- round(edges.netloc$mass.diff, 3)
  }
  
  # ---- Création du graphe ----
  library(igraph)
  library(RColorBrewer)
  
  pal3 <- brewer.pal(8, "BuGn")
  
  net <- graph_from_data_frame(
    d = edges.netloc,
    vertices = nodes.netloc,
    directed = TRUE
  )
  
  
  if (nettype == "CSPP") {
    edge_labels <- as.character(E(net)$conv.type)
    main_title <- "CSPP"
  } else {
    edge_labels <- format(E(net)$mass.diff, nsmall = 3)
    main_title <- "GNPS-like"
  }
  
  # Plot 
  plot(net,
       edge.arrow.size = 0.6,
       vertex.label.cex = 0.7,
       vertex.size = V(net)$X2,
       edge.label = edge_labels,
       edge.label.cex = 0.8,
       edge.label.color = "black",
       edge.label.dist = 0.5,
       edge.color = pal3[pmax(E(net)$dot.rel,1)],
       edge.width = 40*pmax(E(net)$common.nr.rel,0.01),
       main = main_title
  )
  
  return(list(nodes.netloc, edges.netloc))
}