#' @title display the network graph
#' 
#' @description The Local.net_SQL() function will select those nodes and edges 
#' belonging to the local network of the selected COMPID and display 
#' the network graph.
#' 
#' @param net.lst 'list()' list of nodes and edges.
#' 
#' @param dbkey 'number(1)' refers to the compound_id.
#' 
#' @param nettype 'character(1)' the graph type the user want to display, 
#' it could be either 
#' 
#' @param nr_of_seq 'number(1)' size of local network, i.e., number of 
#' subsequent edges starting from the node representing the selected COMPID.
#' 
#' @return list of edges and nodes belonging to the local network of the 
#'         selected COMPID and plot the network.
#' 
#' @import RColorBrewer
#' @import igraph
#' 
#' @author Ahlam Mentag
#' @export
Local.net_SQL <- function(net.lst, dbkey, nettype,
                          nr_of_seq = 2) {
  
  library(igraph)
  library(data.table)
  library(scales)
  
  nodes.net <- as.data.table(net.lst[[1]])
  edges.net <- as.data.table(net.lst[[2]])
  
  # Standardize types 
  nodes.net[, compound_id := as.character(compound_id)]
  edges.net[, compid.sub  := as.character(compid.sub)]
  edges.net[, compid.prod := as.character(compid.prod)]
  
  # Remove invalid rows 
  edges.net <- edges.net[!is.na(compid.sub) & !is.na(compid.prod)]
  
  # Keep only valid nodes 
  valid.nodes <- nodes.net$compound_id
  edges.net   <- edges.net[
    compid.sub %in% valid.nodes & compid.prod %in% valid.nodes
  ]
  
  dbkey <- as.character(dbkey)
  if (!(dbkey %in% valid.nodes)) {
    message("compound_id not found: ", dbkey)
    return(NULL)
  }
  
  current       <- dbkey
  visited_edges <- integer()
  visited_nodes <- dbkey
  
  for (i in seq_len(nr_of_seq)) {
    idx <- which(edges.net$compid.sub %in% current)
    if (length(idx) == 0) break
    visited_edges <- c(visited_edges, idx)
    new_nodes     <- setdiff(unique(edges.net$compid.prod[idx]), visited_nodes)
    if (length(new_nodes) == 0) break
    visited_nodes <- unique(c(visited_nodes, new_nodes))
    current       <- new_nodes
  }
  
  all_idx <- unique(visited_edges)
  if (length(all_idx) == 0) {
    message("no conversions")
    return(NULL)
  }
  
  edges.netloc <- edges.net[all_idx]
  
  # Aggregate multiple conversions
  edges.netloc <- edges.netloc[
    , .(
      common.nr.rel = max(as.numeric(common.nr.rel), na.rm = TRUE),
      dot.rel       = max(as.numeric(dot.rel),       na.rm = TRUE),
      ave.dot       = max(as.numeric(ave.dot),       na.rm = TRUE),
      conv.type     = paste(unique(conv.type), collapse = "/")
    ),
    by = .(compid.sub, compid.prod)
  ]
  
  # Fix NAs
  edges.netloc[is.na(common.nr.rel), common.nr.rel := 0.01]
  edges.netloc[is.na(dot.rel),       dot.rel       := 1]
  edges.netloc[is.na(ave.dot),       ave.dot       := 0]
  
  # Nodes subset
  nodes.ids    <- unique(c(edges.netloc$compid.sub, edges.netloc$compid.prod))
  nodes.netloc <- nodes.net[compound_id %in% nodes.ids]
  
  # Build graph
  net <- graph_from_data_frame(
    d        = as.data.frame(edges.netloc),
    vertices = as.data.frame(nodes.netloc),
    directed = TRUE
  )
  


  set.seed(42)
  lay <- layout_with_kk(net)                        # good global spacing
  lay <- layout_with_fr(                            # local repulsion pass
    net,
    coords     = lay,                               # warm-start from KK
    niter      = 3000,
    start.temp = vcount(net) * 3,
    weights    = rep(0.5, ecount(net))              
  )
  lay <- norm_coords(lay, xmin = -1.8, xmax = 1.8,
                     ymin = -1.8, ymax = 1.8)
  
 
  vertex_size <- if ("X2" %in% names(nodes.netloc)) {
    rescale(nodes.netloc$X2, to = c(68, 72))
  } else {
    rep(42, vcount(net))
  }
  
  v_colors <- ifelse(
    V(net)$name == dbkey,
    rgb(1.0, 0.85, 0.20, 0.95), 
    rgb(1.0, 0.65, 0.75, 0.90)   
  )
  

  edge_width  <- rescale(edges.netloc$ave.dot, to = c(1.5, 6))
  edge_colors <- adjustcolor("steelblue", alpha.f = 0.60)
  

  src_nodes   <- edges.netloc$compid.sub
  edge_curve  <- ave(
    seq_along(src_nodes),
    src_nodes,
    FUN = function(idx) {
      n   <- length(idx)
      seq <- seq(-0.35, 0.35, length.out = n)   # spread from -0.35 to +0.35
      seq[rank(idx)]
    }
  )
  # Hub edges: slightly stronger curve to open up the spoke bundle
  edge_curve[src_nodes == dbkey] <- seq(-0.45, 0.45,
                                        length.out = sum(src_nodes == dbkey))
  

  # Shorten long compound type strings to reduce overlap
  edge_labels <- edges.netloc$conv.type
  # Wrap labels that have "/" separators onto two lines
  edge_labels <- gsub("/", "/\n", edge_labels)
  

  old_par <- par(mar = c(2, 1, 3, 1))
  on.exit(par(old_par))
  
  plot(
    net,
    layout  = lay,
    rescale = FALSE,
    xlim    = c(-2.2, 2.2),
    ylim    = c(-2.2, 2.2),
    
    # Nodes
    vertex.color       = v_colors,
    vertex.frame.color = "deeppink4",
    vertex.size        = vertex_size,
    
    # Labels inside circles
    vertex.label        = V(net)$name,
    vertex.label.cex    = 0.65,
    vertex.label.color  = "black",
    vertex.label.dist   = 0,
    vertex.label.degree = 0,
    
    # Edges
    edge.arrow.size = 0.4,
    edge.color      = edge_colors,
    edge.width      = edge_width,
    edge.curved     = edge_curve,      
    
    # Edge labels: larger + dark for readability
    edge.label       = edge_labels,
    edge.label.cex   = 0.58,           
    edge.label.color = "grey15",       
    
    main = nettype
  )
  
  legend(
    "bottomright",
    legend = c("Query compound", "Downstream compound"),
    pt.bg  = c(rgb(1.0, 0.85, 0.20, 0.95), rgb(1.0, 0.65, 0.75, 0.90)),
    col    = "deeppink4",
    pch    = 21,
    pt.cex = 1.6,
    bty    = "n",
    cex    = 0.75
  )
  
  return(list(nodes.netloc, edges.netloc))
}