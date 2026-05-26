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
                          nr_of_seq    = 2,
                          output_png   = NULL,   # e.g. "network.png"
                          png_width    = 3000,
                          png_height   = 3000,
                          png_res      = 150) {
  
  library(igraph)
  library(data.table)
  library(scales)
  
  # Open PNG device if requested 
  if (!is.null(output_png)) {
    png(output_png, width = png_width, height = png_height, res = png_res)
    on.exit({ dev.off(); message("Saved: ", output_png) }, add = TRUE)
  }
  
  # Data preparation 
  nodes.net <- as.data.table(net.lst[[1]])
  edges.net <- as.data.table(net.lst[[2]])
  
  nodes.net[, compound_id := as.character(compound_id)]
  edges.net[, compid.sub  := as.character(compid.sub)]
  edges.net[, compid.prod := as.character(compid.prod)]
  
  edges.net <- edges.net[!is.na(compid.sub) & !is.na(compid.prod)]
  
  valid.nodes <- nodes.net$compound_id
  edges.net   <- edges.net[
    compid.sub %in% valid.nodes & compid.prod %in% valid.nodes
  ]
  
  dbkey <- as.character(dbkey)
  if (!(dbkey %in% valid.nodes)) {
    message("compound_id not found: ", dbkey)
    return(NULL)
  }
  
  # BFS traversal 
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
  
  #  Build local edge/node tables 
  edges.netloc <- edges.net[all_idx]
  
  edges.netloc <- edges.netloc[
    , .(
      common.nr.rel = max(as.numeric(common.nr.rel), na.rm = TRUE),
      dot.rel       = max(as.numeric(dot.rel),       na.rm = TRUE),
      ave.dot       = max(as.numeric(ave.dot),       na.rm = TRUE),
      conv.type     = paste(unique(conv.type), collapse = "/")
    ),
    by = .(compid.sub, compid.prod)
  ]
  
  edges.netloc[is.na(common.nr.rel), common.nr.rel := 0.01]
  edges.netloc[is.na(dot.rel),       dot.rel       := 1]
  edges.netloc[is.na(ave.dot),       ave.dot       := 0]
  
  nodes.ids    <- unique(c(edges.netloc$compid.sub, edges.netloc$compid.prod))
  nodes.netloc <- nodes.net[compound_id %in% nodes.ids]
  
  # Build igraph object 
  net <- graph_from_data_frame(
    d        = as.data.frame(edges.netloc),
    vertices = as.data.frame(nodes.netloc),
    directed = TRUE
  )
  
  # Layout 
  # graphopt with generous spring.length spreads nodes apart;
  # norm_coords to [-3, 3] keeps units in a range igraph vertex.size understands.
  # Do NOT norm to [-1,1] — that crushes nodes together.
  set.seed(42)
  lay <- layout_with_graphopt(
    net,
    niter           = 10000,
    charge          = 0.06,    # repulsion between nodes
    mass            = 30,
    spring.length   = 150,     # natural edge length
    spring.constant = 0.005    # soft springs → more freedom to spread
  )
  
  # Rescale to a wide coordinate window so vertex.size stays proportional
  lay <- norm_coords(lay,
                     xmin = -3, xmax = 3,
                     ymin = -3, ymax = 3)
  

  
  # Node aesthetics
  vertex_size <- if ("X2" %in% names(nodes.netloc)) {
    rescale(nodes.netloc$X2, to = c(55, 80))
  } else {
    rep(60, vcount(net))
  }
  
  v_colors <- ifelse(
    V(net)$name == dbkey,
    rgb(1.0, 0.85, 0.20, 0.95),   # yellow  — substrate / hub
    rgb(1.0, 0.65, 0.75, 0.90)    # pink    — product nodes
  )
  
  # Edge aesthetics
  ave_dot_vals  <- as.numeric(edges.netloc$ave.dot)
  ave_dot_power <- ave_dot_vals ^ 3          # exaggerates contrast within 0-1
  
  edge_width  <- rescale(ave_dot_power, to = c(1, 14))   # thin=weak, thick=strong
  edge_alpha  <- rescale(ave_dot_power, to = c(0.15, 0.95))
  edge_colors <- sapply(edge_alpha, function(a) adjustcolor("steelblue", alpha.f = a))
  
  #  Edge curvature 
  src_nodes  <- edges.netloc$compid.sub
  edge_curve <- ave(
    seq_along(src_nodes),
    src_nodes,
    FUN = function(idx) {
      n    <- length(idx)
      vals <- if (n == 1) 0.15 else seq(-0.45, 0.45, length.out = n)
      vals[rank(idx, ties.method = "first")]
    }
  )
  hub_mask <- src_nodes == dbkey
  if (sum(hub_mask) > 1) {
    edge_curve[hub_mask] <- seq(-0.55, 0.55, length.out = sum(hub_mask))
  }
  
  # Edge labels
  edge_labels <- gsub("/", "/\n", edges.netloc$conv.type)
  
  # Plot 
  old_par <- par(mar = c(3, 1, 4, 1))
  on.exit(par(old_par), add = TRUE)
  
  plot(
    net,
    layout  = lay,
    rescale = FALSE,
    xlim    = c(-6, 6),    # padding beyond the -3/3 node range
    ylim    = c(-6, 6),
    
    #  Nodes 
    vertex.color       = v_colors,
    vertex.frame.color = "deeppink4",
    vertex.frame.width = 3,
    vertex.size        = vertex_size,
    
    vertex.label        = V(net)$name,
    vertex.label.cex    = 1.1,
    vertex.label.color  = "black",
    vertex.label.dist   = 0,
    vertex.label.font   = 2,
    
    # Edges 
    edge.color       = edge_colors,
    edge.width       = edge_width,
    edge.arrow.size  = 0.7,
    edge.arrow.width = 2.5,
    edge.curved      = edge_curve,
    
    edge.label       = edge_labels,
    edge.label.cex   = 1.1,
    edge.label.color = "black",
    edge.label.font  = 2,
    edge.label.dist  = 2.2, 
    
    main = nettype
  )
  
  # Legend 
  legend(
    "bottomright",
    legend = c("Substrate compound", "Product compound"),
    pt.bg  = c(rgb(1.0, 0.85, 0.20, 0.95), rgb(1.0, 0.65, 0.75, 0.90)),
    col    = "deeppink4",
    pch    = 21,
    pt.cex = 1.6,
    bty    = "n",
    cex    = 0.75,
    inset  = c(0.01, 0.01)
  )
  
  return(list(nodes.netloc, edges.netloc))
}

