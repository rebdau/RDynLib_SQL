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
Local.net_SQL <- function(net.lst,
                          dbkey,
                          nr_of_seq = 2,
                          output_png = NULL,
                          png_width = 6000,
                          png_height = 3000,
                          png_res = 200) {
  
  library(igraph)
  library(data.table)
  library(scales)
  

  ## SAVE BOTH NETWORKS IN SAME PNG

  
  if (!is.null(output_png)) {
    png(output_png,
        width = png_width,
        height = png_height,
        res = png_res)
    
    on.exit({
      dev.off()
      message("Saved: ", output_png)
    }, add = TRUE)
  }
  
  old_par <- par(no.readonly = TRUE)
  on.exit(par(old_par), add = TRUE)
  
  ## TWO PANELS
  par(mfrow = c(1, 2), mar = c(3, 1, 4, 1))
  
  ## LOOP OVER NETWORKS
  
  for (k in seq_along(net.lst)) {
    
    current_net <- net.lst[[k]]
    
    nodes.net <- as.data.table(current_net[[1]])
    edges.net <- as.data.table(current_net[[2]])
    
    if (nrow(edges.net) == 0) {
      plot.new()
      title(main = paste("Network", k, "- empty"))
      next
    }
    
    nodes.net[, compound_id := as.character(compound_id)]
    
    edges.net[, compid.sub  := as.character(compid.sub)]
    edges.net[, compid.prod := as.character(compid.prod)]
    
    edges.net <- edges.net[
      !is.na(compid.sub) & !is.na(compid.prod)
    ]
    
    valid.nodes <- nodes.net$compound_id
    
    edges.net <- edges.net[
      compid.sub %in% valid.nodes &
        compid.prod %in% valid.nodes
    ]
    
    dbkey <- as.character(dbkey)
    
    if (!(dbkey %in% valid.nodes)) {
      plot.new()
      title(main = paste("compound_id not found:", dbkey))
      next
    }

    
    current <- dbkey
    visited_edges <- integer()
    visited_nodes <- dbkey
    
    for (i in seq_len(nr_of_seq)) {
      
      idx <- which(edges.net$compid.sub %in% current)
      
      if (length(idx) == 0)
        break
      
      visited_edges <- c(visited_edges, idx)
      
      new_nodes <- setdiff(
        unique(edges.net$compid.prod[idx]),
        visited_nodes
      )
      
      if (length(new_nodes) == 0)
        break
      
      visited_nodes <- unique(c(visited_nodes, new_nodes))
      current <- new_nodes
    }
    
    all_idx <- unique(visited_edges)
    
    if (length(all_idx) == 0) {
      plot.new()
      title(main = "No conversions")
      next
    }
    
    edges.netloc <- edges.net[all_idx]
    
    ## AGGREGATION
    
    if ("conv.type" %in% names(edges.netloc)) {
      
      edges.netloc <- edges.netloc[
        ,
        .(
          common.nr.rel = max(as.numeric(common.nr.rel), na.rm = TRUE),
          dot.rel       = max(as.numeric(dot.rel), na.rm = TRUE),
          ave.dot       = max(as.numeric(ave.dot), na.rm = TRUE),
          conv.type     = paste(unique(conv.type), collapse = "/")
        ),
        by = .(compid.sub, compid.prod)
      ]
      
      edge_labels <- edges.netloc$conv.type
      
    } else {
      
      edges.netloc <- edges.netloc[
        ,
        .(
          common.nr.rel = max(as.numeric(common.nr.rel), na.rm = TRUE),
          dot.rel       = max(as.numeric(dot.rel), na.rm = TRUE),
          ave.dot       = max(as.numeric(ave.dot), na.rm = TRUE)
        ),
        by = .(compid.sub, compid.prod)
      ]
      
      edge_labels <- round(edges.netloc$ave.dot, 3)
    }
    
    edges.netloc[is.na(common.nr.rel), common.nr.rel := 0.01]
    edges.netloc[is.na(dot.rel),       dot.rel       := 1]
    edges.netloc[is.na(ave.dot),       ave.dot       := 0]
    
    nodes.ids <- unique(c(
      edges.netloc$compid.sub,
      edges.netloc$compid.prod
    ))
    
    nodes.netloc <- nodes.net[
      compound_id %in% nodes.ids
    ]
    
    ## GRAPH
    
    net <- graph_from_data_frame(
      d = as.data.frame(edges.netloc),
      vertices = as.data.frame(nodes.netloc),
      directed = TRUE
    )
    
    set.seed(42)
    
    lay <- layout_with_graphopt(
      net,
      niter = 20000,
      charge = 0.2,
      mass = 50,
      spring.length = 300,
      spring.constant = 0.002
    )
    lay <- norm_coords(
      lay,
      xmin = -8, xmax = 8,
      ymin = -8, ymax = 8
    )
    
    ## NODE STYLE
    
    vertex_size <- if ("X2" %in% names(nodes.netloc)) {
      rescale(nodes.netloc$X2, to = c(77, 100))
    } else {
      rep(150, vcount(net))
    }
    
    v_colors <- ifelse(
      V(net)$name == dbkey,
      rgb(1.0, 0.85, 0.20, 0.95),
      rgb(1.0, 0.65, 0.75, 0.90)
    )
    
    ## EDGE STYLE
    
    ave_dot_vals <- as.numeric(edges.netloc$ave.dot)
    ave_dot_power <- ave_dot_vals ^ 3
    
    edge_width <- rescale(
      ave_dot_power,
      to = c(1, 14)
    )
    
    edge_alpha <- rescale(
      ave_dot_power,
      to = c(0.15, 0.95)
    )
    
    edge_colors <- sapply(
      edge_alpha,
      function(a)
        adjustcolor("steelblue", alpha.f = a)
    )
    
    src_nodes <- edges.netloc$compid.sub
    
    edge_curve <- ave(
      seq_along(src_nodes),
      src_nodes,
      FUN = function(idx) {
        n <- length(idx)
        
        vals <- if (n == 1) {
          0.15
        } else {
          seq(-0.45, 0.45, length.out = n)
        }
        
        vals[rank(idx, ties.method = "first")]
      }
    )
    
    ## TITLE
    
    plot_title <- if (k == 1) {
      "CSPP Network"
    } else {
      "GNPS Network"
    }
    
    ## PLOT
    
    plot(
      net,
      layout = lay,
      rescale = FALSE,
      xlim = c(-10, 10),
      ylim = c(-10, 10),
      
      vertex.color = v_colors,
      vertex.frame.color = "deeppink4",
      vertex.frame.width = 3,
      vertex.size = vertex_size,
      
      vertex.label = V(net)$name,
      vertex.label.cex = 1.1,
      vertex.label.color = "black",
      vertex.label.font = 2,
      edge.label.dist = 1.5,
      
      edge.color = edge_colors,
      edge.width = edge_width,
      edge.arrow.size = 1.2,
      edge.arrow.width = 3,
      edge.curved = edge_curve,
      
      edge.label = edge_labels,
      edge.label.cex = 0.9,
      edge.label.color = "black",
      edge.label.font = 2,

      main = plot_title
    )
    
    ## LEGEND
    
    legend(
      "bottomright",
      legend = c(
        "Seed / substrate compound",
        "Product compounds"
      ),
      
      pt.bg = c(
        rgb(1.0, 0.85, 0.20, 0.95),
        rgb(1.0, 0.65, 0.75, 0.90)
      ),
      
      col = "deeppink4",
      pch = 21,
      pt.cex = 2,
      bty = "n",
      cex = 1
    )}
  
  return(net.lst)
}
