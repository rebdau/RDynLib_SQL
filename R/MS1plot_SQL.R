MS1plot_SQL <- function(dbkey, exp, con, proc, fc, finlist,
                        lc.err, mz.err, prcx) {
  

  query <- "
    SELECT node_name, compound_id, dbkey
    FROM ms_compound
    WHERE expid = ? AND dbkey = ?
  "
  
  dynlb <- DBI::dbGetQuery(con, query,
                           params = list(exp, dbkey))
  
  if (nrow(dynlb) == 0)
    stop("No entry found for this dbkey and experiment.")
  
  nodename <- dynlb$node_name[1]
  compound_id <- dynlb$compound_id[1]
  

  ## Prepare processed MS1 table

  if (nrow(proc) == 0)
    stop("proc table is empty.")
  
  proc <- proc[order(proc[,1]),]
  
  sel <- which(proc[,3] %in% nodename)
  
  if (length(sel) == 0)
    stop("Node not found in processed MS1 data.")
  
  sel <- sel[1]  # ensure single anchor row
  
  rettime <- proc[sel, 1]
  lowRt  <- rettime - lc.err
  highRt <- rettime + lc.err
  

  ## 3. Extract RT window 

  
  rt_filter <- proc[,1] >= lowRt & proc[,1] <= highRt
  node_filter <- proc[,3] %in% nodename
  
  window_rows <- which(rt_filter & node_filter)
  
  if (length(window_rows) == 0)
    stop("No MS1 peaks found in RT window.")
  
  mzMS1  <- proc[window_rows, 2]
  namMS1 <- proc[window_rows, 3]
  
  abunMS1 <- apply(proc[window_rows, 5:ncol(proc), drop = FALSE],
                   1,
                   function(x) mean(as.numeric(x), na.rm = TRUE))

  
  if (length(abunMS1) == 0)
    stop("No intensity values found.")
  
  if (abunMS1[1] == 0)
    stop("Anchor intensity is zero — cannot compute relative abundance.")
  
  

  ##  Calculations
  
  relabun <- round(100 * abunMS1 / abunMS1[1], 1)
  mzdiff  <- round(mzMS1 - mzMS1[1], 4)
  
  compID <- rep(compound_id, length(mzMS1))
  
  

  ## Feature annotation

  
  FeatName <- rep(NA, length(mzMS1))
  
  for (i in seq_along(mzMS1)) {
    
    lowthr  <- abs(mzdiff[i]) - mz.err
    highthr <- abs(mzdiff[i]) + mz.err
    
    sl <- which(finlist[[6]][,1] > lowthr &
                  finlist[[6]][,1] < highthr)
    
    if (length(sl) == 1) {
      FeatName[i] <- finlist[[6]][sl,2]
    }
  }
  
  

  ## Ion-neutral combinations

  
  massMS1 <- as.integer(round(mzMS1 + 1.0073, 2) * 100)
  ion_neut <- outer(massMS1, massMS1, "+")
  
  Combin <- rep(NA, length(mzMS1))
  
  ion_neut_ind <- which(ion_neut %in% massMS1)
  
  if (length(ion_neut_ind) > 0) {
    
    for (idx in ion_neut_ind) {
      
      row <- ceiling(idx / length(massMS1))
      col <- ifelse(idx %% length(massMS1) == 0,
                    length(massMS1),
                    idx %% length(massMS1))
      
      target <- which(massMS1 == ion_neut[row, col])
      
      if (length(target) == 1) {
        Combin[target] <-
          paste(round(mzMS1[row], 2), "/",
                round(mzMS1[col], 2), sep = "")
      }
    }
  }
  

  if (all(is.na(fc))) {
    
    dMS1 <- data.frame(
      mzMS1, mzdiff, relabun,
      FeatName, compID, Combin
    )
    
    colnames(dMS1) <- c(
      "m/z", "m/z diff", "Rel.Int.",
      "Type", "COMPID", "Combin"
    )
    
  } else {
    
    corrcoef <- MS1corr(fc, proc, namMS1)
    
    dMS1 <- data.frame(
      mzMS1, mzdiff, relabun,
      FeatName, compID, Combin,
      corrcoef
    )
    
    colnames(dMS1) <- c(
      "m/z", "m/z diff", "Rel.Int.",
      "Type", "COMPID", "Combin",
      "R2"
    )
  }
  
  

  ## Plot

  
  adj_intens   <- max(relabun) * 1.1
  adj_prod_min <- min(mzMS1) * 0.9
  adj_prod_max <- max(mzMS1) * 1.1
  
  par(cex = prcx)
  
  plot(mzMS1, relabun,
       type = "h",
       xlab = "m/z",
       ylab = "relative intensity",
       main = paste("MS1: m/z", round(mzMS1[1], 4)),
       xlim = c(adj_prod_min, adj_prod_max),
       ylim = c(0, adj_intens))
  
  text(mzMS1, relabun,
       round(mzMS1, 3), pos = 3)
  
  text(mzMS1, relabun,
       round(mzdiff, 3),
       pos = 3, offset = 1.5, col = 2)
  
  return(dMS1)
}
