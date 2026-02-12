MultMatch_SQL <- function(MS2list, unk, minIons) {
  end_ms2 <- length(MS2list)
  ms2_sect <- data.frame(COMPID.start = integer(),
                         COMPID.stop = integer(),
                         IonCount = integer())
  
  # Initial check for the full range
  i <- 1
  j <- end_ms2
  Seq_ms2 <- unlist(MS2list[i:j])  # flatten numeric vectors
  CommonIons <- intersect(unk, Seq_ms2)
  
  if (length(CommonIons) > minIons) {
    ms2_sect <- rbind(ms2_sect, data.frame(i, j, length(CommonIons)))
  }
  
  k <- 1
  while (k <= nrow(ms2_sect)) {
    if (ms2_sect[k, 3] <= minIons) {
      k <- k + 1
      next
    }
    
    i <- ms2_sect[k, 1]
    end_ms2 <- ms2_sect[k, 2]
    if ((end_ms2 - i) <= 0) {
      k <- k + 1
      next
    }
    
    step <- as.integer(((end_ms2 - i) / 2) + 0.5)
    j <- i + step
    
    repeat {
      Seq_ms2 <- unlist(MS2list[i:j])
      CommonIons <- intersect(unk, Seq_ms2)
      
      if (length(CommonIons) > minIons) {
        ms2_sect <- rbind(ms2_sect, data.frame(i, j, length(CommonIons)))
      }
      
      i <- j + 1
      j <- j + step
      if (j > end_ms2) {
        if (i < end_ms2) j <- end_ms2 else break
      }
    }
    
    # Remove the current row after splitting
    ms2_sect <- ms2_sect[-k, , drop = FALSE]
  }
  
  return(ms2_sect)
}
