subMultMatch_SQL <- function(msnSect, unk, MS2list, minIons) {
  ENTRYretained <- integer()
  
  for (i in 1:nrow(msnSect)) {
    ENTRYretained <- c(ENTRYretained, seq(msnSect[i, 1], msnSect[i, 2]))
  }
  ENTRYretained <- sort(unique(ENTRYretained))
  
  i <- 1
  while (i <= length(ENTRYretained)) {
    Seq_ms2 <- MS2list[[as.character(ENTRYretained[i])]]  # already numeric vector
    CommonIons <- intersect(unk, Seq_ms2)
    
    if (length(CommonIons) > minIons) {
      i <- i + 1
    } else {
      ENTRYretained <- ENTRYretained[-i]
    }
  }
  
  return(as.integer(ENTRYretained))
}
