CommonDotProd <- function(ac, bd, tol = 0.01) {
  
  # ac, bd must have columns: mz (or nloss), rint
  
  mz_a <- as.numeric(ac[[1]])
  mz_b <- as.numeric(bd[[1]])
  
  matches <- which(abs(outer(mz_a, mz_b, "-")) <= tol, arr.ind = TRUE)
  
  # Ensure matches is always a matrix
  if (is.null(dim(matches))) {
    matches <- matrix(matches, ncol = 2)
  }
  
  if (nrow(matches) == 0) {
    return(list(0, 0))
  }
  
  matched_pairs <- split(as.data.frame(matches), matches[,1])
  
  v <- numeric(length(matched_pairs))
  w <- numeric(length(matched_pairs))
  
  i <- 1
  for (idx in matched_pairs) {
    
    idx <- as.matrix(idx) 
    
    a_idx <- idx[1, 1]
    b_idx <- idx[, 2]
    
    v[i] <- as.numeric(ac[a_idx, 2])
    w[i] <- as.numeric(max(bd[b_idx, 2]))
    
    i <- i + 1
  }
  
  common <- length(v)
  
  theta <- sum(v * w) / (sqrt(sum(v^2)) * sqrt(sum(w^2)))
  
  return(list(common, theta))
}