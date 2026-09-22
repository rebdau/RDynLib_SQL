LimKeys_SQL <- function(LCal.lst) {
  
  LCal <- LCal.lst[[1]]
  Reg   <- LCal.lst[[2]]
  
  if (is.null(LCal) || nrow(LCal) == 0) {
    return(list(LCal, Reg))
  }
  
  
  Ind  <- as.integer(LCal[, 1])
  RT   <- as.numeric(LCal[, 2])
  Key3 <- as.numeric(LCal[, 3])
  Key4 <- as.numeric(LCal[, 4])
  Key5 <- as.numeric(LCal[, 5])
  Key6 <- as.numeric(LCal[, 6])
  Key7 <- as.integer(LCal[, 7])
  
  
  # Predicted retention time
  t1 <- if (Reg[3] != 0) {
    pmax(RT - Reg[3], 0)
  } else {
    rep(0, length(RT))
  }
  
  t2 <- if (Reg[5] != 0) {
    pmax(RT - Reg[5], 0)
  } else {
    rep(0, length(RT))
  }
  
  
  RegRT <-
    Reg[1] +
    Reg[2] * RT +
    Reg[4] * t1 +
    Reg[6] * t2
  
  
  Key <- data.frame(
    Ind,
    RT,
    Key3,
    Key4,
    Key5,
    Key6,
    Key7,
    RegRT
  )
  
  
  # Remove RT outliers
  diff <- Key$Key4 - Key$RegRT
  
  AveDiff <- mean(diff, na.rm = TRUE)
  SDDiff  <- stats::sd(diff, na.rm = TRUE)
  
  
  if (!is.na(SDDiff) && SDDiff > 0) {
    
    LowDiff  <- AveDiff - 1.96 * SDDiff
    HighDiff <- AveDiff + 1.96 * SDDiff
    
    keep <- diff >= LowDiff &
      diff <= HighDiff
    
    Key <- Key[keep, , drop = FALSE]
  }
  
  
  # Remove duplicated reference matches
  Key <- Key[order(Key$RT), , drop = FALSE]
  
  if (nrow(Key) > 1) {
    
    grp1 <- paste(
      round(Key$RT, 2),
      round(Key$Key3, 2),
      sep = "_"
    )
    
    drift <- abs(Key$Key4 - Key$RegRT)
    
    Key <- Key[
      ave(
        drift,
        grp1,
        FUN = function(x) x == min(x)
      ) == 1,
      ,
      drop = FALSE
    ]
  }
  
  
  # Remove duplicated target matches
  Key <- Key[order(Key$Key4), , drop = FALSE]
  
  if (nrow(Key) > 1) {
    
    grp2 <- paste(
      round(Key$Key4, 2),
      round(Key$Key5, 2),
      sep = "_"
    )
    
    drift <- abs(Key$Key4 - Key$RegRT)
    
    Key <- Key[
      ave(
        drift,
        grp2,
        FUN = function(x) x == min(x)
      ) == 1,
      ,
      drop = FALSE
    ]
  }
  
  
  Key <- Key[order(Key$RT), , drop = FALSE]
  rownames(Key) <- NULL
  
  
  list(Key, Reg)
}