#' @title Piecewise Robust Regression Between FTMS and QTOF Retention Times
#'
#' @description
#' 
#' `RegressionPie_LCalign_SQL()` performs piecewise regression between FTMS and
#' Synapt QTOF retention times obtained after local alignment (`LCal`).  
#' 
#' Because chromatograms may come from different LC systems, deviations from 
#' linearity are common. Therefore the function evaluates:
#' - a model without knots (simple linear regression),
#' - piecewise linear models with one knot,
#' - piecewise linear models with two knots.  
#' 
#' For each model, the adjusted R^2 is computed, and the best model is selected.
#' The selected model is then refitted using robust regression (`rlm`) to reduce
#' sensitivity to outliers.  
#'
#' The `startpoint` parameter specifies the minimum FT retention time from which
#' knot positions can be considered, and iterates only over the unique retention 
#' times actually present in LCal.
#'
#' @param LCal A matrix or data frame representing local alignment results, 
#'        where column 2 contains FTMS retention times and column 4 contains 
#'        QTOF retention times.
#' @param startpoint Numeric. Minimum FTMS retention time at which knots may be 
#'        placed. Defaults to 1.
#'
#' @return A numeric vector containing six regression parameters:  
#'  
#'    - `exp.intercept`: intercept of the robust regression,
#'    - `exp.slope`: main slope before first knot,
#'    - `k1_best`: first knot position (0 if no knot),
#'    - `exp.t1`: slope adjustment after first knot,
#'    - `k2_best`: second knot position (0 if not used),
#'    - `exp.t2`: slope adjustment after second knot.
#'   
#'
#' @author Ahlam Mentag
#'
#' @import MASS
#'
#' @export
RegressionPie_LCalign_SQL <- function(LCal, startpoint = 1) {
  
  x <- as.numeric(LCal[, 2])  # FT retention times
  y <- as.numeric(LCal[, 4])  # Target retention times
  
  # Remove NA
  valid <- which(!is.na(x) & !is.na(y))
  x <- x[valid]
  y <- y[valid]
  
  
  # Helper to extract adjusted R2 without
  # "essentially perfect fit" warnings
  get_adj_r2 <- function(model) {
    suppressWarnings(
      summary(model)$adj.r.squared
    )
  }
  
  
  # Prepare output vectors
  knottime1 <- c()
  knottime2 <- c()
  VarExp <- c()
  
  
  # Initial model: no knots
  out <- lm(y ~ x)
  
  knottime1 <- c(knottime1, 0)
  knottime2 <- c(knottime2, 0)
  VarExp <- c(
    VarExp,
    get_adj_r2(out)
  )
  
  
  unique_x <- sort(
    unique(x[x >= startpoint])
  )
  
  
  # One-knot models
  for (k1 in unique_x) {
    
    t1 <- ifelse(
      x > k1,
      x - k1,
      0
    )
    
    out <- lm(
      y ~ x + t1
    )
    
    knottime1 <- c(
      knottime1,
      k1
    )
    
    knottime2 <- c(
      knottime2,
      0
    )
    
    VarExp <- c(
      VarExp,
      get_adj_r2(out)
    )
  }
  
  
  # Two-knot models
  if (length(unique_x) > 1) {
    
    for (i in seq_len(length(unique_x) - 1)) {
      
      k1 <- unique_x[i]
      
      t1 <- ifelse(
        x > k1,
        x - k1,
        0
      )
      
      
      for (j in (i + 1):length(unique_x)) {
        
        k2 <- unique_x[j]
        
        t2 <- ifelse(
          x > k2,
          x - k2,
          0
        )
        
        out <- lm(
          y ~ x + t1 + t2
        )
        
        
        knottime1 <- c(
          knottime1,
          k1
        )
        
        knottime2 <- c(
          knottime2,
          k2
        )
        
        VarExp <- c(
          VarExp,
          get_adj_r2(out)
        )
      }
    }
  }
  
  
  # Select best model
  best.model <- data.frame(
    knottime1,
    knottime2,
    VarExp
  )
  
  best.model <- best.model[
    order(
      best.model$VarExp,
      decreasing = TRUE
    ),
    ,
    drop = FALSE
  ]
  
  
  k1_best <- best.model$knottime1[1]
  k2_best <- best.model$knottime2[1]
  
  
  t1 <- ifelse(
    x > k1_best,
    x - k1_best,
    0
  )
  
  t2 <- ifelse(
    x > k2_best,
    x - k2_best,
    0
  )
  
  
  # Robust final regression
  if (k1_best == 0 && k2_best == 0) {
    
    out <- MASS::rlm(
      y ~ x
    )
    
    co <- summary(out)$coef
    
    exp.intercept <- co[1, 1]
    exp.slope <- co[2, 1]
    exp.t1 <- 0
    exp.t2 <- 0
    
    
  } else if (k2_best == 0) {
    
    out <- MASS::rlm(
      y ~ x + t1
    )
    
    co <- summary(out)$coef
    
    exp.intercept <- co[1, 1]
    exp.slope <- co[2, 1]
    exp.t1 <- co[3, 1]
    exp.t2 <- 0
    
    
  } else {
    
    out <- MASS::rlm(
      y ~ x + t1 + t2
    )
    
    co <- summary(out)$coef
    
    exp.intercept <- co[1, 1]
    exp.slope <- co[2, 1]
    exp.t1 <- co[3, 1]
    exp.t2 <- co[4, 1]
  }
  
  
  return(
    c(
      exp.intercept,
      exp.slope,
      k1_best,
      exp.t1,
      k2_best,
      exp.t2
    )
  )
}