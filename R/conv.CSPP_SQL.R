#' @title Identify conversion-specific precursor–product pairs and compute CSPP
#'        similarity scores
#'
#' @description
#' conv.CSPP_SQL() identifies candidate precursor–product compound pairs for a
#' given conversion type based on mass differences and retention time constraints,
#' and computes CSPP similarity scores by comparing their MS/MS spectra stored
#' in a SQLite database.
#'
#' @details
#' For each compound in the input data frame, a theoretical product mass is
#' generated using the provided mass difference. Candidate product compounds are
#' selected within a specified mass tolerance and filtered according to their
#' relative elution order.
#'
#' @param inp.x `data.frame()` containing compound information for a single
#'   experiment. Must include the columns \code{compound_id},
#'   \code{mass_measured}, and \code{retention_time}.
#'
#' @param mzdiff `Numeric(1)` mass difference (in Da) corresponding to the
#'   conversion under consideration.
#'
#' @param direc `Integer(1)` elution order constraint between substrate and
#'   product compounds. See Details.
#'
#' @param peakwidth `Numeric(1)` minimum required retention time difference
#'   between substrate and product compounds.
#'
#' @param mzerr `Numeric(1)` mass tolerance (in Da) used to match theoretical
#'   and observed product masses. Default is \code{0.015}.
#'
#' @param data_con `DBIConnection` to an open SQLite database containing
#'   MS/MS spectra and fragment peak information.
#'
#' @return
#' A `data.frame()` containing CSPP similarity metrics for all valid
#' precursor–product pairs identified for the given conversion. Each row
#' corresponds to one candidate conversion and includes fragment ion and
#' neutral loss similarity scores.
#'
#' @author Ahlam Mentag
#'
#' @export
conv.CSPP_SQL <- function(inp.x,
                          mzdiff,
                          direc,
                          peakwidth,
                          mzerr = 0.015,
                          ms2_split) {
  
  Prod.dat <- inp.x[order(inp.x$mass_measured), ]
  Sub.dat  <- Prod.dat
  
  cspp.list <- list()
  idx <- 1
  
  i <- 1
  repeat {
    
    if (i > nrow(Sub.dat)) break
    
    mz.sub  <- Sub.dat$mass_measured[i]
    mz.prod <- mz.sub + mzdiff
    
    prd.low  <- mz.prod - mzerr
    prd.high <- mz.prod + mzerr
    
    j <- 1
    while (j <= nrow(Prod.dat)) {
      
      if (Prod.dat$mass_measured[j] < prd.low) {
        Prod.dat <- Prod.dat[-j, ]
        next
      }
      
      if (Prod.dat$mass_measured[j] <= prd.high) {
        
        out <- targMS2comp_SQL(
          Sub.dat$compound_id[i],
          Prod.dat$compound_id[j],
          ms2_split
        )
        
        if (!is.null(out)) {
          
          rt.sub  <- Sub.dat$retention_time[i]
          rt.prod <- Prod.dat$retention_time[j]
          
          if (
            (direc == 1 && rt.prod < rt.sub - peakwidth) ||
            (direc == 2 && rt.prod > rt.sub + peakwidth) ||
            (direc == 3 && (rt.prod < rt.sub - peakwidth ||
                            rt.prod > rt.sub + peakwidth))
          ) {
            cspp.list[[idx]] <- out
            idx <- idx + 1
          }
        }
        
        j <- j + 1
        next
      }
      
      if (Prod.dat$mass_measured[j] > prd.high) break
    }
    
    if (nrow(Prod.dat) == 0) break
    i <- i + 1
  }
  
  if (length(cspp.list) == 0)
    return(data.frame())
  
  do.call(rbind, cspp.list)
}

