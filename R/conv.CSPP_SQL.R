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
#'@import data.table
#' @author Ahlam Mentag
#'
#'
#' @export
NULL

conv.CSPP_SQL <- function(inp.x,
                          mzdiff,
                          direc,
                          peakwidth,
                          mzerr = 1,
                          ms2_split,
                          IntThres = 0,
                          mz_tol = 0.3,
                          dot_thresh = 0,
                          min_common = 0) {
  
  inp.x <- as.data.table(inp.x)
  inp.x$compound_id <- as.character(inp.x$compound_id)
  
  # ---- SELF JOIN BASE ----
  pairs <- CJ(sub  = inp.x$compound_id,
              prod = inp.x$compound_id)
  
  # ---- JOIN metadata ----
  pairs <- merge(
    pairs,
    data.table(sub      = inp.x$compound_id,
               mass_sub = inp.x$mass_measured,
               rt_sub   = inp.x$retention_time),
    by    = "sub",
    all.x = TRUE
  )
  
  pairs <- merge(
    pairs,
    data.table(prod      = inp.x$compound_id,
               mass_prod = inp.x$mass_measured,
               rt_prod   = inp.x$retention_time),
    by    = "prod",
    all.x = TRUE
  )
  
  # Coerce to atomic before comparison
  pairs[, sub  := as.character(sub)]
  pairs[, prod := as.character(prod)]
  
  # avoid self matches
  pairs <- pairs[sub != prod]
  
  # ---- MASS FILTER ----
  pairs[, diff := (mass_prod - mass_sub) - mzdiff]
  pairs <- pairs[abs(diff) <= mzerr]
  
  if (nrow(pairs) == 0) return(data.frame())
  
  # ---- RT FILTER ----
  if (direc == 1) {
    pairs <- pairs[rt_prod < rt_sub - peakwidth]
  } else if (direc == 2) {
    pairs <- pairs[rt_prod > rt_sub + peakwidth]
  } else {
    pairs <- pairs[abs(rt_prod - rt_sub) > peakwidth]
  }
  
  if (nrow(pairs) == 0) return(data.frame())
  
  # ---- MS2 SCORING ----
  out_list <- vector("list", nrow(pairs))
  
  for (i in seq_len(nrow(pairs))) {
    
    res <- targMS2comp_SQL(
      pairs$sub[i],
      pairs$prod[i],
      ms2_split,
      IntThres = IntThres,
      mz_tol = mz_tol
    )
    
    if (is.null(res) || !is.data.frame(res)) {
      res <- data.frame(
        COMPID.sub = pairs$sub[i],
        COMPID.prod = pairs$prod[i],
        DOT_IONS = NA,
        DOT_LOSS = NA,
        COMMON_IONS = NA,
        COMMON_LOSS = NA,
        FORW_IONS = NA,
        REV_IONS = NA,
        FORW_LOSS = NA,
        REV_LOSS = NA
      )
    }
    
    out_list[[i]] <- res
  }
  
  out <- rbindlist(out_list, fill = TRUE)
  
  if (nrow(out) == 0) return(out)
  
  out <- out[
    DOT_IONS >= dot_thresh &
      COMMON_IONS >= min_common
  ]
  
  return(out)
}