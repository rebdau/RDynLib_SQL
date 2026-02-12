#' @title Compute CSPP scores for multiple experiment and update the compound_add SQL
#'        table
#'
#' @description 
#' MultExp_CSPP_SQL() computes CSPP (Conversion-Specific Product Pair) scores for 
#' all compounds belonging to a set of experiments and updates the corresponding
#' rows in the \code{compound_add} table of a SQLite database.
#'
#' @details 
#' For each compound in the experiment, theoretical conversion products are
#' generated based on the conversion definitions provided in a CSPP configuration
#' file. Candidate precursor/product pairs are matched using mass differences
#' and retention time constraints, and their MS/MS spectra are compared to
#' compute similarity scores.
#'
#' No rows representing the different compound IDs of the experiment need to be
#' pre-created in the \code{compound_add} table. During execution, CSPP scores
#' are written to the columns specified in the CSPP configuration file. If no
#' valid conversions are found for a given compound, the corresponding CSPP
#' fields are left as \code{NA}.
#'
#' @param sql_path `Character(1)` string giving the path to the SQLite database.
#'
#' @param startexp `Integer(1)` the start experiment id. 
#'
#' @param stopexp `Integer(1)` the last experiment id, where the process should
#'        stop. 
#'
#' @param mzerr `Numeric(1)` mass tolerance (in Da) used to match precursor and
#'   product ions. Default it is configured for ftms data
#'   \code{0.01}, for qtof data the user could set it to \code{0.015}.
#'
#' @param cspp `Character(1)` string giving the path to the CSPP configuration 
#'   file. This file defines the conversion types, mass differences,
#'   elution order, and the target columns in \code{compound_add}.
#'
#' @param peakwidth `Numeric(1)` retention time window used to enforce elution 
#'   order constraints between substrate and product compounds. If \code{NULL},
#'   defaults to \code{0.2}.
#'
#' @param IntThres Numeric intensity threshold applied to MS/MS fragment ions
#'   before similarity calculations. Default it is configured for ftms data
#'   \code{100}, for qtof data the user could set it to \code{5}.
#'
#' @return Invisibly returns the updated \code{compound_add} data frame
#'   corresponding to the processed experiment. The primary effect of the
#'   function is the update of the SQL table.
#'
#' @author Ahlam Mentag
#' 
#' @export
MultExp_CSPP_SQL <- function(sql_path,
                             startexp,
                             stopexp,
                             peakwidth = NULL,
                             mzerr = 0.015,
                             IntThres  = 100,
                             cspp = "cspp.txt") {
  

  if (is.null(peakwidth)) peakwidth <- 0.2
  
 
  ## loop over experiments
  for (expid in seq(startexp, stopexp)) {
    
    message("Running CSPP for experiment ", expid)
    
    cspp.tot_SQL(
      sql_path  = sql_path,
      expid     = expid,
      mzerr     = mzerr,
      cspp      = cspp,
      peakwidth = peakwidth
    )
  }
  
  invisible(TRUE)
}
