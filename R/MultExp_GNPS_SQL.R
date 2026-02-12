#' @title feature pairs generation based on highly similar CID spectra.
#'
#' @description 
#' MultExp_GNPS_SQL() function generates feature pairs based on highly similar 
#' CID spectra followed by mass difference computation. This method is similar 
#' to the method used by the Global Natural Products Social Molecular Networking 
#' (GNPS) website (Wang et al., 2016). GNPS-like associations are computed and 
#' stored within the ‘gnps_add’ table in the sql database.
#'
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
#' @param peakwidth `Numeric(1)` half of the retention time window that is 
#' defined for the CSPP ‘substrate’ feature allowing
#' a minimum retention time difference with the CSPP ‘product’ feature,
#'  by default set on \code{0.2}.
#'
#' @param min minimum absolute intensity of product ion to be withheld for CID 
#' spectral matching, by default \code{5} for QTOF, but the user can set it to
#'  \code{100} for FTMS spectra.
#' 
#' @param adduct mass of buffer, by default \code{46.0055 Da} representing 
#'  formic acid. 
#' @param thr1 average of the number of common ions and the number of common 
#' losses, set on \code{3} by default.
#' 
#' @param thr2 ratio of thr1 to the number of product ions in the CID spectrum 
#' of the feature in the feature pair having the lowest number of product ions,
#' by default set on \code{0.1}.
#' 
#' @param thr3 average of the dot products obtained for the common product ions 
#' and the common neutral losses, by default set on \code{0.4}.
#' 
#'@return Invisibly returns \code{TRUE}. The primary effect of the
#' function is the update of the \code{gnps_add} SQL table.
#'
#' @author Ahlam Mentag
#' 
#' @export
MultExp_GNPS_SQL <- function(sql_path, startexp, stopexp, peakwidth = NULL,
                             mzerr = 0.015, min = 5, adduct = NULL,
                             thr1 = NULL, thr2 = NULL, thr3 = NULL) {
  
  if (is.null(peakwidth)) peakwidth <- 0.2
  if (is.null(adduct)) adduct <- 46.0055
  if (is.null(thr1)) thr1 <- 3
  if (is.null(thr2)) thr2 <- 0.1
  if (is.null(thr3)) thr3 <- 0.4
  
  for (expid in seq(startexp, stopexp)) {
    
    message("Running GNPS for experiment ", expid)
    
    conv.GNPS_SQL(
      sql_path  = sql_path,
      expid     = expid,
      peakwidth = peakwidth,
      mzerr     = mzerr,
      min       = min,
      adduct    = adduct,
      thr1      = thr1,
      thr2      = thr2,
      thr3      = thr3,
    )
  }
  
  invisible(TRUE)
}


# gnps_add<-MultExp_GNPS(base.dir,finlist,SubDB="FTneg",startexp=1,stopexp=3)
# write.table(gnps_add,"gnps_add.txt",sep="\t",row.names=F)