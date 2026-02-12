#' @title Compute CSPP scores for one experiment and update the compound_add SQL
#'        table
#'
#' @description 
#' cspp.tot_SQL() computes CSPP (Conversion-Specific Product Pair) scores for all
#' compounds belonging to a given experiment and updates the corresponding rows
#' in the \code{compound_add} table of a SQLite database.
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
#' @param expid `Integer(1)` experiment id. Only compounds belonging to this
#'   experiment are processed.
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
cspp.tot_SQL <- function(sql_path, expid, mzerr = 0.015,
                         cspp = "cspp.txt", peakwidth = NULL, IntThres  = 100) {
  

  data_con <- dbConnect(RSQLite::SQLite(), sql_path)
  on.exit(dbDisconnect(data_con), add = TRUE)
  
  if (is.null(peakwidth)) peakwidth <- 0.2

  inp.x <- dbGetQuery(
    data_con,
    sprintf(
      "SELECT compound_id,
              mass_measured,
              retention_time
       FROM ms_compound
       WHERE expid = %d",
      expid
    )
  )
  
  if (nrow(inp.x) == 0) {
    stop("No compounds found for expid = ", expid)
  }

  comp_add <- dbGetQuery(
    data_con,
    "SELECT * FROM compound_add"
  )
  

  conv.table <- read.table(
    cspp,
    header = TRUE,
    sep = "\t",
    stringsAsFactors = FALSE
  )
  
  conver <- data.frame(
    conv.type  = as.character(conv.table[, 1]),
    conv.mz    = as.numeric(conv.table[, 3]),
    conv.direc = as.integer(conv.table[, 5]),
    conv.colmn = as.integer(conv.table[, 6]),
    stringsAsFactors = FALSE
  )
  

  for (k in seq_len(nrow(conver))) {
    
    cspp.df <- conv.CSPP_SQL(
      inp.x,
      mzdiff    = conver$conv.mz[k],
      direc     = conver$conv.direc[k],
      peakwidth = peakwidth,
      mzerr     = mzerr
    )
    
    comp_add <- rank.cspp(
      cspp.df,
      conver$conv.colmn[k],
      comp_add
    )
  }

  comp_add$compound_id <- inp.x$compound_id[
    seq_len(nrow(comp_add))
  ]
  

  cols_to_update <- setdiff(colnames(comp_add), "compound_id")
  
  dbBegin(data_con)
  
  for (i in seq_len(nrow(comp_add))) {
    
    cid <- comp_add$compound_id[i]
    
    for (col in cols_to_update) {
      
      dbExecute(
        data_con,
        sprintf(
          "UPDATE compound_add
           SET %s = ?
           WHERE compound_id = ?",
          col
        ),
        params = list(comp_add[i, col], cid)
      )
    }
  }
  
  dbCommit(data_con)
  

  invisible(comp_add)
}


# cspp_add<-cspp.tot(base.dir,finlist,SubDB="FTneg",Prod.exp=2)
# write.table(cspp_add,"compound_add.txt",sep="\t",row.names=F)
