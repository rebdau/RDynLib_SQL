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
#' @import DBI
#' @import RSQLite
#' @import data.table
#' 
#' @author Ahlam Mentag
#' 
#' @export
cspp.tot_SQL <- function(sql_path, expid, mzerr = 0.015,
                         cspp = "cspp.txt",
                         peakwidth = NULL,
                         IntThres = 100) {
  
  if (is.null(peakwidth)) peakwidth <- 0.2
  
  data_con <- dbConnect(SQLite(), sql_path)
  on.exit(dbDisconnect(data_con), add = TRUE)
  
  # Load compounds
  inp.x <- dbGetQuery(
    data_con,
    sprintf(
      "SELECT compound_id, mass_measured, retention_time
       FROM ms_compound
       WHERE expid = %d", expid
    )
  )
  
  if (nrow(inp.x) == 0)
    stop("No compounds found for expid = ", expid)
  
  inp.x <- inp.x[order(inp.x$mass_measured), ]
  
  # Load compound_add
  comp_add <- dbReadTable(data_con, "compound_add")
  
  # Load conversion table
  conv.table <- fread(cspp, header = TRUE, sep = "\t")
  
  conver <- data.frame(
    conv.type  = as.character(conv.table[[1]]),
    conv.mz    = as.numeric(conv.table[[3]]),
    conv.direc = as.integer(conv.table[[5]]),
    conv.colmn = as.integer(conv.table[[6]]),
    stringsAsFactors = FALSE
  )
  
  # LOAD ALL MS2 ONCE
  message("Loading all MS2 spectra into memory...")
  
  spec_df <- dbGetQuery(
    data_con,
    "
    SELECT s.spectrum_id,
           s.compound_id,
           s.precursor_mz
    FROM msms_spectrum s
    WHERE s.ms_level = 2
      AND s.spectrum_id = (
          SELECT MIN(s2.spectrum_id)
          FROM msms_spectrum s2
          WHERE s2.compound_id = s.compound_id
            AND s2.ms_level = 2
      )
    "
  )
  
  if (nrow(spec_df) == 0)
    stop("No MS2 spectra found.")
  
  # peak_df <- dbGetQuery(
  #   data_con,
  #   "SELECT spectrum_id, mz, intensity FROM msms_spectrum_peak"
  # )
  
  
  peak_df <- DBI::dbGetQuery(
    data_con,
    sprintf("
      SELECT p.spectrum_id, p.mz, p.intensity
      FROM msms_spectrum_peak p
      INNER JOIN msms_spectrum s
        ON p.spectrum_id = s.spectrum_id
      INNER JOIN ms_compound c
        ON s.compound_id = c.compound_id
      WHERE c.expid = %d
    ", expid)
  )
  
  ms2_df <- merge(peak_df, spec_df, by = "spectrum_id")
  ms2_df <- ms2_df[ms2_df$intensity >= IntThres, ]
  
  ms2_split <- split(ms2_df, ms2_df$compound_id)
  
  message("Loaded MS2 for ", length(ms2_split), " compounds.")
  
  
  for (k in seq_len(nrow(conver))) {
    
    message("Processing conversion: ", k, "/", nrow(conver))
    
    cspp.df <- conv.CSPP_SQL(
      inp.x,
      mzdiff    = conver$conv.mz[k],
      direc     = conver$conv.direc[k],
      peakwidth = peakwidth,
      mzerr     = mzerr,
      ms2_split = ms2_split
    )
    
    
    if (is.null(cspp.df) || nrow(cspp.df) == 0) {
      message("  No matches found for conversion ", conver$conv.type[k])
      next
    }
    
    if (ncol(cspp.df) < 14) {
      warning("  Unexpected column structure for conversion ",
              conver$conv.type[k], ". Skipping.")
      next
    }
    
    comp_add <- rank.cspp(cspp.df,
                          conver$conv.colmn[k],
                          comp_add)
  }
  
  # Ensure correct compound_id alignment
  if (nrow(comp_add) <= nrow(inp.x)) {
    comp_add$compound_id <- inp.x$compound_id[seq_len(nrow(comp_add))]
  }
  
  # Rewrite table
  message("Writing results to database...")
  
  dbExecute(data_con, "DROP TABLE IF EXISTS compound_add_new")
  dbWriteTable(data_con, "compound_add_new", comp_add)
  
  dbExecute(data_con, "DROP TABLE compound_add")
  dbExecute(data_con, "ALTER TABLE compound_add_new RENAME TO compound_add")
  
  message("Done.")
  
  invisible(comp_add)
}



# cspp_add<-cspp.tot(base.dir,finlist,SubDB="FTneg",Prod.exp=2)
# write.table(cspp_add,"compound_add.txt",sep="\t",row.names=F)
