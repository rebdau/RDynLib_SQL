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
cspp.tot_SQL <- function(sql_path, expid,
                         mzerr = 1,
                         cspp = "cspp.txt",
                         peakwidth = 0.2,
                         IntThres = 0,
                         mz_tol = 0.3,
                         dot_thresh = 0,
                         min_common = 0,
                         spectrum_type = "assembled",
                         charge_file = NULL,
                         reset = TRUE,        
                         verbose = TRUE) {
  

  
  log <- function(...) if (verbose) message(sprintf(...))
  
  con <- dbConnect(SQLite(), sql_path)
  on.exit(dbDisconnect(con), add = TRUE)
  
  
  inp.x <- as.data.table(dbGetQuery(con, sprintf("
    SELECT compound_id, mass_measured, retention_time
    FROM ms_compound
    WHERE expid = %d
  ", expid)))
  
  inp.x <- inp.x[!is.na(mass_measured)]
  inp.x[, compound_id := as.character(compound_id)]
  
  log("Loaded MS1 compounds: %d", nrow(inp.x))
  
  

  if (!is.null(charge_file)) {
    
    charge_dt <- as.data.table(readxl::read_excel(charge_file))
    charge_dt[, compound_id := as.character(compound_id)]
    charge_dt[, doublecharged := tolower(trimws(as.character(doublecharged)))]
    
    inp.x <- merge(
      inp.x,
      charge_dt[, .(compound_id, doublecharged)],
      by = "compound_id",
      all.x = TRUE
    )
    
    inp.x[, mass_measured := fifelse(
      doublecharged %in% c("true", "t", "1"),
      (mass_measured * 2) + 1.007276,
      mass_measured
    )]
    
    inp.x[, doublecharged := NULL]
    
    log("Charge correction applied from Excel file.")
  }

  
  
  conv_raw <- fread(cspp)
  
  conv <- data.table(
    type = trimws(conv_raw[[1]]),
    mzdiff = as.numeric(conv_raw[[3]]),
    direc = as.integer(conv_raw[[5]])
  )
  

  conv[, col_name := type]
  
  
  spec_df <- as.data.table(dbGetQuery(con, sprintf("
    SELECT spectrum_id, compound_id, precursor_mz
    FROM msms_spectrum
    WHERE ms_level = 2
      AND spectrum_type = '%s'
  ", spectrum_type)))
  
  peak_df <- as.data.table(dbGetQuery(con, sprintf("
    SELECT p.spectrum_id, p.mz, p.intensity, s.compound_id
    FROM msms_spectrum_peak p
    JOIN msms_spectrum s USING(spectrum_id)
    WHERE p.intensity >= %d
  ", IntThres)))
  
  ms2_df <- peak_df[spec_df, on = "spectrum_id", nomatch = 0]
  ms2_list <- split(ms2_df, ms2_df$compound_id)
  
  log("Loaded MS2 compounds: %d", length(ms2_list))
  
  
  comp_add <- as.data.table(dbReadTable(con, "compound_add"))
  comp_add[, compound_id := as.character(compound_id)]
  

  missing_ids <- setdiff(inp.x$compound_id, comp_add$compound_id)
  
  if (length(missing_ids) > 0) {
    
    log("Adding %d missing compound_ids to compound_add",
        length(missing_ids))
    
    new_rows <- data.table(compound_id = missing_ids)
    
    comp_add <- rbind(
      comp_add,
      new_rows,
      fill = TRUE
    )
  }

  
 
  cspp_cols <- unique(conv$col_name)
  missing_cols <- setdiff(cspp_cols, names(comp_add))
  
  if (length(missing_cols) > 0) {
    log("Adding missing columns to compound_add: %s",
        paste(missing_cols, collapse = ", "))
    
    comp_add[, (missing_cols) := NA_character_]
  }

  
  
  if (reset) {
    for (col in conv$col_name) {
      if (col %in% names(comp_add)) {
        comp_add[, (col) := NA_character_]
      }
    }
  }
  
  success <- character()
  
  
  for (k in seq_len(nrow(conv))) {
    
    log("Conversion %d/%d (%s)", k, nrow(conv), conv$type[k])
    
    res <- conv.CSPP_SQL(
      inp.x,
      mzdiff = conv$mzdiff[k],
      direc = conv$direc[k],
      peakwidth = peakwidth,
      mzerr = mzerr,
      ms2_split = ms2_list,
      IntThres = IntThres,
      mz_tol = mz_tol,
      dot_thresh = dot_thresh,
      min_common = min_common
    )
    
    if (is.null(res) || nrow(res) == 0) next
    
    res <- as.data.table(res)
    
    res <- unique(res, by = c("COMPID.sub", "COMPID.prod"))
    
    res[, val := paste0(
      "!!", COMMON_IONS,
      "!", round((FORW_IONS + REV_IONS)/2, 3),
      "!", round((DOT_IONS + DOT_LOSS)/2, 3),
      "!!", COMPID.prod
    )]
    
    res_agg <- res[, .(val = paste(unique(val), collapse = "|")), by = COMPID.sub]
    res_agg[, sub_id := as.character(COMPID.sub)]
    
    col <- conv$col_name[k]
    if (!(col %in% names(comp_add))) next
    
    comp_add[res_agg,
             on = .(compound_id = sub_id),
             (col) := {
               
               old <- get(col)
               
               combined <- ifelse(
                 is.na(old) | old == "",
                 val,
                 paste0(old, "|", val)
               )
               
               # remove duplicates safely
               sapply(strsplit(combined, "\\|"), function(x) {
                 paste(unique(x), collapse = "|")
               })
             }
    ]
    
    
    success <- c(success, col)
  }
  
  dbWriteTable(con, "compound_add", comp_add, overwrite = TRUE)
  
  log("Successful conversions: %s",
      if (length(success) == 0) "none" else paste(unique(success), collapse = ", "))
  
  invisible(comp_add)
}