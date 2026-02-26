#' @title Convert and rank GNPS-like MS/MS relationships and store results
#'  in SQLite
#'
#' @description
#' 
#' This function iterates the GNPS-like feature pair generation procedure over
#' all compound_id belonging to an experiment. the features (cfr. compound_id)   
#' are first sorted with increasing retention times. 
#' For each feature, CID spectral similarities are computed with all features 
#' that elute later. Only those feature pairs for which the CID spectral
#' similarity passes certain thresholds (thr1, thr2 and thr3) are retained. 
#' Subsequently, for a particular feature, all generated feature pairs are
#' ranked, via the rank.GNPS() function.
#'
#' @param sql_path `Character(1)` string giving the path to the SQLite database.
#' 
#' @param expid Integer. Experiment ID used to select compounds from the database.
#' 
#' @param peakwidth `Numeric(1)` half of the retention time window that is 
#' defined for the CSPP ‘substrate’ feature allowing
#' a minimum retention time difference with the CSPP ‘product’ feature,
#' by default set on \code{0.2}. 
#' 
#' @param mzerr `Numeric(1)` mass tolerance (in Da) used to match precursor and
#'   product ions. Default it is configured for ftms data
#'   \code{0.01}, for qtof data the user could set it to \code{0.015}.
#' 
#' @param min minimum absolute intensity of product ion to be withheld for CID 
#' spectral matching, by default \code{5} for QTOF, but the user can set it to
#'  \code{100} for FTMS spectra.
#'  
#' @param adduct mass of buffer, by default \code{46.0055 Da} representing 
#'  formic acid. 
#'  
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
#' @param IntThres Numeric intensity threshold applied to MS/MS fragment ions
#'   before similarity calculations. Default it is configured for ftms data
#'   \code{100}, for qtof data the user could set it to \code{5}.
#'    
#' @return
#' Invisibly returns a data frame corresponding to the updated GNPS SQL table
#' for the processed experiment.
#'
#' @author Ahlam Mentag
#' 
#' @export
conv.GNPS_SQL <- function(sql_path, expid,
                          peakwidth = 0.2,
                          mzerr = 0.015,
                          min = 5,
                          adduct = 46.0055,
                          thr1 = 3,
                          thr2 = 0.1,
                          thr3 = 0.4,
                          IntThres = 5) {
  
  con <- DBI::dbConnect(RSQLite::SQLite(), sql_path)
  on.exit(DBI::dbDisconnect(con), add = TRUE)
  
  # Load compounds 
  inp.x <- DBI::dbGetQuery(con, sprintf("
      SELECT compound_id, mass_measured, retention_time
      FROM ms_compound
      WHERE expid = %d", expid))
  if (nrow(inp.x) == 0) stop("No compounds found for expid = ", expid)
  
  inp.x$compound_id    <- as.integer(inp.x$compound_id)
  inp.x$mass_measured  <- as.numeric(inp.x$mass_measured)
  inp.x$retention_time <- as.numeric(inp.x$retention_time)
  inp.x <- inp.x[order(inp.x$retention_time), ]
  
  #Load MS2 spectra 
  # spec_df <- DBI::dbGetQuery(con, sprintf("
  #     SELECT s.spectrum_id, s.compound_id, s.precursor_mz
  #     FROM msms_spectrum s
  #     INNER JOIN ms_compound c ON s.compound_id = c.compound_id
  #     WHERE s.ms_level = 2 AND c.expid = %d
  #       AND s.spectrum_id = (SELECT MIN(s2.spectrum_id)
  #                            FROM msms_spectrum s2
  #                            WHERE s2.compound_id = s.compound_id
  #                              AND s2.ms_level = 2)
  # ", expid))
  
  spec_df <- DBI::dbGetQuery(con, sprintf("
    SELECT s.spectrum_id, s.compound_id, s.precursor_mz
    FROM msms_spectrum s
    INNER JOIN ms_compound c ON s.compound_id = c.compound_id
    WHERE s.ms_level = 2 
      AND c.expid = %d
", expid))
  
  peak_df <- DBI::dbGetQuery(con, sprintf("
      SELECT p.spectrum_id, p.mz, p.intensity
      FROM msms_spectrum_peak p
      INNER JOIN msms_spectrum s ON p.spectrum_id = s.spectrum_id
      INNER JOIN ms_compound c ON s.compound_id = c.compound_id
      WHERE c.expid = %d
  ", expid))
  
  ms2_df <- merge(peak_df, spec_df, by = "spectrum_id")
  ms2_df$mz           <- as.numeric(ms2_df$mz)
  ms2_df$intensity    <- as.numeric(ms2_df$intensity)
  ms2_df$precursor_mz <- as.numeric(ms2_df$precursor_mz)
  
  # Apply intensity threshold 
  ms2_df <- ms2_df[ms2_df$intensity >= IntThres, ]
  
  #Split MS2 by compound
  ms2_split <- split(ms2_df, ms2_df$compound_id)
  ms2_split <- lapply(ms2_split, function(x) if (nrow(x) == 0) NULL else x)
  
  gnps_add <- tryCatch({
    DBI::dbReadTable(con, "gnps_add")
  }, error = function(e) {
    data.frame(compound_id = inp.x$compound_id,
               GNPS_1 = NA_character_,
               GNPS_2 = NA_character_,
               GNPS_3 = NA_character_,
               GNPS_4 = NA_character_,
               GNPS_5 = NA_character_,
               stringsAsFactors = FALSE)
  })
  
  missing_ids <- setdiff(inp.x$compound_id, gnps_add$compound_id)
  if (length(missing_ids) > 0) {
    gnps_add <- rbind(gnps_add,
                      data.frame(compound_id = missing_ids,
                                 GNPS_1 = NA_character_,
                                 GNPS_2 = NA_character_,
                                 GNPS_3 = NA_character_,
                                 GNPS_4 = NA_character_,
                                 GNPS_5 = NA_character_,
                                 stringsAsFactors = FALSE))
  }
  
  #  Main loop over compounds 
  for (i in seq_len(nrow(inp.x))) {
    sub_id <- inp.x$compound_id[i]
    sub_mz <- inp.x$mass_measured[i]
    sub_rt <- inp.x$retention_time[i] + peakwidth
    
    sub_lowmz  <- sub_mz - mzerr
    sub_highmz <- sub_mz + 1.0034 + mzerr
    forb_lowmz  <- sub_mz + adduct - mzerr
    forb_highmz <- sub_mz + adduct + 1.0034 + mzerr
    
    # Candidate products with higher RT
    valid_j <- which(inp.x$retention_time > sub_rt)
    if (length(valid_j) == 0) next
    
    gnps_list <- list()
    k <- 1
    
    for (j in valid_j) {
      prod_id <- inp.x$compound_id[j]
      mzj <- inp.x$mass_measured[j]
      
      # Exclude isotopes/adducts
      if ((mzj < sub_lowmz || mzj > sub_highmz) &&
          (mzj < forb_lowmz || mzj > forb_highmz)) {
        
        # Skip if MS2 missing 
        if (is.null(ms2_split[[as.character(prod_id)]]) ||
            is.null(ms2_split[[as.character(sub_id)]])) next
        
        out <- targMS2comp_SQL(sub_id, prod_id, ms2_split, IntThres = min)
        if (is.null(out) || nrow(out) == 0) next
        
        # thresholds
        min_ions <- min(out$IONS.sub, out$IONS.prod)
        if (is.na(min_ions) || min_ions == 0) next
        
        thresh1 <- (out$COMMON_IONS + out$COMMON_LOSS)/2
        thresh2 <- thresh1 / min_ions
        thresh3 <- (out$DOT_IONS + out$DOT_LOSS)/2
        
        if (!is.na(thresh1) && !is.na(thresh2) && !is.na(thresh3) &&
            thresh1 > thr1 && thresh2 > thr2 && thresh3 > thr3) {
          gnps_list[[k]] <- out
          k <- k + 1
        }
      }
    }
    
    if (length(gnps_list) > 0) {
      gnps.df <- do.call(rbind, gnps_list)
      gnps_add <- rank.GNPS_SQL(gnps.df, gnps_add, inp.x)
    }
  }
  
  # Write back to database 
  DBI::dbExecute(con, "DROP TABLE IF EXISTS gnps_add_new")
  DBI::dbWriteTable(con, "gnps_add_new", gnps_add)
  DBI::dbExecute(con, "DROP TABLE IF EXISTS gnps_add")
  DBI::dbExecute(con, "ALTER TABLE gnps_add_new RENAME TO gnps_add")
  
  return(gnps_add)
}