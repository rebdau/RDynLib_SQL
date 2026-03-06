#' @title Merge two SQLite spectral databases.
#'
#' @description
#' Merges two SQLite databases containing mass spectrometry data into a
#' single output database. The function aligns tables present in both
#' databases, shifts identifier columns to avoid ID conflicts,
#' and combines the records. Tables that exist in only
#' one database are copied directly to the output database.
#'
#' Identifier columns (`expid`, `compound_id`, `spectrum_id`, `peak_id`)
#' are automatically shifted to ensure uniqueness across the merged data.
#' Retention times (`rtime`, `retention_time`) from the additional database
#' are converted from seconds to minutes before merging.
#'
#' @param main_db 'Character(1)' Path to the primary SQLite database.
#' 
#' @param add_db 'Character(1)' Path to the second SQLite database that would be 
#' added to the main database.
#' 
#' @param output_db Character. Path where the merged SQLite database will
#'   be written.
#'
#' @return No value is returned. The function writes the merged database
#'   to `output_db` and prints progress messages during execution.
#'
#' @author Ahlam Mentag
#' @export
merge_databases <- function(main_db, add_db, output_db) {
  
  # Connect to databases
  con_main <- dbConnect(SQLite(), main_db)
  con_add  <- dbConnect(SQLite(), add_db)
  con_out  <- dbConnect(SQLite(), output_db)
  
  tables_main <- setdiff(dbListTables(con_main), "sqlite_sequence")
  tables_add  <- setdiff(dbListTables(con_add), "sqlite_sequence")
  
  # Intersect the columns of both databases
  common_tables <- intersect(tolower(tables_main), tolower(tables_add))
  
  # ID shifting function
  shift_ids <- function(df, exp_shift, spectrum_shift, peak_shift, compound_shift) {
    if ("expid" %in% names(df))        df$expid       <- as.integer(df$expid) + exp_shift
    if ("spectrum_id" %in% names(df))  df$spectrum_id <- as.integer(df$spectrum_id) + spectrum_shift
    if ("peak_id" %in% names(df))      df$peak_id     <- as.integer(df$peak_id) + peak_shift
    if ("compound_id" %in% names(df))  df$compound_id <- as.character(as.integer(df$compound_id) + compound_shift)
    if ("subsid" %in% names(df))       df$subsid      <- as.integer(df$subsid)
    return(df)
  }
  
  # Compute ID shift
  compute_shift <- function(con_main, con_add, table, key) {
    if (!(table %in% tolower(tables_main)) || !(table %in% tolower(tables_add))) return(0)
    main_tbl <- dbReadTable(con_main, tables_main[tolower(tables_main) == table])
    add_tbl  <- dbReadTable(con_add, tables_add[tolower(tables_add) == table])
    if (!key %in% names(main_tbl) || !key %in% names(add_tbl)) return(0)
    main_vals <- suppressWarnings(as.integer(main_tbl[[key]]))
    add_vals  <- suppressWarnings(as.integer(add_tbl[[key]]))
    max_main <- if(any(!is.na(main_vals))) max(main_vals, na.rm=TRUE) else 0
    min_add  <- if(any(!is.na(add_vals)))  min(add_vals, na.rm=TRUE) else 0
    return(max_main + 1 - min_add)
  }
  
  exp_shift      <- compute_shift(con_main, con_add, "experiment", "expid")
  compound_shift <- compute_shift(con_main, con_add, "ms_compound", "compound_id")
  spectrum_shift <- compute_shift(con_main, con_add, "msms_spectrum", "spectrum_id")
  peak_shift     <- compute_shift(con_main, con_add, "msms_spectrum_peak", "peak_id")
  
  # Fix column types
  fix_column_types <- function(df, table_name) {
    if (nrow(df) == 0) return(df)
    
    if (table_name == "msms_spectrum") {
      df$spectrum_id <- as.integer(df$spectrum_id)
      df$compound_id <- as.character(df$compound_id)
      df$ms_level <- as.integer(df$ms_level)
      df$polarity <- as.integer(df$polarity)
      df$spectrum_type <- as.character(df$spectrum_type)
      df$precursor_mz <- as.numeric(df$precursor_mz)
      df$precursorIntensity <- as.numeric(df$precursorIntensity)
      df$precursorCharge <- as.integer(df$precursorCharge)
      df$collision_energy <- as.character(df$collision_energy)
      df$isolationWindowLowerMz <- as.numeric(df$isolationWindowLowerMz)
      df$isolationWindowTargetMz <- as.numeric(df$isolationWindowTargetMz)
      df$isolationWindowUpperMz <- as.numeric(df$isolationWindowUpperMz)
      df$peaks_count <- as.integer(df$peaks_count)
      df$totIonCurrent <- as.numeric(df$totIonCurrent)
      df$basePeakMZ <- as.numeric(df$basePeakMZ)
      df$basePeakIntensity <- as.numeric(df$basePeakIntensity)
      df$ionisationEnergy <- as.numeric(df$ionisationEnergy)
      df$lowMZ <- as.numeric(df$lowMZ)
      df$highMZ <- as.numeric(df$highMZ)
      df$mergedScan <- as.integer(df$mergedScan)
      df$mergedResultScanNum <- as.integer(df$mergedResultScanNum)
      df$mergedResultStartScanNum <- as.integer(df$mergedResultStartScanNum)
      df$mergedResultEndScanNum <- as.integer(df$mergedResultEndScanNum)
      df$injectionTime <- as.numeric(df$injectionTime)
      df$filterString <- as.character(df$filterString)
      df$spectrumId <- as.integer(df$spectrumId)
      df$ionMobilityDriftTime <- as.numeric(df$ionMobilityDriftTime)
      df$scanWindowLowerLimit <- as.numeric(df$scanWindowLowerLimit)
      df$scanWindowUpperLimit <- as.numeric(df$scanWindowUpperLimit)
      df$electronBeamEnergy <- as.numeric(df$electronBeamEnergy)
      df$originalPrecursorMz <- as.numeric(df$originalPrecursorMz)
      df$precursorPurity <- as.numeric(df$precursorPurity)
      df$chromPeakRT <- as.numeric(df$chromPeakRT)
      df$chromPeakMz <- as.numeric(df$chromPeakMz)
      df$chromPeakId <- as.character(df$chromPeakId)
      df$rtime <- as.numeric(df$rtime)
      df$scanIndex <- as.integer(df$scanIndex)
      df$dataStorage <- as.character(df$dataStorage)
      df$centroided <- as.integer(df$centroided)
      df$smoothed <- as.integer(df$smoothed)
      df$instrument <- as.character(df$instrument)
      df$splash <- as.character(df$splash)
      df$instrument_type <- as.character(df$instrument_type)
      df$acquisitionNum <- as.integer(df$acquisitionNum)
      df$precScanNum <- as.integer(df$precScanNum)
      df$predicted <- as.numeric(df$predicted)
      df$dataOrigin <- as.character(df$dataOrigin)
      df$original_id <- as.character(df$original_id)
      
    } else if (table_name == "ms_compound") {
      df$compound_id <- as.character(df$compound_id)
      df$expid <- as.integer(df$expid)
      df$exactmass <- as.numeric(df$exactmass)
      df$subsid <- as.integer(df$subsid)
      df$retention_time <- as.numeric(df$retention_time)
      df$mass_measured <- as.numeric(df$mass_measured)
      df$wavelen <- as.numeric(df$wavelen)
      df$isotope_ratio <- as.numeric(df$isotope_ratio)
      df$drift_time <- as.numeric(df$drift_time)
    } else if (table_name == "msms_spectrum_peak") {
      df$peak_id <- as.integer(df$peak_id)
      df$spectrum_id <- as.integer(df$spectrum_id)
      df$mz <- as.numeric(df$mz)
      df$intensity <- as.numeric(df$intensity)
    }
    return(df)
  }
  
  # Merge common tables
  for(tbl_name in common_tables) {
    tbl_main <- dbReadTable(con_main, tables_main[tolower(tables_main) == tbl_name])
    tbl_add  <- dbReadTable(con_add, tables_add[tolower(tables_add) == tbl_name])
    
    if(nrow(tbl_add) > 0) {
      tbl_add <- shift_ids(tbl_add, exp_shift, spectrum_shift, peak_shift, compound_shift)
      
      # Convert times to minutes in add_db
      if ("rtime" %in% names(tbl_add)) {
        tbl_add$rtime <- as.numeric(tbl_add$rtime) / 60
      }
      if ("retention_time" %in% names(tbl_add)) {
        tbl_add$retention_time <- as.numeric(tbl_add$retention_time) / 60
      }
    }
    
    all_cols <- union(names(tbl_main), names(tbl_add))
    for(col in setdiff(all_cols, names(tbl_main))) tbl_main[[col]] <- rep(NA, nrow(tbl_main))
    for(col in setdiff(all_cols, names(tbl_add))) tbl_add[[col]] <- rep(NA, nrow(tbl_add))
    
    tbl_main <- tbl_main[, all_cols]
    tbl_add  <- tbl_add[, all_cols]
    
    tbl_main[] <- lapply(tbl_main, as.character)
    tbl_add[]  <- lapply(tbl_add, as.character)
    
    merged_tbl <- bind_rows(tbl_main, tbl_add)
    merged_tbl <- fix_column_types(merged_tbl, tbl_name)
    
    dbWriteTable(con_out, tbl_name, merged_tbl, overwrite = TRUE)
    cat("Merged:", tbl_name, "\n")
  }
  
  # Copy unique tables
  for(tbl_name in setdiff(tables_main, tables_add)) {
    tbl <- dbReadTable(con_main, tbl_name)
    tbl <- fix_column_types(tbl, tbl_name)
    dbWriteTable(con_out, tbl_name, tbl, overwrite = TRUE)
  }
  
  for(tbl_name in setdiff(tables_add, tables_main)) {
    tbl <- dbReadTable(con_add, tbl_name)
    if(nrow(tbl) > 0) tbl <- shift_ids(tbl, exp_shift, spectrum_shift, peak_shift, compound_shift)
    
    # Convert times to minutes in add_db
    if ("rtime" %in% names(tbl)) {
      tbl$rtime <- as.numeric(tbl$rtime) / 60
    }
    if ("retention_time" %in% names(tbl)) {
      tbl$retention_time <- as.numeric(tbl$retention_time) / 60
    }
    
    tbl <- fix_column_types(tbl, tbl_name)
    dbWriteTable(con_out, tbl_name, tbl, overwrite = TRUE)
  }
  
  dbDisconnect(con_main)
  dbDisconnect(con_add)
  dbDisconnect(con_out)
  
  cat("Database merge completed. Output saved to:", output_db, "\n")
}