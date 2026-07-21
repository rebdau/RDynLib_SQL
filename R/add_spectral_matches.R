#' @title adding names from a similarity table to the 'name' column in a sql 
#'         database.
#'         
#' @param con_merged 'character(1)' connection to the sqlite database.
#' 
#' @param sim_table 'character(1)' path to the similarity table with the matches 
#'  of some compounds.
#'  
#'  @return fill the 'name' column in the ms_compounds table the names from 
#'  similar compounds.
#'  
#'  @import DBI
#'  @import dplyr
#'  @import tools
#'  
#'  @author Ahlam Mentag
#'  
#'  @export
add_spectral_matches <- function(con_merged, sim_table) {
  
  
  # Build lookup table directly from similarity results
  sim_lookup <- sim_table %>%
    transmute(
      
      dataOrigin,
      acquisitionNum,
      msLevel,
      spectrum_type = query_spectrum_type,
      
      matches = paste0(
        "!",
        matched_peaks_count,
        "!",
        score,
        "!",
        target_name
      ),
      
      match_name = target_name,
      
      match_compound_id = target_compound_id,
      
      match_precursor_mz = target_precursorMz,
      
      match_score = score,
      
      match_common = matched_peaks_count,
      
      match_mslevel = target_msLevel,
      
      match_spectrum_type = target_spectrum_type
    ) %>%
    distinct()
  
  
  
  # Columns to create
  new_columns <- c(
    matches = "TEXT",
    match_name = "TEXT",
    match_compound_id = "INTEGER",
    match_precursor_mz = "REAL",
    match_score = "REAL",
    match_common = "INTEGER",
    match_mslevel = "INTEGER",
    match_spectrum_type = "TEXT"
  )
  
  
  existing_cols <- dbListFields(
    con_merged,
    "msms_spectrum"
  )
  
  
  for (col in names(new_columns)) {
    
    if (!(col %in% existing_cols)) {
      
      dbExecute(
        con_merged,
        sprintf(
          "ALTER TABLE msms_spectrum ADD COLUMN %s %s",
          col,
          new_columns[[col]]
        )
      )
    }
  }
  
  
  
  tryCatch({
    
    dbExecute(
      con_merged,
      "BEGIN TRANSACTION"
    )
    
    
    dbWriteTable(
      con_merged,
      "tmp_matches",
      sim_lookup,
      temporary = TRUE,
      overwrite = TRUE
    )
    
    
    dbExecute(
      con_merged,
      "
      UPDATE msms_spectrum

      SET

      matches = (
        SELECT matches
        FROM tmp_matches t
        WHERE t.dataOrigin = msms_spectrum.dataOrigin
        AND t.acquisitionNum = msms_spectrum.acquisitionNum
        AND t.msLevel = msms_spectrum.ms_level
        AND t.spectrum_type = msms_spectrum.spectrum_type
      ),

      match_name = (
        SELECT match_name
        FROM tmp_matches t
        WHERE t.dataOrigin = msms_spectrum.dataOrigin
        AND t.acquisitionNum = msms_spectrum.acquisitionNum
        AND t.msLevel = msms_spectrum.ms_level
        AND t.spectrum_type = msms_spectrum.spectrum_type
      ),

      match_compound_id = (
        SELECT match_compound_id
        FROM tmp_matches t
        WHERE t.dataOrigin = msms_spectrum.dataOrigin
        AND t.acquisitionNum = msms_spectrum.acquisitionNum
        AND t.msLevel = msms_spectrum.ms_level
        AND t.spectrum_type = msms_spectrum.spectrum_type
      ),

      match_precursor_mz = (
        SELECT match_precursor_mz
        FROM tmp_matches t
        WHERE t.dataOrigin = msms_spectrum.dataOrigin
        AND t.acquisitionNum = msms_spectrum.acquisitionNum
        AND t.msLevel = msms_spectrum.ms_level
        AND t.spectrum_type = msms_spectrum.spectrum_type
      ),

      match_score = (
        SELECT match_score
        FROM tmp_matches t
        WHERE t.dataOrigin = msms_spectrum.dataOrigin
        AND t.acquisitionNum = msms_spectrum.acquisitionNum
        AND t.msLevel = msms_spectrum.ms_level
        AND t.spectrum_type = msms_spectrum.spectrum_type
      ),

      match_common = (
        SELECT match_common
        FROM tmp_matches t
        WHERE t.dataOrigin = msms_spectrum.dataOrigin
        AND t.acquisitionNum = msms_spectrum.acquisitionNum
        AND t.msLevel = msms_spectrum.ms_level
        AND t.spectrum_type = msms_spectrum.spectrum_type
      ),

      match_mslevel = (
        SELECT match_mslevel
        FROM tmp_matches t
        WHERE t.dataOrigin = msms_spectrum.dataOrigin
        AND t.acquisitionNum = msms_spectrum.acquisitionNum
        AND t.msLevel = msms_spectrum.ms_level
        AND t.spectrum_type = msms_spectrum.spectrum_type
      ),

      match_spectrum_type = (
        SELECT match_spectrum_type
        FROM tmp_matches t
        WHERE t.dataOrigin = msms_spectrum.dataOrigin
        AND t.acquisitionNum = msms_spectrum.acquisitionNum
        AND t.msLevel = msms_spectrum.ms_level
        AND t.spectrum_type = msms_spectrum.spectrum_type
      )

      WHERE EXISTS (

        SELECT 1
        FROM tmp_matches t

        WHERE t.dataOrigin = msms_spectrum.dataOrigin
        AND t.acquisitionNum = msms_spectrum.acquisitionNum
        AND t.msLevel = msms_spectrum.ms_level
        AND t.spectrum_type = msms_spectrum.spectrum_type
      )
      "
    )
    
    
    dbExecute(
      con_merged,
      "COMMIT"
    )
    
    
  }, error = function(e) {
    
    dbExecute(
      con_merged,
      "ROLLBACK"
    )
    
    stop(e)
    
  })
  
  
  invisible(nrow(sim_lookup))
}