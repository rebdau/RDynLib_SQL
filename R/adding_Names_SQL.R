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
adding_Names_SQL <- function(con_merged, sim_table) {
  
  # Read tables
  ms_compound <- dbReadTable(con_merged, "ms_compound")
  msms_spectrum <- dbReadTable(con_merged, "msms_spectrum")
  
  # Keep only relevant columns from sim_table
  sim_sub <- sim_table %>%
    select(
      feature_id,
      dataOrigin,
      matched_peaks_count,
      score,
      target_name
    ) %>%
    distinct() %>%
    # Extract just the filename
    mutate(file_name = basename(dataOrigin))
  
  # extract just the filename in msms_spectrum
  msms_spectrum <- msms_spectrum %>%
    mutate(file_name = basename(dataOrigin))
  
  # Join sim_table with msms_spectrum using only the file name
  sim_msms <- sim_sub %>%
    inner_join(
      msms_spectrum %>% select(compound_id, file_name),
      by = "file_name"
    )
  
  # Join with ms_compound via nodename = feature_id
  sim_compound <- sim_msms %>%
    inner_join(
      ms_compound %>% select(compound_id, nodename),
      by = "compound_id"
    ) %>%
    filter(nodename == feature_id) %>%
    mutate(
      # Build the name in the desired format
      name = paste0(
        "!",
        matched_peaks_count,
        "!",
        score,
        "!",
        target_name
      )
    ) %>%
    select(compound_id, name) %>%
    distinct()
  
  # Update ms_compound 
  tryCatch({
    
    dbExecute(con_merged, "BEGIN TRANSACTION")
    
    dbWriteTable(
      con_merged,
      "tmp_update_name",
      sim_compound,
      temporary = TRUE,
      overwrite = TRUE
    )
    
    dbExecute(con_merged, "
      UPDATE ms_compound
      SET name = (
        SELECT tmp_update_name.name
        FROM tmp_update_name
        WHERE tmp_update_name.compound_id = ms_compound.compound_id
      )
      WHERE compound_id IN (
        SELECT compound_id FROM tmp_update_name
      )
    ")
    
    dbExecute(con_merged, "COMMIT")
    
  }, error = function(e) {
    dbExecute(con_merged, "ROLLBACK")
    stop(e)
  })
  
  invisible(nrow(sim_compound))
}
