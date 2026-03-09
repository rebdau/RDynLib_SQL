#' @title Extracting the aligned compound IDs with their database names of given 
#'        compounds.
#'
#' @param compound_id 'vector()' given compound_id, single value a vector of 
#'        compounds, to extract their aligned compounds.
#'        
#' @param database_name 'character(1)' database name of the given compound_id.
#'
#' @param alignment_file 'character(1)' the file path that contains the 
#' alignment results.
#'
#' @return a 'data.frame' with compid and database.
#'
#' @author Ahlam Mentag
#'
#' @export
getAlignedCompids_SQL <- function(compound_id, database_name, alignment_file) {
  
  # Read alignment file
  assoc <- read.table(alignment_file, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
  
  # compound_id must be a vector
  compound_id <- as.character(compound_id)
  
  res_list <- list()
  
  for (cid in compound_id) {
    
    # Case 1: compound is in reference database
    ref_matches <- assoc[
      assoc$ref_compid == cid & assoc$ref_database == database_name,
      c("target_compid", "target_database")
    ]
    
    if (nrow(ref_matches) > 0) {
      res_list[[length(res_list) + 1]] <- data.frame(
        query_compid = cid,
        compid = ref_matches$target_compid,
        database = ref_matches$target_database,
        stringsAsFactors = FALSE
      )
    }
    
    # Case 2: compound is in target database
    target_matches <- assoc[
      assoc$target_compid == cid & assoc$target_database == database_name,
      c("ref_compid", "ref_database")
    ]
    
    if (nrow(target_matches) > 0) {
      res_list[[length(res_list) + 1]] <- data.frame(
        query_compid = cid,
        compid = target_matches$ref_compid,
        database = target_matches$ref_database,
        stringsAsFactors = FALSE
      )
    }
  }
  
  if (length(res_list) == 0) {
    return(data.frame())
  }
  
  result <- do.call(rbind, res_list)
  rownames(result) <- NULL
  
  return(result)
}