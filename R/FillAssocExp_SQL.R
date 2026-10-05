#' @title Associate compounds between experiments in a DynLib database
#'
#' @description
#' Identifies and records associations between compounds from two experiments
#' stored in the same DynLib SQL database. For each compound in the reference
#' experiment, candidate compounds from the target experiment are searched
#' based on measured mass and retention-time criteria using
#' `Find_pos_compound()`. Retained associations are added to an association
#' table together with the corresponding database information.
#'
#' @param con A DBI connection to the DynLib SQLite database containing the
#'   reference and target experiments.
#' @param Assoc An optional data frame containing existing compound
#'   associations. If `NULL`, a new association table is initialized.
#' @param ref_expnr Integer. Experiment identifier (`expid`) of the reference
#'   experiment.
#' @param target_expnr Integer. Experiment identifier (`expid`) of the target
#'   experiment.
#' @param rg Numeric. Retention-time range parameter passed to
#'   `Find_pos_compound()`.
#' @param lc.err Numeric. Retention-time tolerance used for matching compounds
#'   between experiments.
#' @param err Numeric. Mass tolerance used for matching compounds between
#'   experiments.
#' @param db_path Character. Path to the SQLite database. The database file
#'   name is stored in the resulting association table.
#'
#' @return A data frame containing the compound associations with the columns
#'   `ref_compid`, `target_compid`, `ref_database`, and `target_database`.
#'   Duplicate associations are removed.
#'
#' @author Ahlam Mentag
#'
#' @export
FillAssocExp_SQL <- function(
    con,
    Assoc = NULL,
    ref_expnr,
    target_expnr,
    rg,
    lc.err,
    err,
    db_path
) {
  
  # Initialize association table
  if (is.null(Assoc)) {
    Assoc <- data.frame(
      ref_compid = integer(0),
      target_compid = integer(0),
      ref_database = character(0),
      target_database = character(0),
      stringsAsFactors = FALSE
    )
  }
  
  
  # Extract compounds from reference experiment
  ref.exp <- DBI::dbGetQuery(
    con,
    sprintf(
      "
      SELECT
        retention_time,
        mass_measured,
        compound_id,
        expid
      FROM ms_compound
      WHERE expid = %d
      ",
      ref_expnr
    )
  )
  
  
  # Extract compounds from target experiment
  target.exp <- DBI::dbGetQuery(
    con,
    sprintf(
      "
      SELECT
        retention_time,
        mass_measured,
        compound_id,
        expid
      FROM ms_compound
      WHERE expid = %d
      ",
      target_expnr
    )
  )
  
  
  if (nrow(ref.exp) == 0 || nrow(target.exp) == 0) {
    warning("One or both experiments contain no compounds.")
    return(Assoc)
  }
  
  
  # Order by measured mass
  ref.o <- ref.exp[
    order(ref.exp$mass_measured),
    ,
    drop = FALSE
  ]
  
  target.o <- target.exp[
    order(target.exp$mass_measured),
    ,
    drop = FALSE
  ]
  
  
  # Same database for reference and target
  db_name <- basename(db_path)
  
  
  # Search target compound for each reference compound
  for (i in seq_len(nrow(ref.o))) {
    
    mass.ref <- ref.o$mass_measured[i]
    time.ref <- ref.o$retention_time[i]
    compid.ref <- ref.o$compound_id[i]
    
    
    pres <- Find_pos_compound(
      mass.ref,
      time.ref,
      target.o,
      err,
      lc.err,
      rg
    )
    
    
    if (is.null(pres) || nrow(pres) == 0) {
      next
    }
    
    if (ncol(pres) < 3) {
      next
    }
    
    
    new_row <- data.frame(
      ref_compid = as.integer(compid.ref),
      target_compid = as.integer(pres[1, 3]),
      ref_database = db_name,
      target_database = db_name,
      stringsAsFactors = FALSE
    )
    
    
    Assoc <- rbind(
      Assoc,
      new_row
    )
  }
  
  
  # Remove duplicated associations
  Assoc <- unique(Assoc)
  
  
  return(Assoc)
}