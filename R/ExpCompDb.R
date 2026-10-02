#' @title Add experiment information to a CompDb object
#'
#' @description
#' Helper function that creates a `CompDb` object from a database and adds a
#' join definition between the `ms_compound` and `experiment` tables using
#' the `expid` column.
#'
#' @param x file path to the sql database passed to `CompDb()` to create 
#'      the compound database object.
#'
#' @return A `CompDb` object with an additional join definition between the
#'   `ms_compound` and `experiment` tables.
#'
#' @author Ahlam Mentag
#'
#' @noRd
ExpCompDb <- function(x) {
  cdb <- suppressMessages(CompDb(x))
  CompoundDb::addJoinDefinition(cdb, "ms_compound", "experiment", "expid", "expid")
}