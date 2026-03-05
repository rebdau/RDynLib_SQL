ExpCompDb <- function(x) {
  cdb <- suppressMessages(CompDb(x))
  CompoundDb::addJoinDefinition(cdb, "ms_compound", "experiment", "expid", "expid")
}