#' @title Interactively annotate and update compound names in an SQLite database
#'
#' @description
#'
#' `EnterName_SQL()` interactively constructs a structured compound name based on
#' user input and updates the corresponding entry in the \code{ms_compound}
#' table of an SQLite database.
#' 
#' @details
#'
#' The function reads the \code{ms_compound} table from a given SQLite database
#' and retrieves the compound associated with the provided compound identifier.
#' The current annotation metadata are displayed to the user, who is then guided
#' through a sequence of interactive prompts to construct a standardized compound
#' name.
#'
#' The name construction includes:
#' \itemize{
#'   \item annotator initials,
#'   \item structural confidence level (ID, AN, KA, PU),
#'   \item buffer adducts and isotopic labeling,
#'   \item in-source fragments and odd-electron ions,
#'   \item ion–neutral complexes (dimers or heteromers),
#'   \item indication whether the name refers to a molecular substructure.
#' }
#'
#' Optionally, the user may also update additional metadata such as the molecular
#' formula, SMILES string, and ppm mass deviation. After user confirmation, the
#' function updates the database directly using SQL queries.
#' Only fields explicitly provided by the user are modified.
#'
#' @param data_path `character(1)`
#'   Path to the SQLite database file containing the \code{ms_compound} table.
#'
#' @param CoMPID `numeric(1)`
#'   Unique compound identifier corresponding to
#'   \code{ms_compound.compound_id}.
#'
#' @return
#' Invisibly returns \code{TRUE} if the database entry was successfully updated.
#' If the user aborts the confirmation step, the function returns \code{NULL}
#' and no changes are written to the database.
#'
#' @importFrom DBI dbConnect
#'
#' @importFrom DBI dbDisconnect
#'
#' @importFrom DBI dbGetQuery
#'
#' @importFrom DBI dbExecute
#'
#' @importFrom RSQLite SQLite
#'
#' @importFrom utils readline
#'
#' @author Ahlam Mentag
#'
#' @export
EnterName_SQL <- function(data_path, CoMPID){
  
  # connect to database
  data_con <- dbConnect(SQLite(), data_path)
  on.exit(dbDisconnect(data_con), add = TRUE)
  
  # read compound
  cmp_data <- dbGetQuery(
    data_con,
    paste0("SELECT * FROM ms_compound WHERE compound_id = ", as.integer(CoMPID))
  )
  
  if (nrow(cmp_data) == 0) {
    stop("No compound found for compound_id = ", CoMPID)
  }
  
  CmpCSV <- cmp_data[, c(
    "retention_time", "mass_measured", "compound_id",
    "expid", "name", "smiles", "nodename", "formula",
    "ppm_deviation"
  )]
  
  dbkey <- 1L
  ## print current information
  cat("Currently, this info for your compound is present:\n")
  cat("COMPID: ", CmpCSV[dbkey, "compound_id"], "\n")
  cat("COMPNAME: ", CmpCSV[dbkey, "name"], "\n")
  cat("FORMULA: ", CmpCSV[dbkey, "formula"], "\n")
  cat("MASS_MEASURED: ", CmpCSV[dbkey, "mass_measured"], "\n")
  cat("PPM_DEVIATION: ", CmpCSV[dbkey, "ppm_deviation"], "\n")
  cat("RETENTION_TIME: ", CmpCSV[dbkey, "retention_time"], "\n")
  cat("EXPID: ", CmpCSV[dbkey, "expid"], "\n")
  cat("SMILES: ", CmpCSV[dbkey, "smiles"], "\n\n")
  
  # user initials
  rm_sel <- ""
  while (rm_sel == "")
    rm_sel <- toupper(readline("Initials first and last name: "))
  
  if (nchar(rm_sel) != 2) {
    stop("Initials must be exactly 2 characters.")
  }
  
  # confidence level
  cat(
    "Info regarding structural confidence of the given name:\n",
    "ID = identified via NMR or spiking of standard\n",
    "AN = annotated\n",
    "KA = characterized\n",
    "PU = putative\n"
  )
  
  rm_sel1 <- ""
  while (rm_sel1 == "")
    rm_sel1 <- toupper(readline("Structural confidence level: "))
  
  if (nchar(rm_sel1) != 2) {
    stop("Structural confidence level must be exactly 2 characters.")
  }
  rm_sel1 <- paste0(rm_sel1, " ")
  
  # compound name
  cat("For the following name, replace comma by underscore and primes by 'pr'\n")
  rm_sel2 <- ""
  while (rm_sel2 == "")
    rm_sel2 <- tolower(readline("Name: "))
  
  rm_sel2 <- paste0(rm_sel2, " ")
  
  # modifiers
  rm_sel3 <- substr(toupper(readline("Adduct (enter if none): ")), 1, 4)
  if (rm_sel3 != "") rm_sel3 <- paste0(rm_sel3, " ")
  
  rm_sel4 <- readline("Isotope (enter if none): ")
  if (rm_sel4 != "") rm_sel4 <- paste0(rm_sel4, " ")
  
  rm_sel5 <- readline("In-source fragment? (enter if no): ")
  if (rm_sel5 != "") rm_sel5 <- "ISF "
  
  rm_sel6 <- readline("Odd electron ion? (enter if no): ")
  if (rm_sel6 != "") rm_sel6 <- "OE "
  
  rm_sel8 <- readline("Ion-neutral complex of molecule? (enter if no): ")
  if (rm_sel8 != "") rm_sel8 <- "DIM "
  
  rm_sel9 <- readline("Ion-neutral complex? (enter if no): ")
  if (rm_sel9 != "") rm_sel9 <- "ADD "
  
  ## part-of-molecule flag
  rm_sel7 <- readline("Is the name referring to part of the molecule? (y/n): ")
  
  if (tolower(rm_sel7) == "n") {
    rmsl <- paste0(
      rm_sel2, rm_sel6, rm_sel5, rm_sel3,
      rm_sel9, rm_sel8, rm_sel4, rm_sel1, rm_sel
    )
  } else {
    rmsl <- paste0(
      "[", rm_sel2, rm_sel6, rm_sel5, rm_sel3,
      rm_sel9, rm_sel8, rm_sel4, rm_sel1, rm_sel, "]"
    )
  }
  
  # optional metadata
  formul <- gsub(" ", "", toupper(readline("Formula (enter if none): ")), fixed = TRUE)
  smiles <- readline("SMILES (enter if none): ")
  ppm    <- readline("PPM deviation (enter if none): ")
  
  # summary
  cat("\nNew data to be entered:\n")
  cat("COMPNAME: ", rmsl, "\n")
  cat("FORMULA: ", formul, "\n")
  cat("PPM_DEVIATION: ", ppm, "\n")
  cat("SMILES: ", smiles, "\n")
  
  SATIS <- readline("Are you satisfied? Type 'yes' to confirm: ")
  
  if (tolower(SATIS) != "yes") {
    message(
      "Oops! Looks like you made a mistake.\n",
      "Call EnterName_SQL(data_path, CoMPID) again."
    )
    return(invisible(NULL))
  }
  
  # build SQL update
  update_fields <- c("name = :name")
  params <- list(name = rmsl, compound_id = as.integer(CoMPID))
  
  if (formul != "") {
    update_fields <- c(update_fields, "formula = :formula")
    params$formula <- formul
  }
  if (ppm != "") {
    update_fields <- c(update_fields, "ppm_deviation = :ppm")
    params$ppm <- ppm
  }
  if (smiles != "") {
    update_fields <- c(update_fields, "smiles = :smiles")
    params$smiles <- smiles
  }
  
  sql <- paste0(
    "UPDATE ms_compound SET ",
    paste(update_fields, collapse = ", "),
    " WHERE compound_id = :compound_id"
  )
  
  dbExecute(data_con, sql, params = params)
  
  message("ms_compound successfully updated for compound_id = ", CoMPID)
  
  invisible(TRUE)
}
