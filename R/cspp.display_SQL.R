#' @title Display compound similarity results from an SQLite database.
#' 
#' @description
#' Extracts and formats compound similarity information from the
#' `compound_add` table in an SQLite database for a given `compound_id`.
#' The function parses encoded similarity values stored in the table
#' columns and returns them in a structured data frame.
#' 
#' @param sql_path `Character(1)` string giving the path to the SQLite database.
#' 
#' @param dbkey 'numeric(1)' the compound_id used to select the corresponding 
#' row in the `compound_add` table.
#' 
#' @param nr_col 'numeric(1)' maximum number of columns to scan in the table.
#'   Default is 35.
#' 
#' @return A `data.frame` containing similarity information with the columns:
#' \describe{
#'   \item{conv.type}{Conversion type (column name in the database).}
#'   \item{compid.prod}{Product compound identifier.}
#'   \item{ions.prod}{Number of product ions.}
#'   \item{ave.common}{Average number of common ions.}
#'   \item{ave.dot}{Average dot-product similarity score.}
#' }
#'
#' Returns an empty `data.frame` if no matching `compound_id` is found.
#'
#' @author Ahlam Mentag
#' @export
cspp.display_SQL <- function(sql_path, dbkey, nr_col = 35) {
  
  con <- dbConnect(RSQLite::SQLite(), sql_path)
  on.exit(dbDisconnect(con), add = TRUE)
  
  compound_add <- dbGetQuery(con, "SELECT * FROM compound_add")
  
  # Filter by compound_id
  row_idx <- which(compound_add$compound_id == dbkey)
  if (length(row_idx) == 0) {
    warning("No matching compound_id found")
    return(data.frame())
  }
  # take the first match if multiple
  row_data <- compound_add[row_idx[1], ]  
  
  cspp.res <- data.frame(
    conv.type = character(),
    compid.prod = integer(),
    ions.prod = integer(),
    ave.common = numeric(),
    ave.dot = numeric(),
    stringsAsFactors = FALSE
  )
  
  i <- 3
  j <- 1
  
  repeat {
    cell_value <- row_data[[i]]
    
    if (is.na(cell_value) || cell_value == "") {
      if (i == nr_col) break
      i <- i + 1
      next
    }
    
    parts <- strsplit(cell_value, "!!")[[1]]
    if (length(parts) < 3) {
      if (i == nr_col) break
      i <- i + 1
      next
    }
    
    compid.prod <- as.integer(parts[3])
    if (is.na(compid.prod)) {
      if (i == nr_col) break
      i <- i + 1
      next
    }
    
    subparts <- strsplit(parts[2], "!")[[1]]
    if (length(subparts) < 3) {
      if (i == nr_col) break
      i <- i + 1
      next
    }
    
    cspp.res[j, "conv.type"] <- colnames(compound_add)[i]
    cspp.res[j, "compid.prod"] <- compid.prod
    cspp.res[j, "ions.prod"] <- as.integer(subparts[1])
    cspp.res[j, "ave.common"] <- as.numeric(subparts[2])
    cspp.res[j, "ave.dot"] <- as.numeric(subparts[3])
    
    if (i == nr_col) break
    
    i <- i + 1
    j <- j + 1
  }
  
  return(cspp.res)
}
