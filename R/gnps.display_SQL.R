#' @title Display GNPS similarity results from an SQLite database.
#' 
#' @description
#'
#' Extracts and formats GNPS-based similarity information from the
#' `gnps_add` table in an SQLite database for a given `compound_id`.
#' Encoded similarity values stored in the table columns are parsed
#' and returned as a structured data frame.
#' 
#' @param sql_path `Character(1)` string giving the path to the SQLite database.
#' 
#' @param dbkey 'numeric(1)' the compound_id used to select the corresponding 
#' row in the `compound_add` table.
#' 
#' @param nr_col 'numeric(1)' maximum number of columns to scan in the table.
#' Default is 6.
#' 
#' @return A `data.frame` containing GNPS similarity results with the columns:
#' \describe{
#'   \item{massdiff}{Mass difference associated with the match.}
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
gnps.display_SQL <- function(sql_path, dbkey, nr_col2 = 6) {
  
  con <- dbConnect(RSQLite::SQLite(), sql_path)
  on.exit(dbDisconnect(con), add = TRUE)
  
  gnps_add <- dbGetQuery(con, "SELECT * FROM gnps_add")
  
  # Filter by compound_id
  row_idx <- which(gnps_add$compound_id == dbkey)
  if (length(row_idx) == 0) {
    warning("No matching compound_id found")
    return(data.frame())
  }
  
  row_data <- gnps_add[row_idx[1], ]  # take the first match if multiple
  
  gnps.res <- data.frame(
    massdiff = numeric(),
    compid.prod = integer(),
    ions.prod = integer(),
    ave.common = numeric(),
    ave.dot = numeric(),
    stringsAsFactors = FALSE
  )
  
  i <- 2
  j <- 1
  
  repeat {
    cell_value <- row_data[[i]]
    
    if (is.na(cell_value) || cell_value == "") {
      if (i == nr_col2) break
      i <- i + 1
      next
    }
    
    parts <- strsplit(cell_value, "!!")[[1]]
    if (length(parts) < 3) {
      if (i == nr_col2) break
      i <- i + 1
      next
    }
    
    subparts <- strsplit(parts[2], "!")[[1]]
    if (length(subparts) < 4) {
      if (i == nr_col2) break
      i <- i + 1
      next
    }
    
    compid.prod <- as.integer(subparts[4])
    if (is.na(compid.prod)) {
      if (i == nr_col2) break
      i <- i + 1
      next
    }
    
    gnps.res[j, "massdiff"] <- as.numeric(parts[3])
    gnps.res[j, "compid.prod"] <- compid.prod
    gnps.res[j, "ions.prod"] <- as.integer(subparts[1])
    gnps.res[j, "ave.common"] <- as.numeric(subparts[2])
    gnps.res[j, "ave.dot"] <- as.numeric(parts[3])
    
    if (i == nr_col2) break
    
    i <- i + 1
    j <- j + 1
  }
  
  return(gnps.res)
}
