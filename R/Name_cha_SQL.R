#' @title Retrieve compound names for an association entry
#'
#' @description Name_cha_SQL() function
#' extracts compound names corresponding to a single association
#' entry by matching compound identifiers against in-memory compound tables.
#' It supports explicit reference and target compound selection.
#'
#' @param sel Integer vector indicating which compound roles to retrieve
#'   (for reference and target).
#' @param Names_lst A list of compound tables as returned by
#'   \code{Files_TransferNames_SQL()}.
#' @param Assoc Association data frame containing compound ID mappings.
#' @param i Integer index specifying the row of \code{Assoc} to process.
#'
#' @return A character vector of compound names corresponding to the selected
#'   association entry. Missing names are returned as \code{NA}.
#'
#' @author Ahlam Mentag
#' @export
Name_cha_SQL <- function(sel, Names_lst, Assoc, i) {
  
  name.cha <- character(0)
  
  for (j in sel) {
    
    if (j == 1) {
      ## reference compound
      compid <- Assoc$ref_compid[i]
      df <- Names_lst[[1]]
    } else if (j == 2) {
      ## target compound
      compid <- Assoc$target_compid[i]
      df <- Names_lst[[2]]
    } else {
      stop("Unsupported selection index")
    }
    
    name <- df$COMPNAME[df$COMPID == as.character(compid)]
    
    if (length(name) == 0) {
      name <- NA_character_
    }
    
    name.cha <- c(name.cha, name[1])
  }
  
  return(name.cha)
}
