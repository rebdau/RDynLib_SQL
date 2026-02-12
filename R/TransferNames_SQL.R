#' @title Transfer compound names between aligned compounds
#'
#' @description TransferNames_SQL() function
#' resolves and propagates compound names across aligned
#' compounds defined in an association table. It operates exclusively on
#' compound IDs and updates in-memory compound tables produced by
#' \code{Files_TransferNames_SQL()}.
#'
#' @param Names_lst A list as returned by \code{Files_TransferNames_SQL()},
#'   containing reference, target, and synthetic compound tables, as well as
#'   the association data frame.
#'
#' @return A list with two elements:
#' \describe{
#'   \item{Names_lst}{Updated compound tables with transferred names applied.}
#'   \item{Mult.frame}{A data.frame reporting compound IDs involved in multiple
#'     name transfer events, or the string \code{"empty"} if none are detected.}
#' }
#' 
#' @author Ahlam Mentag
#' @export
TransferNames_SQL <- function(Names_lst) {
  
  Assoc <- Names_lst[[4]]
  No.name <- c("NULL", "!")
  
  Mult.1 <- Mult.2 <- Mult.3 <- Mult.4 <- c()
  
  for (i in seq_len(nrow(Assoc))) {
    sel <- c(1, 2)
    
    ## get names
    name.cha <- Name_cha_SQL(sel, Names_lst, Assoc, i)
    HasName  <- Has_Name(name.cha, No.name, sel)
    
    if (sum(HasName) == 0)
      next
    
    Assoc_line <- c(
      Assoc$ref_compid[i],
      Assoc$target_compid[i]
    )
    
    full.lst <- Real_Name(
      HasName = HasName,
      sel = sel,
      Assoc_line = Assoc_line,
      Names_lst = Names_lst
    )
    
    Names_lst <- full.lst[[1]]
    
    if (length(full.lst[[2]]) != 1) {
      Mult.1 <- c(Mult.1, full.lst[[2]][1])
      Mult.2 <- c(Mult.2, full.lst[[2]][2])
      Mult.3 <- c(Mult.3, full.lst[[2]][3])
      Mult.4 <- c(Mult.4, full.lst[[2]][4])
    }
  }
  
  if (length(Mult.1) > 0) {
    Mult.frame <- data.frame(
      COMPID = as.integer(Mult.1),
      REF    = as.character(Mult.2),
      TARGET = as.character(Mult.3),
      SYN    = as.character(Mult.4),
      stringsAsFactors = FALSE
    )
  } else {
    Mult.frame <- "empty"
  }
  
  return(list(Names_lst, Mult.frame))
}
