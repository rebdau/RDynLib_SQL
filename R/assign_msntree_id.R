#' Assign fragmentation tree IDs to MSn spectra
#'
#' This function assigns a fragmentation tree ID (`MSntreeID`) to each spectrum
#' in a \linkS4class{Spectra} object based on precursor relationships within
#' each sample (`dataOrigin`). MS2 spectra always start a new tree, and downstream
#' MS3/MS4 spectra inherit the tree ID from their parent.
#'
#' Additionally, the MS2-level precursor m/z is propagated to all downstream
#' spectra (MS3, MS4, …) via the `precursorMz.MS2` variable, with `NA` for MS2 spectra.
#'
#' @param x A \linkS4class{Spectra} object containing MS2–MSn spectra.
#'
#' @return A \linkS4class{Spectra} object with two additional columns:
#'   \itemize{
#'     \item \code{MSntreeID} — integer ID of the fragmentation tree each spectrum belongs to
#'     \item \code{precursorMz.MS2} — MS2-level precursor m/z propagated to MS3/MS4 (and beyond), \code{NA} for MS2
#'   }
#'
#' @details
#' The function iterates over each sample (dataOrigin) and assigns tree IDs
#' hierarchically: MS2 spectra start new trees, and higher-level spectra
#' (MS3, MS4, etc.) inherit the tree ID from their parent MS level using the
#' \code{precScanNum} field. Precursor m/z at MS2 level is propagated to all
#' descendant spectra for easier feature linking.
#'
#' @examples
#' ## Assume 'msn' is a Spectra object containing MS2–MSn spectra
#' msn <- assign_msntree_id(msn)
#' msn$MSntreeID
#' msn$precursorMz.MS2
#'
#' @export
#' @importFrom Spectra msLevel scanIndex precScanNum dataOrigin precursorMz
assign_msntree_id <- function(x) {
  ## Order x by MS level: ensures we first assign an ID to the MS2.
  o <- order(msLevel(x))
  x <- x[o]
  scan_index <- scanIndex(x)
  prec_scan_num <- precScanNum(x)
  ms_level <- msLevel(x)
  data_origin <- dataOrigin(x)
  
  ## store precursorMz values
  precursor_mz <- precursorMz(x)
  
  tree_id <- rep(NA_integer_, length(scan_index))
  ## initialize output vector for MS2-level precursor m/z
  precursor_mz_MS2 <- rep(NA_real_, length(scan_index))
  
  global_counter <- 0
  
  for (origin in unique(data_origin)) {
    idx <- which(data_origin == origin)
    origin_scan_index <- scan_index[idx]
    origin_prec_scan_num <- prec_scan_num[idx]
    origin_ms_level <- ms_level[idx]
    ## per-origin precursorMz vector
    origin_precursor_mz <- precursor_mz[idx]
    origin_tree_id <- rep(NA_integer_, length(idx))
    ## per-origin precursorMz.MS2 vector
    origin_precursor_mz_MS2 <- rep(NA_real_, length(idx))
    
    scan_idx_map <- setNames(seq_along(origin_scan_index), origin_scan_index)
    
    for (i in seq_along(idx)) {
      level <- origin_ms_level[i]
      if (level == 2) {
        ## If MS2 hasn't been assigned yet
        if (is.na(origin_tree_id[i])) {
          global_counter <- global_counter + 1
          origin_tree_id[i] <- global_counter
        }
        ## MS2 spectra have no MS2-level precursor
        origin_precursor_mz_MS2[i] <- NA_real_
        
      } else {
        ## Assign ID of the related precursor ID
        parent_scan <- origin_prec_scan_num[i]
        parent_pos <- which(origin_scan_index == parent_scan)
        if (length(parent_pos) == 1 && origin_ms_level[parent_pos] == (level - 1)) {
          origin_tree_id[i] <- origin_tree_id[parent_pos]
          
          ## assign precursorMz.MS2 depending on depth
          if (level == 3) {
            origin_precursor_mz_MS2[i] <- origin_precursor_mz[i]
          } else {
            ancestor_pos <- parent_pos
            while (!is.na(ancestor_pos) && origin_ms_level[ancestor_pos] > 3) {
              ancestor_scan <- origin_prec_scan_num[ancestor_pos]
              ancestor_pos <- scan_idx_map[as.character(ancestor_scan)]
            }
            if (!is.na(ancestor_pos) && origin_ms_level[ancestor_pos] == 3) {
              origin_precursor_mz_MS2[i] <- origin_precursor_mz[ancestor_pos]
            } else {
              origin_precursor_mz_MS2[i] <- NA_real_
            }
          }
          
        } else {
          stop("MS", level, " with ", length(parent_pos), " precursors")
        }
      }
    }
    tree_id[idx] <- origin_tree_id
    ## store per-origin results
    precursor_mz_MS2[idx] <- origin_precursor_mz_MS2
  }
  
  ## Return in original order
  tree_id <- tree_id[order(o)]
  ## reorder new variable to match
  precursor_mz_MS2 <- precursor_mz_MS2[order(o)]
  
  x <- x[order(o)]
  
  ## add both to Spectra object
  x$MSntreeID <- tree_id
  x$precursorMz.MS2 <- precursor_mz_MS2
  
  return(x)
}
