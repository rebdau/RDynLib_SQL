# take care that compound_add.txt is present in the working directory
# ave.cnt is the average of common neutrals and common ions. If half of the
# product ions are common and the other half yield neutral losses in common,
# than the ave.cnt is 0.5! Noticeably, 0.5 would also be obtained if that half
# of the product ions that are in common are also responsible for the common
# neutral losses.
# ave.dot refers to the average dot product obtained via either the common
# product ions or the common neutral losses.
# Comment: a better dot product might be based on all product ions of a parti-
# cular spectrum whether present or not in the other spectrum???? This would
# take the count into account in the dot product. Two possibilities are then a
# forward and reverse search.
# CSPPs are ordered first on decreasing number of ions in common, then on ave.cnt
# and finally on ave.dot. For each substrate COMPID, only one product COMPID is
# retained, i.e. the one that has the "highest similarity", namely that is the
# highest on the ordered CSPP.df list.
# Other CSPPs for the particular conversion type, e.g. methylation,
# that contain either the same substrate COMPID or the same product COMPID
# are removed from the cspp.df list before continuation. If you want to keep a
# complete record of all generated CSPPs, write the initial cspp.df file away to
# a .txt file.
# the conversion type allows to determine the column (indicated by conv.col) of
# comp_app (and finally of compound_add.txt) to which the results need to be written. 

rank.cspp <- function(cspp.df, conv.col, comp_add) {
  
  if (is.null(cspp.df) || nrow(cspp.df) == 0) return(comp_add)
  
  cspp.df <- as.data.frame(cspp.df)
  
  # safety check
  required <- c("COMPID.sub", "COMPID.prod",
                "FORW_IONS", "REV_IONS",
                "DOT_IONS", "DOT_LOSS",
                "COMMON_IONS")
  
  if (!all(required %in% names(cspp.df))) {
    stop("cspp.df missing required columns")
  }
  
  cspp.df$ave.cnt <- (cspp.df$FORW_IONS + cspp.df$REV_IONS) / 2
  cspp.df$ave.dot <- (cspp.df$DOT_IONS + cspp.df$DOT_LOSS) / 2
  
  cspp.df <- cspp.df[
    order(-cspp.df$COMMON_IONS,
          -cspp.df$ave.cnt,
          -cspp.df$ave.dot),
  ]
  
  rownames(comp_add) <- comp_add$compound_id
  
  while (nrow(cspp.df) > 0) {
    
    sub_id  <- cspp.df$COMPID.sub[1]
    prod_id <- cspp.df$COMPID.prod[1]
    
    cspp.name <- paste0(
      "!!", cspp.df$COMMON_IONS[1],
      "!", cspp.df$ave.cnt[1],
      "!", cspp.df$ave.dot[1],
      "!!", prod_id
    )
    
    comp_add[as.character(sub_id), conv.col] <- cspp.name
    
    rm_idx <- which(
      cspp.df$COMPID.sub == sub_id |
        cspp.df$COMPID.prod == prod_id
    )
    
    cspp.df <- cspp.df[-rm_idx, , drop = FALSE]
  }
  
  comp_add
}