#' This script converts the *Assoc* text files from the original to the new
#' *long* format which is more flexible and expandable.
#'
#' @author Ahlam Mentag
convert_compid_assoc <- function(input_path, output_path) {

  df <- read.table(input_path, header = TRUE, sep = "\t",
                   stringsAsFactors = FALSE)

  ref_col <- names(df)[1]
  target_cols <- names(df)[-1]

  # Keep rows where reference ID is not NA
  df <- df[!is.na(df[[ref_col]]), ]

  # change the databases name
  map_database_name <- function(col_name) {
    if (col_name == "COMPIDposSynapt") return("DyLib_QTOF_pos.sqlite")
    if (col_name == "COMPIDnegSynapt") return("DYL_FLAX_QTOFneg.sqlite")
    if (col_name == "COMPIDnegFT") return("DYL_FLAX_FTMSneg.sqlite")
    if (col_name == "COMPIDposFT") return("DyLib_FTMS_pos.sqlite")
    return(col_name)
  }

  long_rows <- lapply(seq_len(nrow(df)), function(i) {

    ref_compid <- df[[ref_col]][i]
    ref_database <- map_database_name(ref_col)

    # Extract target values as a named vector
    target_vec <- df[i, target_cols]
    target_vec <- unlist(target_vec, use.names = TRUE)
    target_vec <- target_vec[!is.na(target_vec)]

    # Skip rows with no targets
    if (length(target_vec) == 0) return(NULL)

    data.frame(
      ref_compid = ref_compid,
      target_compid = as.numeric(target_vec),
      ref_database = ref_database,
      target_database = sapply(names(target_vec), map_database_name),
      stringsAsFactors = FALSE
    )
  })

  long_df <- do.call(rbind, long_rows)

  write.table(long_df, output_path, sep = "\t", quote = FALSE, row.names = FALSE)

  long_df
}

## Run the conversion.
setwd("C:/Users/amentag/Desktop/these/These/clones/new_clones/RDynLib/New_Dynlib")


# Process each file and store the results
ftneg_df  <- convert_compid_assoc("data/DynLib_subDB_alignment/DynLibDBassociation.txt",
                                  "data/New_assoc/ftneg.txt")

ftpos_df  <- convert_compid_assoc("data/DynLib_subDB_alignment/DynLibDBassociation2.txt",
                                  "data/New_assoc/ftpos.txt")

qtofneg_df <- convert_compid_assoc("data/DynLib_subDB_alignment/DynLibDBassociation3.txt",
                                   "data/New_assoc/qtofneg.txt")

# Combine all three into one data.frame
combined_df <- rbind(ftneg_df, ftpos_df, qtofneg_df)

write.table(combined_df, "data/combined_assoc.txt",
            sep = "\t", quote = FALSE, row.names = FALSE)

# Check the result
head(combined_df)
