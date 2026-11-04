#' @title Add Names to Compounds Based on Spectral Matches
#'
#' @description
#' Updates the `name` column of the `ms_compound` table using spectral
#' similarity matches stored in the `msms_spectrum` table. Compound names are
#' constructed from `match_common`, `match_score`, and `match_name`.
#'
#' By default, only matches associated with MS2 spectra (`msLevel = 2`) are
#' considered. 
#'
#' The function can optionally be restricted to a selected set of compounds
#' by providing their compound identifiers through the `compound_id` argument.
#' If `compound_id = NULL` (default), all compounds are considered.
#'
#' By default, only compounds without an existing name are updated.
#' Existing names can be replaced by setting `overwrite = TRUE`.
#'
#' Each name addition or modification is recorded as a new row in the
#' `annotation_history` table if it exists, otherwise the table will be created.
#' If the `annotation_events` column does not
#' already exist in `annotation_history`, it is created automatically.
#' The annotation event stores information about the spectra involved in the
#' match, including the matched compound accession, spectrum types, and MS
#' levels.
#'
#' @param con_merged A `DBIConnection` to the SQLite DynLib database.
#' 
#' @param spectrum_type `character` the spectrum type of the compounds where 
#'   the names will be added.
#'
#' @param msLevel `integer(1)`. MS level of the query spectra to consider when
#'   selecting spectral matches. The default is `2`, corresponding to MS2
#'   spectra. Set to `NULL` to consider matches from all MS levels.
#'
#' @param compound_id `integer` or `NULL`. Optional vector of compound
#'   identifiers specifying which compounds should be considered for name
#'   assignment. If `NULL` (default), all compounds are considered.
#'
#' @param overwrite `logical(1)`. If `FALSE` (default), names are added only
#'   when the `name` field in `ms_compound` is `NULL` or empty. If `TRUE`,
#'   existing compound names can also be replaced by the name derived from
#'   the spectral match.
#'
#' @return Invisibly returns the number of compound names that were added or
#'   updated. The function modifies the `ms_compound` and `annotation_history`
#'   tables in the database.
#'
#' @import DBI
#' @import dplyr
#'
#' @author Ahlam Mentag
#' @export

adding_Names_SQL <- function(
    con_merged,
    spectrum_type,
    msLevel = 2,
    compound_id = NULL,
    overwrite = FALSE
) {
  
  # Check required tables
  required_tables <- c(
    "ms_compound",
    "msms_spectrum"
  )
  
  missing_tables <- setdiff(
    required_tables,
    DBI::dbListTables(con_merged)
  )
  
  if (length(missing_tables) > 0) {
    stop(
      "Missing table(s): ",
      paste(missing_tables, collapse = ", ")
    )
  }
  
  # Check spectrum_type argument
  if (missing(spectrum_type) ||
      !is.character(spectrum_type) ||
      length(spectrum_type) != 1L ||
      is.na(spectrum_type) ||
      !nzchar(spectrum_type)) {
    
    stop(
      "`spectrum_type` must be provided as a single non-empty character value."
    )
  }
  
  # Check msLevel argument
  if (!is.null(msLevel)) {
    
    if (length(msLevel) != 1L || is.na(msLevel)) {
      stop("`msLevel` must be a single MS level or NULL.")
    }
  }
  
  # Create annotation_history table if it does not exist
  DBI::dbExecute(
    con_merged,
    '
    CREATE TABLE IF NOT EXISTS annotation_history (
      annotation_id INTEGER PRIMARY KEY AUTOINCREMENT,
      expid INTEGER,
      compound_id INTEGER NOT NULL,
      date TEXT NOT NULL,
      type TEXT,
      description TEXT,
      "references" TEXT,
      annotation_events TEXT,

      FOREIGN KEY(compound_id)
        REFERENCES ms_compound(compound_id),

      UNIQUE(compound_id, date)
    )
    '
  )
  
  # Check existing annotation_history columns
  annotation_fields <- DBI::dbListFields(
    con_merged,
    "annotation_history"
  )
  
  # Rename reference to references if needed
  if ("reference" %in% annotation_fields &&
      !"references" %in% annotation_fields) {
    
    DBI::dbExecute(
      con_merged,
      '
      ALTER TABLE annotation_history
      RENAME COLUMN "reference" TO "references"
      '
    )
  }
  
  # Refresh column names
  annotation_fields <- DBI::dbListFields(
    con_merged,
    "annotation_history"
  )
  
  # Add references column if missing
  if (!"references" %in% annotation_fields) {
    
    DBI::dbExecute(
      con_merged,
      '
      ALTER TABLE annotation_history
      ADD COLUMN "references" TEXT
      '
    )
  }
  
  # Add annotation_events column if missing
  if (!"annotation_events" %in% annotation_fields) {
    
    DBI::dbExecute(
      con_merged,
      "
      ALTER TABLE annotation_history
      ADD COLUMN annotation_events TEXT
      "
    )
  }
  
  # Check for duplicate annotations
  duplicates <- DBI::dbGetQuery(
    con_merged,
    "
    SELECT compound_id, date, COUNT(*) AS n
    FROM annotation_history
    GROUP BY compound_id, date
    HAVING COUNT(*) > 1
    "
  )
  
  if (nrow(duplicates) > 0) {
    stop(
      "Cannot create UNIQUE constraint: ",
      "duplicate (compound_id, date) pairs exist ",
      "in annotation_history."
    )
  }
  
  # Ensure UNIQUE constraint for existing databases
  DBI::dbExecute(
    con_merged,
    "
    CREATE UNIQUE INDEX IF NOT EXISTS
      idx_annotation_compound_date
    ON annotation_history(compound_id, date)
    "
  )
  
  # Get available matches
  candidates <- DBI::dbGetQuery(
    con_merged,
    "
    SELECT
      c.compound_id,
      c.compound_accession,
      c.expid,
      c.name AS current_name,

      s.spectrum_id,
      s.ms_level,
      s.spectrum_type,

      s.match_name,
      s.match_compound_id,
      s.match_score,
      s.match_common,
      s.match_mslevel,
      s.match_spectrum_type,
      s.match_spectrum_id

    FROM ms_compound AS c

    INNER JOIN msms_spectrum AS s
      ON c.compound_id = s.compound_id

    WHERE s.match_name IS NOT NULL
      AND s.match_common IS NOT NULL
      AND s.match_score IS NOT NULL

    ORDER BY
      c.compound_id,
      s.match_score DESC
    "
  )
  
  if (nrow(candidates) == 0) {
    message("No matches available.")
    return(invisible(0))
  }
  
  # Check that requested spectrum type exists
  available_types <- unique(candidates$spectrum_type)
  available_types <- available_types[!is.na(available_types)]
  
  if (!spectrum_type %in% available_types) {
    
    stop(
      "`spectrum_type` '",
      spectrum_type,
      "' was not found. Available spectrum types: ",
      paste(available_types, collapse = ", ")
    )
  }
  
  # Filter by spectrum type
  candidates <- candidates[
    !is.na(candidates$spectrum_type) &
      candidates$spectrum_type == spectrum_type,
    ,
    drop = FALSE
  ]
  
  # Filter by MS level
  if (!is.null(msLevel)) {
    
    candidates <- candidates[
      !is.na(candidates$ms_level) &
        candidates$ms_level == msLevel,
      ,
      drop = FALSE
    ]
  }
  
  # Filter by selected compound IDs
  if (!is.null(compound_id)) {
    
    candidates <- candidates[
      candidates$compound_id %in% compound_id,
      ,
      drop = FALSE
    ]
  }
  
  if (nrow(candidates) == 0) {
    message("No matches available for the selected filters.")
    return(invisible(0))
  }
  
  # Keep only unnamed compounds unless overwrite = TRUE
  if (!overwrite) {
    
    candidates <- candidates[
      is.na(candidates$current_name) |
        candidates$current_name == "",
      ,
      drop = FALSE
    ]
  }
  
  if (nrow(candidates) == 0) {
    message("No compound names need to be updated.")
    return(invisible(0))
  }
  
  # Construct new name
  candidates$new_name <- paste0(
    "!",
    candidates$match_common,
    "!",
    candidates$match_score,
    "!",
    candidates$match_name
  )
  
  # Define annotation type
  candidates$annotation_type <- ifelse(
    is.na(candidates$current_name) |
      candidates$current_name == "",
    "name addition",
    "name update"
  )
  
  # Create annotation description
  candidates$annotation_description <- candidates$new_name
  
  # Get accession of matched compounds
  accessions <- DBI::dbGetQuery(
    con_merged,
    "
    SELECT
      compound_id AS match_compound_id,
      compound_accession AS match_compound_accession
    FROM ms_compound
    "
  )
  
  candidates <- merge(
    candidates,
    accessions,
    by = "match_compound_id",
    all.x = TRUE,
    sort = FALSE
  )
  
  # Restore original match ordering after merge
  candidates <- candidates[
    order(
      candidates$compound_id,
      -candidates$match_score
    ),
    ,
    drop = FALSE
  ]
  
  # Create annotation_events
  candidates$annotation_events <- paste0(
    "ms_level=",
    candidates$ms_level,
    "; spectrum_type=",
    candidates$spectrum_type,
    "; match_mslevel=",
    candidates$match_mslevel,
    "; match_spectrum_type=",
    candidates$match_spectrum_type,
    "; match_compound_accession=",
    candidates$match_compound_accession,
    "; match_spectrum_id=",
    candidates$match_spectrum_id
  )
  
  # Update database + annotation history
  DBI::dbWithTransaction(con_merged, {
    
    for (i in seq_len(nrow(candidates))) {
      
      # Update ms_compound$name
      DBI::dbExecute(
        con_merged,
        "
        UPDATE ms_compound
        SET name = ?
        WHERE compound_id = ?
        ",
        params = list(
          candidates$new_name[i],
          candidates$compound_id[i]
        )
      )
      
      # Insert or update annotation_history
      DBI::dbExecute(
        con_merged,
        '
        INSERT INTO annotation_history (
          expid,
          compound_id,
          date,
          type,
          description,
          annotation_events
        )
        VALUES (?, ?, ?, ?, ?, ?)

        ON CONFLICT(compound_id, date)
        DO UPDATE SET
          expid = excluded.expid,
          type = excluded.type,
          description = excluded.description,
          annotation_events = excluded.annotation_events
        ',
        params = list(
          candidates$expid[i],
          candidates$compound_id[i],
          format(Sys.Date(), "%d/%m/%Y"),
          candidates$annotation_type[i],
          candidates$annotation_description[i],
          candidates$annotation_events[i]
        )
      )
    }
  })
  
  message(
    nrow(candidates),
    " compound name update(s) processed and recorded in annotation_history."
  )
  
  invisible(nrow(candidates))
}
