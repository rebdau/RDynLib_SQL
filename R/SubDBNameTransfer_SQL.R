#' @title Transfer compound metadata across SQLite sub-databases
#'
#' @description
#' Transfers `name`, `subsid`, and `smiles` between associated compounds
#' from two SQLite databases. Only missing values are filled; existing values
#' are not overwritten and conflicting values are reported.
#'
#' @param sqlite_dir Character string. Directory containing the SQLite databases.
#' @param assoc_path Character string. Path to the tab-separated association file.
#'
#' @return Invisibly returns a list containing:
#'   `Mult.frame`, compounds with multiple mappings, and `conflicts`,
#'   existing values that differ from the proposed transferred values.
#'
#' @author Ahlam Mentag
#' @export
SubDBNameTransfer_SQL <- function(sqlite_dir, assoc_path) {
  
  # Read association file
  
  Assoc <- read.table(
    assoc_path,
    header = TRUE,
    sep = "\t",
    stringsAsFactors = FALSE
  )
  
  ref_db <- unique(Assoc$ref_database)
  target_db <- unique(Assoc$target_database)
  
  if (length(ref_db) != 1 || length(target_db) != 1) {
    stop(
      "Association file must contain exactly one ref_database and one target_database"
    )
  }
  
  ref_db <- ref_db[1]
  target_db <- target_db[1]
  
  cat("Reference database :", ref_db, "\n")
  cat("Target database    :", target_db, "\n")
  
  
  # Transfer metadata
  
  Names_lst <- Files_TransferNames_SQL(
    sqlite_dir = sqlite_dir,
    assoc_path = assoc_path
  )
  
  transfertNam <- TransferNames_SQL(Names_lst)
  
  Updated_Names <- transfertNam[[1]]
  Mult.frame <- transfertNam[[2]]
  
  
  
  is_missing <- function(x) {
    is.na(x) | as.character(x) == ""
  }
  
  
  # Update one database
  
  update_database <- function(db_name, proposed) {
    
    sqlite_path <- file.path(
      sqlite_dir,
      db_name
    )
    
    cat(
      "\nProcessing :",
      sqlite_path,
      "\n"
    )
    
    
    # Check database
    
    if (!file.exists(sqlite_path)) {
      stop(
        "Database file not found: ",
        sqlite_path
      )
    }
    
    
    # Open database
    
    con <- DBI::dbConnect(
      RSQLite::SQLite(),
      sqlite_path
    )
    
    on.exit(
      DBI::dbDisconnect(con),
      add = TRUE
    )
    
  
    
    if (!"ms_compound" %in% DBI::dbListTables(con)) {
      stop(
        "Table 'ms_compound' not found in ",
        db_name
      )
    }
    
    
    # Check transferred columns
    
    required_cols <- c(
      "compound_id",
      "name",
      "subsid",
      "smiles"
    )
    
    missing_cols <- setdiff(
      required_cols,
      names(proposed)
    )
    
    if (length(missing_cols) > 0) {
      stop(
        "Missing columns in transferred data: ",
        paste(
          missing_cols,
          collapse = ", "
        )
      )
    }
    
    
    # Read ms_compound ONCE
    
    current <- DBI::dbGetQuery(
      con,
      "
      SELECT
        compound_id,
        name,
        subsid,
        smiles
      FROM ms_compound
      "
    )
    
    
    # Keep only required columns
    
    proposed <- proposed[
      ,
      c(
        "compound_id",
        "name",
        "subsid",
        "smiles"
      ),
      drop = FALSE
    ]
    
    
    # Rename proposed columns
    names(proposed) <- c(
      "compound_id",
      "proposed_name",
      "proposed_subsid",
      "proposed_smiles"
    )
    
    
    # Match proposed compounds to database
    
    idx <- match(
      proposed$compound_id,
      current$compound_id
    )
    
    not_found <- is.na(idx)
    
    
    if (any(not_found)) {
      warning(
        sum(not_found),
        " compound_id(s) were not found in ",
        db_name
      )
    }
    
    
    # Keep only existing compound IDs
    proposed <- proposed[
      !not_found,
      ,
      drop = FALSE
    ]
    
    idx <- idx[!not_found]
    
    
    # Add current database values

    
    proposed$current_name <-
      current$name[idx]
    
    proposed$current_subsid <-
      current$subsid[idx]
    
    proposed$current_smiles <-
      current$smiles[idx]
    
    

    # Determine which missing values can be filled

    
    fill_name <-
      is_missing(proposed$current_name) &
      !is_missing(proposed$proposed_name)
    
    
    fill_subsid <-
      is_missing(proposed$current_subsid) &
      !is_missing(proposed$proposed_subsid)
    
    
    fill_smiles <-
      is_missing(proposed$current_smiles) &
      !is_missing(proposed$proposed_smiles)
    
    
    needs_update <-
      fill_name |
      fill_subsid |
      fill_smiles
    

    #
    # Existing values are NEVER overwritten
    
    conflict_name <-
      !is_missing(proposed$current_name) &
      !is_missing(proposed$proposed_name) &
      proposed$current_name != proposed$proposed_name
    
    
    conflict_subsid <-
      !is_missing(proposed$current_subsid) &
      !is_missing(proposed$proposed_subsid) &
      proposed$current_subsid != proposed$proposed_subsid
    
    
    conflict_smiles <-
      !is_missing(proposed$current_smiles) &
      !is_missing(proposed$proposed_smiles) &
      proposed$current_smiles != proposed$proposed_smiles
    
    
    # Build conflict table
    
    conflicts_list <- list()
    
    
    # NAME
    if (any(conflict_name)) {
      
      conflicts_list[["name"]] <- data.frame(
        database = db_name,
        compound_id =
          proposed$compound_id[conflict_name],
        column = "name",
        existing_value =
          proposed$current_name[conflict_name],
        proposed_value =
          proposed$proposed_name[conflict_name],
        stringsAsFactors = FALSE
      )
    }
    
    
    # SUBSID
    if (any(conflict_subsid)) {
      
      conflicts_list[["subsid"]] <- data.frame(
        database = db_name,
        compound_id =
          proposed$compound_id[conflict_subsid],
        column = "subsid",
        existing_value =
          proposed$current_subsid[conflict_subsid],
        proposed_value =
          proposed$proposed_subsid[conflict_subsid],
        stringsAsFactors = FALSE
      )
    }
    
    
    # SMILES
    if (any(conflict_smiles)) {
      
      conflicts_list[["smiles"]] <- data.frame(
        database = db_name,
        compound_id =
          proposed$compound_id[conflict_smiles],
        column = "smiles",
        existing_value =
          proposed$current_smiles[conflict_smiles],
        proposed_value =
          proposed$proposed_smiles[conflict_smiles],
        stringsAsFactors = FALSE
      )
    }
    
    
    # Combine conflicts
    if (length(conflicts_list) > 0) {
      
      conflicts <- do.call(
        rbind,
        conflicts_list
      )
      
      rownames(conflicts) <- NULL
      
    } else {
      
      conflicts <- data.frame(
        database = character(),
        compound_id = integer(),
        column = character(),
        existing_value = character(),
        proposed_value = character(),
        stringsAsFactors = FALSE
      )
    }
    
    
    # Prepare updates
    
    to_update <- proposed[
      needs_update,
      ,
      drop = FALSE
    ]
    
    
    n_values_added <-
      sum(fill_name) +
      sum(fill_subsid) +
      sum(fill_smiles)
    
    

    # If Nothing to update
    
    if (nrow(to_update) == 0) {
      
      cat(
        "Updated 0 compound(s) in",
        db_name,
        "\n"
      )
      
      cat(
        "Added 0 previously missing value(s).\n"
      )
      
      if (nrow(conflicts) > 0) {
        cat(
          "Found",
          nrow(conflicts),
          "existing value(s) that were NOT overwritten.\n"
        )
      }
      
      return(conflicts)
    }
    

    
    to_update$final_name <-
      to_update$current_name
    
    to_update$final_subsid <-
      to_update$current_subsid
    
    to_update$final_smiles <-
      to_update$current_smiles
    
    
    # Fill missing name
    x <-
      is_missing(to_update$final_name) &
      !is_missing(to_update$proposed_name)
    
    to_update$final_name[x] <-
      to_update$proposed_name[x]
    
    
    # Fill missing subsid
    x <-
      is_missing(to_update$final_subsid) &
      !is_missing(to_update$proposed_subsid)
    
    to_update$final_subsid[x] <-
      to_update$proposed_subsid[x]
    
    
    # Fill missing smiles
    x <-
      is_missing(to_update$final_smiles) &
      !is_missing(to_update$proposed_smiles)
    
    to_update$final_smiles[x] <-
      to_update$proposed_smiles[x]
    
  
    
    updates <- data.frame(
      compound_id =
        as.integer(to_update$compound_id),
      name =
        to_update$final_name,
      subsid =
        to_update$final_subsid,
      smiles =
        to_update$final_smiles,
      stringsAsFactors = FALSE
    )
    
    
    
    DBI::dbBegin(con)
    
    
    tryCatch({
      
      # Remove old temporary table if necessary
      DBI::dbExecute(
        con,
        "
        DROP TABLE IF EXISTS temp_metadata_updates
        "
      )
      
      
      # Write all updates at once
      DBI::dbWriteTable(
        con,
        "temp_metadata_updates",
        updates,
        temporary = TRUE,
        overwrite = TRUE,
        row.names = FALSE
      )
      
      
      # Index temporary table
      DBI::dbExecute(
        con,
        "
        CREATE INDEX idx_temp_metadata_compound
        ON temp_metadata_updates(compound_id)
        "
      )
      
      
      DBI::dbExecute(
        con,
        "
        UPDATE ms_compound
        SET
          name = (
            SELECT u.name
            FROM temp_metadata_updates AS u
            WHERE u.compound_id = ms_compound.compound_id
          ),
          
          subsid = (
            SELECT u.subsid
            FROM temp_metadata_updates AS u
            WHERE u.compound_id = ms_compound.compound_id
          ),
          
          smiles = (
            SELECT u.smiles
            FROM temp_metadata_updates AS u
            WHERE u.compound_id = ms_compound.compound_id
          )
          
        WHERE compound_id IN (
          SELECT compound_id
          FROM temp_metadata_updates
        )
        "
      )
      
      
      DBI::dbCommit(con)
      
      
    }, error = function(e) {
      
      if (DBI::dbIsValid(con)) {
        
        try(
          DBI::dbRollback(con),
          silent = TRUE
        )
      }
      
      stop(e)
    })
    
    
    
    cat(
      "Updated",
      nrow(to_update),
      "compound(s) in",
      db_name,
      "\n"
    )
    
    
    cat(
      "Added",
      n_values_added,
      "previously missing value(s).\n"
    )
    
    
    if (nrow(conflicts) > 0) {
      
      cat(
        "Found",
        nrow(conflicts),
        "existing value(s) that were NOT overwritten.\n"
      )
    }
    
    
    return(conflicts)
  }
  
  
  # Process reference database
  
  conflicts_ref <- update_database(
    ref_db,
    Updated_Names[[1]]
  )
  
  
  # Process target database
  
  conflicts_target <- update_database(
    target_db,
    Updated_Names[[2]]
  )
  
  
  # Combine conflicts
  
  conflicts <- rbind(
    conflicts_ref,
    conflicts_target
  )
  
  rownames(conflicts) <- NULL
  
  
  # Display conflicts
  
  if (nrow(conflicts) > 0) {
    
    cat(
      "\nExisting values were NOT overwritten.\n",
      "The following differences require manual review:\n\n"
    )
    
    print(
      conflicts,
      row.names = FALSE
    )
    
  } else {
    
    cat(
      "\nNo conflicts with existing metadata were found.\n"
    )
  }
  
  
  # Return results
  
  invisible(
    list(
      Mult.frame = Mult.frame,
      conflicts = conflicts
    )
  )
}