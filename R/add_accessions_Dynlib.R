add_accessions_to_db <- function(dbfile) {
  
  con <- DBI::dbConnect(RSQLite::SQLite(), dbfile)
  on.exit(DBI::dbDisconnect(con), add = TRUE)
  
  ## Check required tables
  required_tables <- c("experiment", "ms_compound")
  
  missing_tables <- setdiff(
    required_tables,
    DBI::dbListTables(con)
  )
  
  if (length(missing_tables) > 0) {
    stop(
      "Missing table(s): ",
      paste(missing_tables, collapse = ", ")
    )
  }
  
  ## Add experiment_accession column if necessary
  experiment_fields <- DBI::dbListFields(
    con,
    "experiment"
  )
  
  if (!"experiment_accession" %in% experiment_fields) {
    
    DBI::dbExecute(
      con,
      "
      ALTER TABLE experiment
      ADD COLUMN experiment_accession TEXT
      "
    )
  }
  
  ## Add compound_accession column if necessary
  compound_fields <- DBI::dbListFields(
    con,
    "ms_compound"
  )
  
  if (!"compound_accession" %in% compound_fields) {
    
    DBI::dbExecute(
      con,
      "
      ALTER TABLE ms_compound
      ADD COLUMN compound_accession TEXT
      "
    )
  }
  
  ## Helper function for random chain
  generate_random_chain <- function(n = 6) {
    
    paste(
      sample(
        c(LETTERS, 0:9),
        n,
        replace = TRUE
      ),
      collapse = ""
    )
  }
  
  
  ## Get experiments, including user
  experiments <- DBI::dbGetQuery(
    con,
    "
    SELECT expid, user, date, experiment_accession
    FROM experiment
    ORDER BY CAST(expid AS INTEGER)
    "
  )
  
  ## Generate accession for each experiment
  for (i in seq_len(nrow(experiments))) {
    
    expid <- experiments$expid[i]
    user <- experiments$user[i]
    experiment_date <- experiments$date[i]
    
    current_accession <-
    experiments$experiment_accession[i]
    ## Convert experiment date from DD/MM/YYYY to YYYYMMDD
    parsed_date <- as.Date(
      experiment_date,
      format = "%d/%m/%Y"
    )
    
    if (is.na(parsed_date)) {
      stop(
        "Could not parse date '",
        experiment_date,
        "' for experiment ",
        expid
      )
    }
    
    accession_date <- format(
      parsed_date,
      "%Y%m%d"
    )
    
    ## Determine user initials
    
    if (expid %in% c(266, 267)) {
      
      user_initials <- "RD"
      
    } else {
      
      if (is.na(user) || !nzchar(trimws(user))) {
        stop(
          "Missing user for experiment ",
          expid,
          ". Cannot generate initials."
        )
      }
      
      ## Remove leading/trailing spaces
      user_clean <- trimws(user)
      
      ## Take first two characters and convert to upper case
      user_initials <- toupper(
        substr(user_clean, 1, 2)
      )
    }
    
    ## Only generate accession if missing
    
    if (is.na(current_accession) ||
        current_accession == "") {
      
      repeat {
        
        random_chain <- generate_random_chain()
        
        experiment_accession <- paste(
          user_initials,
          accession_date,
          random_chain,
          sep = "-"
        )
        
        existing <- DBI::dbGetQuery(
          con,
          "
          SELECT expid
          FROM experiment
          WHERE experiment_accession = ?
          ",
          params = list(
            experiment_accession
          )
        )
        
        if (nrow(existing) == 0L) {
          break
        }
      }
      
      ## Update experiment
      DBI::dbExecute(
        con,
        "
        UPDATE experiment
        SET experiment_accession = ?
        WHERE expid = ?
        ",
        params = list(
          experiment_accession,
          expid
        )
      )
      
    } else {
      
      ## Keep existing accession
      experiment_accession <-
        current_accession
    }
    
    ## Generate compound accessions
    
    DBI::dbExecute(
      con,
      "
      UPDATE ms_compound
      SET compound_accession =
          ? || '-' || compound_id
      WHERE expid = ?
        AND (
          compound_accession IS NULL
          OR compound_accession = ''
        )
      ",
      params = list(
        experiment_accession,
        expid
      )
    )
    
    message(
      "Experiment ",
      expid,
      " | user = ",
      user,
      " | initials = ",
      user_initials,
      " | accession = ",
      experiment_accession
    )
  }
  
  ## Add unique indexes
  DBI::dbExecute(
    con,
    "
    CREATE UNIQUE INDEX IF NOT EXISTS
    idx_experiment_accession
    ON experiment(experiment_accession)
    "
  )
  
  DBI::dbExecute(
    con,
    "
    CREATE UNIQUE INDEX IF NOT EXISTS
    idx_compound_accession
    ON ms_compound(compound_accession)
    "
  )
  
  message(
    "Accessions successfully added to database."
  )
  
  invisible(TRUE)
}