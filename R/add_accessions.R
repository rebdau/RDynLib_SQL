#' @title add_accessions()
#'
#' @details 
#' this function adds unique identifier for each experiment based on the user 
#' initials, the date stored in experiment table, and random chain.
#' If the user enter a given expid then the function will add those identifiers 
#' just for that experiment other ways it processed the whole database;
#' The function adds as well identifiers for the compounds using the experiment 
#' accessor concatenated with the compound_id, the compound accessors aims 
#' to keep track of the original compound_ids.
#' 
#' @param  dbfile path to the sql database file
#' @param expid experiment id, if NULL the accessors will be created for the 
#' whole database
#' @param initials_overrides initials of the user if Null, the first two letters
#'        from the user column will be used as initials
#' @example 
#' add_accessions(dbfile, 
#'                initials_overrides = c("266" = "RD",
#'                                    "267" = "RD"))
#' 
#' @returns accessors added to the experiment and ms_compound table
#' 
#' @author Ahlam Mentag
#' @export 
add_accessions <- function(
    dbfile = NULL,
    con = NULL,
    expid = NULL,
    initials_overrides = NULL,
    random_length = 6) {
  
  close_connection <- FALSE
  
  if (is.null(con)) {
    
    if (is.null(dbfile)) {
      stop("Provide either 'dbfile' or 'con'.")
    }
    
    con <- DBI::dbConnect(
      RSQLite::SQLite(),
      dbfile
    )
    
    close_connection <- TRUE
  }
  
  if (close_connection) {
    on.exit(
      DBI::dbDisconnect(con),
      add = TRUE
    )
  }
  
  
  ## Check required tables
  
  required_tables <- c(
    "experiment",
    "ms_compound"
  )
  
  missing_tables <- setdiff(
    required_tables,
    DBI::dbListTables(con)
  )
  
  if (length(missing_tables) > 0L) {
    stop(
      "Missing table(s): ",
      paste(
        missing_tables,
        collapse = ", "
      )
    )
  }
  
  
  ## Add experiment_accession column
  
  experiment_fields <- DBI::dbListFields(
    con,
    "experiment"
  )
  
  if (!"experiment_accession" %in%
      experiment_fields) {
    
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
  
  if (!"compound_accession" %in%
      compound_fields) {
    
    DBI::dbExecute(
      con,
      "
      ALTER TABLE ms_compound
      ADD COLUMN compound_accession TEXT
      "
    )
  }
  
  
  ## Helper function for random chain
  
  generate_random_chain <- function(n) {
    
    paste(
      sample(
        c(LETTERS, 0:9),
        n,
        replace = TRUE
      ),
      collapse = ""
    )
  }
  
  
  
  if (is.null(expid)) {
    
    experiments <- DBI::dbGetQuery(
      con,
      "
      SELECT
        expid,
        date,
        user,
        experiment_accession
      FROM experiment
      ORDER BY CAST(expid AS INTEGER)
      "
    )
    
  } else {
    
    experiments <- DBI::dbGetQuery(
      con,
      "
      SELECT
        expid,
        date,
        user,
        experiment_accession
      FROM experiment
      WHERE expid = ?
      ",
      params = list(expid)
    )
  }
  
  
  if (nrow(experiments) == 0L) {
    
    if (is.null(expid)) {
      stop("No experiments found in the database.")
    } else {
      stop(
        "Experiment ",
        expid,
        " was not found."
      )
    }
  }
  
  
  ## Generate accessions
  
  for (i in seq_len(nrow(experiments))) {
    
    current_expid <-
      experiments$expid[i]
    
    experiment_date <-
      experiments$date[i]
    
    experiment_user <-
      experiments$user[i]
    
    current_accession <-
      experiments$experiment_accession[i]
    
    
    ## If accession already exists, keep it
    
    if (!is.na(current_accession) &&
        nzchar(trimws(current_accession))) {
      
      experiment_accession <-
        current_accession
      
    } else {
      
      
      ## Read date FROM experiment table
      
      if (is.na(experiment_date) ||
          !nzchar(trimws(experiment_date))) {
        
        stop(
          "Missing date for experiment ",
          current_expid,
          ". Cannot generate accession."
        )
      }
      
      
      experiment_date <-
        trimws(
          as.character(
            experiment_date
          )
        )
      
      
      ## Expected database format:
      ## DD/MM/YYYY
      
      parsed_date <- as.Date(
        experiment_date,
        format = "%d/%m/%Y"
      )
      
      
      if (is.na(parsed_date)) {
        
        stop(
          "Could not parse date '",
          experiment_date,
          "' for experiment ",
          current_expid,
          ". Expected DD/MM/YYYY."
        )
      }
      
      
      accession_date <- format(
        parsed_date,
        "%Y%m%d"
      )
      
      
      ## Determine initials
      
      current_initials <- NULL
      
      
      ## If initials are explicitly provided
      if (!is.null(initials_overrides)) {
        
 
        
        if (length(initials_overrides) == 1L &&
            is.null(names(initials_overrides))) {
          
          current_initials <- initials_overrides
          
        } else {
          
  
          
          override_names <- names(
            initials_overrides
          )
          
          if (!is.null(override_names) &&
              as.character(current_expid) %in%
              override_names) {
            
            current_initials <-
              initials_overrides[
                as.character(current_expid)
              ]
          }
        }
      }
      
      
      ## If initials were not provided,
      ## derive them from experiment$user
      
      if (is.null(current_initials)) {
        
        if (is.na(experiment_user) ||
            !nzchar(trimws(experiment_user))) {
          
          stop(
            "Missing user for experiment ",
            current_expid,
            ". Cannot determine initials."
          )
        }
        
        current_initials <- substr(
          trimws(experiment_user),
          1,
          2
        )
      }
      
      
      ## Clean initials
      
      current_initials <- toupper(
        trimws(
          as.character(
            current_initials
          )
        )
      )
      
      
      if (!nzchar(current_initials)) {
        
        stop(
          "Could not determine initials for experiment ",
          current_expid,
          "."
        )
      }
      
      
      ## Generate unique experiment accession
      
      repeat {
        
        random_chain <-
          generate_random_chain(
            random_length
          )
        
        
        experiment_accession <- paste(
          current_initials,
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
      
      
      ## Store experiment accession
      
      DBI::dbExecute(
        con,
        "
        UPDATE experiment
        SET experiment_accession = ?
        WHERE expid = ?
        ",
        params = list(
          experiment_accession,
          current_expid
        )
      )
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
        current_expid
      )
    )
    
    
    message(
      "Experiment ",
      current_expid,
      " | date = ",
      experiment_date,
      " | accession = ",
      experiment_accession
    )
  }
  
  
  ## Unique indexes

  
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
    "Accessions successfully added."
  )
  
  
  invisible(TRUE)
}