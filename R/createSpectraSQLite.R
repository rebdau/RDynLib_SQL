#' @title Create new SQL database from a given spectra object
#'
#' @details
#' From a Spectra object, createSpectraSQLite() creates an SQLite database
#' containing the spectral data of the object given as input.
#'
#' The resulting SQL structure is compatible with 'MsBackendCompDb'
#' from the 'RforMassSpectrometry' ecosystem, allowing the SQL data
#' to be converted easily to a Spectra object.
#'
#' An experiment accession is automatically generated using the date on which
#' the experiment is added to the database:
#'
#'   USER_INITIALS-YYYYMMDD-RANDOM_CHAIN
#'
#' For example:
#'
#'   RD-20260907-Y1D6K8
#'
#' The experimental date stored in the `date` column is provided separately
#' by the user and does not affect the experiment accession.
#'
#' Each compound receives a compound accession combining the experiment
#' accession and its compound identifier:
#'
#'   RD-20260907-Y1D6K8-101280
#' If  No expid is passed: add_accessions() checks the WHOLE database.
#'
#' @param sps A Spectra object containing the spectra to store.
#'
#' @param dbfile `character(1)` Path to the output SQLite database.
#'
#' @param date `character(1)` Experimental acquisition date provided by
#' the user.
#'
#' @param user `character(1)` User associated with the experiment.
#'
#' @param user_initials `character(1)` User initials used to generate the
#' experiment accession.
#'
#' @param machine `character(1)` Instrument used for acquisition.
#'
#' @param mode `character(1)` Ionisation mode.
#'
#' @param tissue `character(1)` Tissue or sample material.
#'
#' @param mstype, column, buffera, bufferb, gradient_time, source, species,
#' ce, meta Optional experimental metadata describing the acquisition
#' conditions, sample origin, and experimental setup.
#'
#' @return Invisibly returns the newly created experiment ID.
#'
#' @importFrom DBI dbGetQuery
#' @importFrom DBI dbConnect
#' @importFrom DBI dbDisconnect
#' @importFrom DBI dbExecute
#' @importFrom DBI dbWriteTable
#' @importFrom RSQLite SQLite
#' @import tidyr
#' @import tibble
#' @import dplyr
#'
#' @author Ahlam Mentag
#'
#' @export
createSpectraSQLite <- function(
    sps,
    dbfile,
    date,
    user,
    user_initials,
    machine,
    mode,
    tissue,
    mstype = NULL,
    column = NULL,
    buffera = NULL,
    bufferb = NULL,
    gradient_time = NULL,
    source = NULL,
    species = NULL,
    ce = NULL,
    meta = NULL) {
  
  
  
  con <- DBI::dbConnect(
    RSQLite::SQLite(),
    dbfile
  )
  
  on.exit(
    DBI::dbDisconnect(con),
    add = TRUE
  )
  
  
  
  safe_scalar <- function(x) {
    
    if (is.null(x)) {
      NA_character_
    } else {
      as.character(x)
    }
  }
  
  
  
  
  user_initials <- toupper(
    trimws(
      user_initials
    )
  )
  
  
  if (length(user_initials) != 1L ||
      is.na(user_initials) ||
      !nzchar(user_initials)) {
    
    stop(
      "'user_initials' must contain one non-empty value."
    )
  }
  
  
  ## TABLE experiment
  
  
  DBI::dbExecute(
    con,
    "
    CREATE TABLE IF NOT EXISTS experiment (

      expid INTEGER PRIMARY KEY AUTOINCREMENT,

      experiment_accession TEXT UNIQUE,

      date TEXT,
      user TEXT,
      machine TEXT,
      column TEXT,
      mstype TEXT,
      buffera TEXT,
      bufferb TEXT,
      gradient_time TEXT,
      source TEXT,
      mode TEXT,
      species TEXT,
      tissue TEXT,
      ce TEXT,
      meta TEXT

    );
    "
  )
  
  
  ## Add experiment
  
  
  experiment_df <- data.frame(
    
    experiment_accession =
      NA_character_,
    
    date =
      safe_scalar(date),
    
    user =
      safe_scalar(user),
    
    machine =
      safe_scalar(machine),
    
    column =
      safe_scalar(column),
    
    mstype =
      safe_scalar(mstype),
    
    buffera =
      safe_scalar(buffera),
    
    bufferb =
      safe_scalar(bufferb),
    
    gradient_time =
      safe_scalar(gradient_time),
    
    source =
      safe_scalar(source),
    
    mode =
      safe_scalar(mode),
    
    species =
      safe_scalar(species),
    
    tissue =
      safe_scalar(tissue),
    
    ce =
      safe_scalar(ce),
    
    meta =
      safe_scalar(meta),
    
    stringsAsFactors = FALSE
  )
  
  
  DBI::dbWriteTable(
    con,
    "experiment",
    experiment_df,
    append = TRUE,
    row.names = FALSE
  )
  
  
  
  ## Get newly created expid
  
  
  expid <- DBI::dbGetQuery(
    con,
    "
    SELECT last_insert_rowid() AS expid
    "
  )$expid[[1]]
  
  
  ## TABLE ms_compound
  
  DBI::dbExecute(
    con,
    "
    CREATE TABLE IF NOT EXISTS ms_compound (

      compound_id TEXT PRIMARY KEY,

      compound_accession TEXT UNIQUE,

      expid INTEGER,

      nodename TEXT,
      retention_time REAL,
      mass_measured REAL,
      name TEXT,
      formula TEXT,
      exactmass REAL,
      ppm_deviation REAL,
      subsid INTEGER,
      conversion TEXT,
      wavelen REAL,
      smiles TEXT,
      isotope_ratio REAL,
      drift_time REAL,
      composition TEXT,
      inchi TEXT,
      inchikey TEXT,

      FOREIGN KEY (expid)
        REFERENCES experiment(expid)

    );
    "
  )
  
  
  ## Extract compound information
  
  sp_data <- Spectra::spectraData(
    sps
  )
  
  sp_df <- as.data.frame(
    sp_data
  )
  
  
  cmp <- sp_df |>
    dplyr::select(
      feature_id,
      feature_rtmed,
      feature_mzmed
    ) |>
    dplyr::filter(
      !is.na(feature_id),
      !is.na(feature_rtmed),
      !is.na(feature_mzmed)
    ) |>
    dplyr::distinct(
      feature_id,
      .keep_all = TRUE
    )
  
  
  n <- nrow(
    cmp
  )
  
  
  ## Generate compound IDs
  
  last_id <- DBI::dbGetQuery(
    con,
    "
    SELECT
      MAX(
        CAST(
          compound_id AS INTEGER
        )
      ) AS max_id
    FROM ms_compound
    "
  )$max_id[[1]]
  
  
  start_id <- if (is.na(last_id)) {
    
    1L
    
  } else {
    
    as.integer(
      last_id
    ) + 1L
  }
  
  
  compound_ids <- as.character(
    seq(
      from = start_id,
      length.out = n
    )
  )
  
  
  ## Accessions are created later by add_accessions()
  
  compound_accessions <- rep(
    NA_character_,
    n
  )
  
  
  ## Build compound table
  
  ms_compound_df <- data.frame(
    
    compound_id =
      compound_ids,
    
    compound_accession =
      compound_accessions,
    
    expid =
      expid,
    
    nodename =
      cmp$feature_id,
    
    retention_time =
      cmp$feature_rtmed,
    
    mass_measured =
      cmp$feature_mzmed,
    
    name =
      NA_character_,
    
    formula =
      NA_character_,
    
    exactmass =
      NA_real_,
    
    ppm_deviation =
      NA_real_,
    
    subsid =
      NA_integer_,
    
    conversion =
      NA_character_,
    
    wavelen =
      NA_real_,
    
    smiles =
      NA_character_,
    
    isotope_ratio =
      NA_real_,
    
    drift_time =
      NA_real_,
    
    composition =
      NA_character_,
    
    inchi =
      NA_character_,
    
    inchikey =
      NA_character_,
    
    stringsAsFactors = FALSE
  )
  
  
  DBI::dbWriteTable(
    con,
    "ms_compound",
    ms_compound_df,
    append = TRUE,
    row.names = FALSE
  )
  
  
  accession_overrides <- stats::setNames(
    user_initials,
    as.character(expid)
  )
  
  
  add_accessions(
    con = con,
    initials_overrides =
      accession_overrides
  )
  
  
  ## Retrieve accession of newly created experimen
  
  experiment_accession <- DBI::dbGetQuery(
    con,
    "
    SELECT experiment_accession
    FROM experiment
    WHERE expid = ?
    ",
    params = list(expid)
  )$experiment_accession[[1]]
  
  
  ## TABLE msms_spectrum
  
  DBI::dbExecute(
    con,
    "
    CREATE TABLE IF NOT EXISTS msms_spectrum (

      spectrum_id INTEGER
        PRIMARY KEY AUTOINCREMENT,

      compound_id TEXT,

      ms_level INTEGER,
      polarity INTEGER,
      spectrum_type TEXT,
      precursor_mz REAL,
      precursorIntensity REAL,
      precursorCharge INTEGER,
      collision_energy TEXT,
      isolationWindowLowerMz REAL,
      isolationWindowTargetMz REAL,
      isolationWindowUpperMz REAL,
      peaks_count INTEGER,
      totIonCurrent REAL,
      basePeakMZ REAL,
      basePeakIntensity REAL,
      ionisationEnergy REAL,
      lowMZ REAL,
      highMZ REAL,
      mergedScan INTEGER,
      mergedResultScanNum INTEGER,
      mergedResultStartScanNum INTEGER,
      mergedResultEndScanNum INTEGER,
      injectionTime REAL,
      filterString TEXT,
      spectrumId INTEGER,
      ionMobilityDriftTime REAL,
      scanWindowLowerLimit REAL,
      scanWindowUpperLimit REAL,
      electronBeamEnergy REAL,
      originalPrecursorMz REAL,
      precursorPurity REAL,
      chromPeakRT REAL,
      chromPeakMz REAL,
      chromPeakId TEXT,
      rtime REAL,
      scanIndex INTEGER,
      dataStorage TEXT,
      centroided INTEGER,
      smoothed INTEGER,
      instrument TEXT,
      splash TEXT,
      instrument_type TEXT,
      acquisitionNum INTEGER,
      precScanNum INTEGER,
      predicted REAL,
      dataOrigin TEXT,
      original_id TEXT,

      FOREIGN KEY (compound_id)
        REFERENCES ms_compound(compound_id)

    );
    "
  )
  
  
  ## Spectra data
  
  sp <- Spectra::spectraData(
    sps
  )
  
  
  getcol <- function(x) {
    
    if (x %in% names(sp)) {
      
      sp[[x]]
      
    } else {
      
      rep(
        NA,
        nrow(sp)
      )
    }
  }
  
  
  msms_df <- data.frame(
    
    compound_id =
      NA_character_,
    
    ms_level =
      getcol("msLevel"),
    
    polarity =
      getcol("polarity"),
    
    spectrum_type =
      getcol("spectrum.type"),
    
    precursor_mz =
      getcol("precursorMz"),
    
    precursorIntensity =
      getcol("precursorIntensity"),
    
    precursorCharge =
      getcol("precursorCharge"),
    
    collision_energy =
      getcol("collisionEnergy"),
    
    isolationWindowLowerMz =
      getcol("isolationWindowLowerMz"),
    
    isolationWindowTargetMz =
      getcol("isolationWindowTargetMz"),
    
    isolationWindowUpperMz =
      getcol("isolationWindowUpperMz"),
    
    peaks_count =
      getcol("peaksCount"),
    
    totIonCurrent =
      getcol("totIonCurrent"),
    
    basePeakMZ =
      getcol("basePeakMZ"),
    
    basePeakIntensity =
      getcol("basePeakIntensity"),
    
    ionisationEnergy =
      getcol("ionisationEnergy"),
    
    lowMZ =
      getcol("lowMZ"),
    
    highMZ =
      getcol("highMZ"),
    
    mergedScan =
      getcol("mergedScan"),
    
    mergedResultScanNum =
      getcol("mergedResultScanNum"),
    
    mergedResultStartScanNum =
      getcol(
        "mergedResultStartScanNum"
      ),
    
    mergedResultEndScanNum =
      getcol(
        "mergedResultEndScanNum"
      ),
    
    injectionTime =
      getcol("injectionTime"),
    
    filterString =
      getcol("filterString"),
    
    spectrumId =
      getcol("spectrumId"),
    
    ionMobilityDriftTime =
      getcol(
        "ionMobilityDriftTime"
      ),
    
    scanWindowLowerLimit =
      getcol(
        "scanWindowLowerLimit"
      ),
    
    scanWindowUpperLimit =
      getcol(
        "scanWindowUpperLimit"
      ),
    
    electronBeamEnergy =
      getcol(
        "electronBeamEnergy"
      ),
    
    originalPrecursorMz =
      getcol(
        "originalPrecursorMz"
      ),
    
    precursorPurity =
      getcol("precursorPurity"),
    
    chromPeakRT =
      getcol("chromPeakRT"),
    
    chromPeakMz =
      getcol("chromPeakMz"),
    
    chromPeakId =
      getcol("chromPeakId"),
    
    rtime =
      getcol("rtime"),
    
    scanIndex =
      getcol("scanIndex"),
    
    dataStorage =
      getcol("dataStorage"),
    
    centroided =
      getcol("centroided"),
    
    smoothed =
      getcol("smoothed"),
    
    instrument =
      getcol("instrument"),
    
    splash =
      getcol("splash"),
    
    instrument_type =
      getcol("instrument_type"),
    
    acquisitionNum =
      getcol("acquisitionNum"),
    
    precScanNum =
      getcol("precScanNum"),
    
    predicted =
      getcol("predicted"),
    
    dataOrigin =
      getcol("dataOrigin"),
    
    original_id =
      getcol("original_id"),
    
    stringsAsFactors = FALSE
  )
  
  
  ## Link spectra to compounds
  
  msms_df$compound_id <-
    ms_compound_df$compound_id[
      match(
        sp$feature_id,
        ms_compound_df$nodename
      )
    ]
  
  
  ## Determine first spectrum ID
  
  last_spectrum_id <- DBI::dbGetQuery(
    con,
    "
    SELECT MAX(spectrum_id) AS max_id
    FROM msms_spectrum
    "
  )$max_id[[1]]
  
  
  if (is.na(last_spectrum_id)) {
    
    first_spectrum_id <- 1L
    
  } else {
    
    first_spectrum_id <-
      as.integer(
        last_spectrum_id
      ) + 1L
  }
  
  
  ## Store spectra
  
  DBI::dbWriteTable(
    con,
    "msms_spectrum",
    msms_df,
    append = TRUE,
    row.names = FALSE
  )
  
  
  spectrum_ids <- seq(
    from =
      first_spectrum_id,
    length.out =
      nrow(msms_df)
  )
  
  
  ## TABLE msms_spectrum_peak
  
  DBI::dbExecute(
    con,
    "
    CREATE TABLE IF NOT EXISTS msms_spectrum_peak (

      peak_id INTEGER
        PRIMARY KEY AUTOINCREMENT,

      spectrum_id INTEGER,

      mz REAL,
      intensity REAL,

      FOREIGN KEY (spectrum_id)
        REFERENCES msms_spectrum(spectrum_id)

    );
    "
  )
  
  
  ## Peak data
  
  pd <- Spectra::peaksData(
    sps
  )
  
  
  peaks_list <- lapply(
    seq_along(pd),
    function(i) {
      
      if (!is.null(pd[[i]]) &&
          nrow(pd[[i]]) > 0L) {
        
        data.frame(
          
          spectrum_id =
            spectrum_ids[[i]],
          
          mz =
            pd[[i]][, 1],
          
          intensity =
            pd[[i]][, 2]
        )
        
      } else {
        
        NULL
      }
    }
  )
  
  
  peaks_list <- Filter(
    Negate(is.null),
    peaks_list
  )
  
  
  if (length(peaks_list) > 0L) {
    
    peaks_df <- do.call(
      rbind,
      peaks_list
    )
    
    
    DBI::dbWriteTable(
      con,
      "msms_spectrum_peak",
      peaks_df,
      append = TRUE,
      row.names = FALSE
    )
  }
  
  
  ## TABLE Synonym
  
  DBI::dbExecute(
    con,
    "
    CREATE TABLE IF NOT EXISTS Synonym (
      id INTEGER PRIMARY KEY
    );
    "
  )
  
  
  ## TABLE metadata
  
  if (exists("make_metadata")) {
    
    md <- make_metadata(
      
      source =
        "Flax FTMS neg",
      
      url =
        NA_character_,
      
      source_version =
        "1.0.0",
      
      source_date =
        as.character(
          Sys.Date()
        )
    )
    
    
    DBI::dbWriteTable(
      con,
      "metadata",
      md,
      overwrite = TRUE,
      row.names = FALSE
    )
  }
  
  
  message(
    "Experiment created: ",
    experiment_accession
  )
  
  
  invisible(expid)
}