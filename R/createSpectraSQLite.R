#' @title Create new SQL database from a given spectra object
#'
#' @details
#' From a spectra object createSpectraSQLite() function creates an sqlite 
#' database with the spectra data of the object given in the input of the 
#' function, the resulting sql structure is compatible with 'MsBackendCompDb'
#' from 'RforMassSpecrtrometry' ecosystem, this allows to convert the sql 
#' data easily to a spectra object.  
#' 
#' @param sps `character(1)` the path to the spectra object to convert.
#'
#' @param dbfile `character(1)` the path to the output sql database.
#' 
#' @param date, user, machine, mode, tissue, mstype, column, buffera,
#' bufferb, gradient_time, source, species, ce, meta
#' Experimental metadata provided by the user. These fields describe
#' the acquisition conditions, sample origin, and experimental setup
#' associated with the spectra.
#'
#'
#' @return SQL database containing the spectral information of the given object
#'
#' @importFrom DBI dbGetQuery
#'
#' @importFrom DBI dbConnect
#'
#' @importFrom DBI dbDisconnect
#'
#' @importFrom RSQLite SQLite
#'
#' @import tidyr
#'
#' @import tibble
#' 
#' @import dplyr
#' 
#' @author Ahlam Mentag
#'
#' @export
createSpectraSQLite <- function(sps, dbfile, date, user, machine, mode, tissue,
                                mstype = NULL, column = NULL, buffera = NULL,
                                bufferb = NULL, gradient_time = NULL, 
                                source = NULL, species = NULL, ce = NULL,
                                meta = NULL) {
  
  con <- dbConnect(SQLite(), dbfile)
  on.exit(dbDisconnect(con))
  
  ## TABLE experiment
  
  dbExecute(con, "
    CREATE TABLE IF NOT EXISTS experiment (
      expid INTEGER PRIMARY KEY AUTOINCREMENT,
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
  ")
  
  safe_scalar <- function(x) {
    if (is.null(x)) NA_character_ else as.character(x)
  }
  
  experiment_df <- data.frame(
    date = safe_scalar(date),
    user = safe_scalar(user),
    machine = safe_scalar(machine),
    column = safe_scalar(column),
    mstype = safe_scalar(mstype),
    buffera = safe_scalar(buffera),
    bufferb = safe_scalar(bufferb),
    gradient_time = safe_scalar(gradient_time),
    source = safe_scalar(source),
    mode = safe_scalar(mode),
    species = safe_scalar(species),
    tissue = safe_scalar(tissue),
    ce = safe_scalar(ce),
    meta = safe_scalar(meta),
    stringsAsFactors = FALSE
  )
  
  dbWriteTable(con, "experiment", experiment_df, append = TRUE, row.names = FALSE)
  
  expid <- dbGetQuery(con, "
    SELECT last_insert_rowid() AS expid
  ")$expid
  
  ## TABLE ms_compound
  
  dbExecute(con, "
    CREATE TABLE IF NOT EXISTS ms_compound (
      compound_id TEXT PRIMARY KEY,
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
      FOREIGN KEY (expid) REFERENCES experiment(expid)
    );
  ")
  
  sp_data <- spectraData(sps)
  sp_df <- as.data.frame(sp_data)
  
  cmp <- sp_df %>%
    select(feature_id, feature_rtmed, feature_mzmed) %>%
    filter(!is.na(feature_id),
           !is.na(feature_rtmed),
           !is.na(feature_mzmed)) %>%
    distinct(feature_id, .keep_all = TRUE)
  
  n <- nrow(cmp)
  
  last_id <- dbGetQuery(con, "
    SELECT MAX(CAST(compound_id AS INTEGER)) AS max_id
    FROM ms_compound
  ")$max_id
  
  start_id <- if (is.na(last_id)) 1 else last_id + 1
  
  ms_compound_df <- data.frame(
    compound_id = as.character(seq(start_id, length.out = n)),
    expid = expid,
    nodename = cmp$feature_id,
    retention_time = cmp$feature_rtmed,
    mass_measured = cmp$feature_mzmed,
    name = NA_character_,
    formula = NA_character_,
    exactmass = NA_real_,
    ppm_deviation = NA_real_,
    subsid = NA_integer_,
    conversion = NA_character_,
    wavelen = NA_real_,
    smiles = NA_character_,
    isotope_ratio = NA_real_,
    drift_time = NA_real_,
    composition = NA_character_,
    inchi = NA_character_,
    inchikey = NA_character_,
    stringsAsFactors = FALSE
  )
  
  dbWriteTable(con, "ms_compound", ms_compound_df, append = TRUE, row.names = FALSE)
  
  ## TABLE msms_spectrum
  
  dbExecute(con, "DROP TABLE IF EXISTS msms_spectrum")
  
  dbExecute(con, "
    CREATE TABLE msms_spectrum (
      spectrum_id INTEGER PRIMARY KEY AUTOINCREMENT,
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
  ")
  
  sp <- spectraData(sps)
  
  getcol <- function(x) {
    if (x %in% names(sp)) sp[[x]] else rep(NA, nrow(sp))
  }
  
  msms_df <- data.frame(
    compound_id = NA,
    ms_level = getcol("msLevel"),
    polarity = getcol("polarity"),
    spectrum_type = getcol("spectrum.type"),
    precursor_mz = getcol("precursorMz"),
    precursorIntensity = getcol("precursorIntensity"),
    precursorCharge = getcol("precursorCharge"),
    collision_energy = getcol("collisionEnergy"),
    isolationWindowLowerMz = getcol("isolationWindowLowerMz"),
    isolationWindowTargetMz = getcol("isolationWindowTargetMz"),
    isolationWindowUpperMz = getcol("isolationWindowUpperMz"),
    peaks_count = getcol("peaksCount"),
    totIonCurrent = getcol("totIonCurrent"),
    basePeakMZ = getcol("basePeakMZ"),
    basePeakIntensity = getcol("basePeakIntensity"),
    ionisationEnergy = getcol("ionisationEnergy"),
    lowMZ = getcol("lowMZ"),
    highMZ = getcol("highMZ"),
    mergedScan = getcol("mergedScan"),
    mergedResultScanNum = getcol("mergedResultScanNum"),
    mergedResultStartScanNum = getcol("mergedResultStartScanNum"),
    mergedResultEndScanNum = getcol("mergedResultEndScanNum"),
    injectionTime = getcol("injectionTime"),
    filterString = getcol("filterString"),
    spectrumId = getcol("spectrumId"),
    ionMobilityDriftTime = getcol("ionMobilityDriftTime"),
    scanWindowLowerLimit = getcol("scanWindowLowerLimit"),
    scanWindowUpperLimit = getcol("scanWindowUpperLimit"),
    electronBeamEnergy = getcol("electronBeamEnergy"),
    originalPrecursorMz = getcol("originalPrecursorMz"),
    precursorPurity = getcol("precursorPurity"),
    chromPeakRT = getcol("chromPeakRT"),
    chromPeakMz = getcol("chromPeakMz"),
    chromPeakId = getcol("chromPeakId"),
    rtime = getcol("rtime"),
    scanIndex = getcol("scanIndex"),
    dataStorage = getcol("dataStorage"),
    centroided = getcol("centroided"),
    smoothed = getcol("smoothed"),
    instrument = getcol("instrument"),
    splash = getcol("splash"),
    instrument_type = getcol("instrument_type"),
    acquisitionNum = getcol("acquisitionNum"),
    precScanNum = getcol("precScanNum"),
    predicted = getcol("predicted"),
    dataOrigin = getcol("dataOrigin"),
    original_id = getcol("original_id"),
    stringsAsFactors = FALSE
  )
  
  msms_df$compound_id <-
    ms_compound_df$compound_id[
      match(sp$feature_id, ms_compound_df$nodename)
    ]
  
  dbWriteTable(con, "msms_spectrum", msms_df, append = TRUE, row.names = FALSE)
  
  ## msms_spectrum_peak
  
  dbExecute(con, "
    CREATE TABLE IF NOT EXISTS msms_spectrum_peak (
      peak_id INTEGER PRIMARY KEY AUTOINCREMENT,
      spectrum_id INTEGER,
      mz REAL,
      intensity REAL,
      FOREIGN KEY (spectrum_id)
        REFERENCES msms_spectrum(spectrum_id)
    );
  ")
  
  pd <- peaksData(sps)
  
  peaks_df <- do.call(rbind, lapply(seq_along(pd), function(i) {
    if (!is.null(pd[[i]]) && nrow(pd[[i]]) > 0) {
      data.frame(
        spectrum_id = i,
        mz = pd[[i]][,1],
        intensity = pd[[i]][,2]
      )
    }
  }))
  
  dbWriteTable(con, "msms_spectrum_peak", peaks_df, append = TRUE, row.names = FALSE)
  

  ## Synonym

  
  dbExecute(con, "
    CREATE TABLE IF NOT EXISTS Synonym (
      id INTEGER PRIMARY KEY
    )
  ")
  
  ## metadata
  
  if (exists("make_metadata")) {
    md <- make_metadata(
      source = "Flax FTMS neg",
      url = NA_character_,
      source_version = "1.0.0",
      source_date = as.character(Sys.Date())
    )
    
    dbWriteTable(con, "metadata", md, overwrite = TRUE, row.names = FALSE)
  }
  
  invisible(expid)
}