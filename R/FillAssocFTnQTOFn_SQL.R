#' @title Associate FT and QTOF features using SQL and MS2 matching
#'
#' @description
#' `FillAssocFTnQTOFn_SQL()` finds associations between FTMS and QTOFMS
#' compounds based on retention time alignment, MS2 peak similarity, and
#' alignment rules. It updates and returns the association table (`Assoc`).
#'
#' @details
#' The function:
#' - loads FT and QTOF compounds from SQL databases,
#' - computes expected RT alignment using regression parameters,
#' - finds candidate matches using SQL-assisted search,
#' - filters candidates using MS2 peak similarity,
#' - appends validated associations to `Assoc`.
#'
#' @param FT_con `DBIConnection`  
#'   A DBI connection object to the FT SQLite database.
#'
#' @param QTOF_con `DBIConnection`  
#'   A DBI connection object to the QTOF SQLite database.
#'
#' @param Assoc `data.frame` A table storing previously detected associations. 
#'   Will be filled with the new association if they are not already 
#'   in the assoc table.
#'
#' @param FT_expnr `numeric(1)`  
#'   The reference experiment number to load from the FTMS SQL database.
#'
#' @param QTOF_expnr `numeric(1)`  
#'   Experiment number to load from the QTOF SQL database.
#'
#' @param cutoff `numeric(1)`  
#'   Minimum retention time to consider for FT compounds.
#'
#' @param rg `numeric(6)`  
#'   Regression coefficients for retention-time alignment.
#'
#' @param lc.err `numeric(1)`  
#'   Allowed retention-time error window after regression transformation.
#'
#' @param err `numeric(1)`  
#'   Error threshold for candidate match search (`Find_cand_matches_SQL`).
#'
#' @param minIon `numeric(1)`  
#'    The minimum matching ions between spectra of peak pairs, which are found
#'    after applying the regression model, and which are the final matching 
#'    written in the Assoc table. (default `0.02`).
#'
#' @param polarity_ft `integer(1)`  
#'   Polarity filter for FTMS MS2 spectra (0 or 1).
#'
#' @param polarity_qtof `integer(1)`  
#'   Polarity filter for QTOF MS2 spectra (0 or 1).
#'
#' @param FT_path `character(1)`  
#'   File path to the FTMS database .
#'
#' @param QTOF_path `character(1)`  
#'   File path to the QTOF database.
#'
#' @param spectrum_type_FT `character(1)` Spectrum type to use for the FT-MS
#'        spectra. If `NULL` and the database does not contain a `spectrum_type`
#'        column, all spectra are used.
#'
#' @param spectrum_type_QTOF `character(1)` Spectrum type to use for the QTOF
#'        spectra. If `NULL` and the database does not contain a `spectrum_type`
#'        column, all spectra are used.
#'
#' @return `data.frame`  
#'   The updated association table with new FTMS–QTOF matches appended.
#'
#' @importFrom DBI dbConnect
#'
#' @importFrom DBI dbDisconnect
#'
#' @importFrom RSQLite SQLite
#'
#' @author Ahlam Mentag
#'
#' @export

FillAssocFTnQTOFn_SQL <- function(
    FT_con,
    QTOF_con,
    Assoc,
    FT_expnr,
    QTOF_expnr,
    cutoff,
    rg,
    lc.err,
    err,
    minIon = 0.02,
    polarity_ft = 0,
    polarity_qtof = 0,
    FT_path,
    QTOF_path,
    spectrum_type_FT = NULL,
    spectrum_type_QTOF = NULL
) {
  
  # Load FT compounds
  ft.exp <- DBI::dbGetQuery(
    FT_con,
    "
    SELECT
      retention_time,
      mass_measured,
      compound_id,
      expid
    FROM ms_compound
    WHERE expid = ?
    ",
    params = list(FT_expnr)
  )
  
  
  # Load QTOF compounds
  syn.exp <- DBI::dbGetQuery(
    QTOF_con,
    "
    SELECT
      retention_time,
      mass_measured,
      compound_id,
      expid
    FROM ms_compound
    WHERE expid = ?
    ",
    params = list(QTOF_expnr)
  )
  
  
  if (nrow(ft.exp) == 0) {
    warning(
      "No compounds found for FT experiment ",
      FT_expnr,
      "."
    )
    return(Assoc)
  }
  
  if (nrow(syn.exp) == 0) {
    warning(
      "No compounds found for QTOF experiment ",
      QTOF_expnr,
      "."
    )
    return(Assoc)
  }
  
  
  # Check spectrum_type columns
  ft_cols <- DBI::dbListFields(
    FT_con,
    "msms_spectrum"
  )
  
  qtof_cols <- DBI::dbListFields(
    QTOF_con,
    "msms_spectrum"
  )
  
  has_type_FT <- "spectrum_type" %in% ft_cols
  has_type_QTOF <- "spectrum_type" %in% qtof_cols
  
  
  # Check spectrum type arguments
  if (has_type_FT && is.null(spectrum_type_FT)) {
    stop(
      "`spectrum_type_FT` must be provided because ",
      "`spectrum_type` exists in the FT database."
    )
  }
  
  if (has_type_QTOF && is.null(spectrum_type_QTOF)) {
    stop(
      "`spectrum_type_QTOF` must be provided because ",
      "`spectrum_type` exists in the QTOF database."
    )
  }
  
  
  # FT MS2 spectra
  
  ft_query <- "
    SELECT
      s.compound_id,
      p.mz
    FROM msms_spectrum s
    JOIN msms_spectrum_peak p
      USING(spectrum_id)
    JOIN ms_compound c
      USING(compound_id)
    WHERE s.ms_level = ?
      AND s.polarity = ?
      AND c.expid = ?
  "
  
  ft_params <- list(
    2,
    polarity_ft,
    FT_expnr
  )
  
  
  # Filter by spectrum_type only if column exists
  if (has_type_FT) {
    
    ft_query <- paste0(
      ft_query,
      " AND s.spectrum_type = ?"
    )
    
    ft_params <- c(
      ft_params,
      list(spectrum_type_FT)
    )
  }
  
  
  ft_query <- paste0(
    ft_query,
    " ORDER BY s.compound_id, p.mz"
  )
  
  
  ms2_ft <- DBI::dbGetQuery(
    FT_con,
    ft_query,
    params = ft_params
  )
  
  
  # Split FT peaks by compound
  ms2_ft_list <- split(
    round(ms2_ft$mz),
    ms2_ft$compound_id
  )
  
  ms2_ft_list <- lapply(
    ms2_ft_list,
    unique
  )
  
  
  # QTOF MS2 spectra
  
  qtof_query <- "
    SELECT
      s.compound_id,
      p.mz
    FROM msms_spectrum s
    JOIN msms_spectrum_peak p
      USING(spectrum_id)
    JOIN ms_compound c
      USING(compound_id)
    WHERE s.ms_level = ?
      AND s.polarity = ?
      AND c.expid = ?
  "
  
  qtof_params <- list(
    2,
    polarity_qtof,
    QTOF_expnr
  )
  
  
  # Filter by spectrum_type only if column exists
  if (has_type_QTOF) {
    
    qtof_query <- paste0(
      qtof_query,
      " AND s.spectrum_type = ?"
    )
    
    qtof_params <- c(
      qtof_params,
      list(spectrum_type_QTOF)
    )
  }
  
  
  qtof_query <- paste0(
    qtof_query,
    " ORDER BY s.compound_id, p.mz"
  )
  
  
  ms2_qtof <- DBI::dbGetQuery(
    QTOF_con,
    qtof_query,
    params = qtof_params
  )
  
  
  # Split QTOF peaks by compound
  ms2_qtof_list <- split(
    round(ms2_qtof$mz),
    ms2_qtof$compound_id
  )
  
  ms2_qtof_list <- lapply(
    ms2_qtof_list,
    unique
  )
  
  
  # Loop over FT compounds
  
  for (i in seq_len(nrow(ft.exp))) {
    
    COMPID <- ft.exp$compound_id[i]
    x.tR <- ft.exp$retention_time[i]
    
    
    if (is.na(x.tR) || x.tR < cutoff) {
      next
    }
    
    
    # Predict QTOF retention time
    t1.tR <- ifelse(
      rg[3] != 0 && x.tR > rg[3],
      x.tR - rg[3],
      0
    )
    
    t2.tR <- ifelse(
      rg[5] != 0 && x.tR > rg[5],
      x.tR - rg[5],
      0
    )
    
    
    y.tR <- rg[1] +
      rg[2] * x.tR +
      rg[4] * t1.tR +
      rg[6] * t2.tR
    
    
    y.l <- y.tR - lc.err
    y.h <- y.tR + lc.err
    
    
    # Candidate matches based on precursor mass
    pres <- Find_cand_matches_SQL(
      COMPID,
      err,
      syn.exp,
      ft.exp
    )
    
    
    if (is.null(pres) || nrow(pres) == 0) {
      next
    }
    
    
    # FT MS2 spectrum
    ms2ion <- ms2_ft_list[[as.character(COMPID)]]
    
    
    if (is.null(ms2ion) ||
        length(ms2ion) == 0) {
      next
    }
    
    
    pres1 <- data.frame()
    w.diff <- numeric()
    
    
    # Test candidate QTOF compounds
    for (w in seq_len(nrow(pres))) {
      
      rt_w <- pres$retention_time[w]
      target_compid <- pres$compound_id[w]
      
      
      # Retention-time window
      if (is.na(rt_w) ||
          rt_w <= y.l ||
          rt_w >= y.h) {
        next
      }
      
      
      # QTOF MS2 spectrum
      msmsion <- ms2_qtof_list[[as.character(target_compid)]]
      
      
      if (is.null(msmsion) ||
          length(msmsion) == 0) {
        next
      }
      
      
      # Number of common MS2 fragment ions
      same_ion <- length(
        intersect(
          ms2ion,
          msmsion
        )
      )
      
      
      # Fraction of FT ions found in QTOF
      matching_ratio <- same_ion /
        length(ms2ion)
      
      
      if (matching_ratio >= minIon) {
        
        pres1 <- rbind(
          pres1,
          pres[w, , drop = FALSE]
        )
        
        w.diff <- c(
          w.diff,
          abs(y.tR - rt_w)
        )
      }
    }
    
    
    # No valid candidate
    if (nrow(pres1) == 0) {
      next
    }
    
    
    # Select candidate closest to predicted RT
    pres1 <- pres1[
      order(w.diff),
      ,
      drop = FALSE
    ]
    
    
    # Create association
    new_row <- data.frame(
      ref_compid = COMPID,
      target_compid = pres1$compound_id[1],
      ref_database = basename(FT_path),
      target_database = basename(QTOF_path),
      stringsAsFactors = FALSE
    )
    
    
    # Check whether association already exists
    is_new <- !any(
      Assoc$ref_compid == new_row$ref_compid &
        Assoc$target_compid == new_row$target_compid &
        Assoc$ref_database == new_row$ref_database &
        Assoc$target_database == new_row$target_database
    )
    
    
    if (is_new) {
      Assoc <- rbind(
        Assoc,
        new_row
      )
    }
  }
  
  
  return(Assoc)
}