match_FTexp_SQL <- function(
    LCal,
    con,
    expnr1,
    expnr2,
    thr1 = 1,
    thr2 = 0.9,
    IntThres = 100,
    tol = 0.01,
    spectrum_type_ref = NULL,
    spectrum_type_target = NULL
) {
  
  # Check LCal
  if (is.null(LCal) || nrow(LCal) == 0) {
    return(LCal)
  }
  
  
  # Check spectrum_type column
  has_spectrum_type <- "spectrum_type" %in%
    DBI::dbListFields(
      con,
      "msms_spectrum"
    )
  
  
  # If spectrum_type exists, types must be provided
  if (has_spectrum_type) {
    
    if (is.null(spectrum_type_ref)) {
      stop(
        "`spectrum_type_ref` must be provided because ",
        "`spectrum_type` exists in the database."
      )
    }
    
    if (is.null(spectrum_type_target)) {
      stop(
        "`spectrum_type_target` must be provided because ",
        "`spectrum_type` exists in the database."
      )
    }
  }
  
  
  # --------------------------------------------------
  # Reference MS2 query
  # --------------------------------------------------
  
  query_ref <- "
    SELECT
      s.compound_id,
      p.mz,
      p.intensity
    FROM msms_spectrum s
    JOIN msms_spectrum_peak p
      USING(spectrum_id)
    JOIN ms_compound c
      USING(compound_id)
    WHERE s.ms_level = ?
      AND c.expid = ?
  "
  
  params_ref <- list(
    2,
    expnr1
  )
  
  
  if (has_spectrum_type) {
    
    query_ref <- paste0(
      query_ref,
      " AND s.spectrum_type = ?"
    )
    
    params_ref <- c(
      params_ref,
      list(spectrum_type_ref)
    )
  }
  
  
  # Extract reference MS2 spectra
  ms2_exp1 <- DBI::dbGetQuery(
    con,
    query_ref,
    params = params_ref
  )
  
  
  # --------------------------------------------------
  # Target MS2 query
  # --------------------------------------------------
  
  query_target <- "
    SELECT
      s.compound_id,
      p.mz,
      p.intensity
    FROM msms_spectrum s
    JOIN msms_spectrum_peak p
      USING(spectrum_id)
    JOIN ms_compound c
      USING(compound_id)
    WHERE s.ms_level = ?
      AND c.expid = ?
  "
  
  params_target <- list(
    2,
    expnr2
  )
  
  
  if (has_spectrum_type) {
    
    query_target <- paste0(
      query_target,
      " AND s.spectrum_type = ?"
    )
    
    params_target <- c(
      params_target,
      list(spectrum_type_target)
    )
  }
  
  
  # Extract target MS2 spectra
  ms2_exp2 <- DBI::dbGetQuery(
    con,
    query_target,
    params = params_target
  )
  
  
  # Check reference spectra
  if (nrow(ms2_exp1) == 0) {
    
    if (has_spectrum_type) {
      
      warning(
        "No MS2 spectra found for experiment ",
        expnr1,
        " with spectrum_type = '",
        spectrum_type_ref,
        "'."
      )
      
    } else {
      
      warning(
        "No MS2 spectra found for experiment ",
        expnr1,
        "."
      )
    }
    
    return(
      LCal[0, , drop = FALSE]
    )
  }
  
  
  # Check target spectra
  if (nrow(ms2_exp2) == 0) {
    
    if (has_spectrum_type) {
      
      warning(
        "No MS2 spectra found for experiment ",
        expnr2,
        " with spectrum_type = '",
        spectrum_type_target,
        "'."
      )
      
    } else {
      
      warning(
        "No MS2 spectra found for experiment ",
        expnr2,
        "."
      )
    }
    
    return(
      LCal[0, , drop = FALSE]
    )
  }
  
  
  # Intensity filtering
  ms2_exp1 <- ms2_exp1[
    !is.na(ms2_exp1$intensity) &
      ms2_exp1$intensity >= IntThres,
    ,
    drop = FALSE
  ]
  
  ms2_exp2 <- ms2_exp2[
    !is.na(ms2_exp2$intensity) &
      ms2_exp2$intensity >= IntThres,
    ,
    drop = FALSE
  ]
  
  
  # Check peaks after intensity filtering
  if (nrow(ms2_exp1) == 0 ||
      nrow(ms2_exp2) == 0) {
    
    warning(
      "No MS2 peaks remained after intensity filtering."
    )
    
    return(
      LCal[0, , drop = FALSE]
    )
  }
  
  
  # Split spectra by compound
  sp1 <- split(
    ms2_exp1[
      ,
      c("mz", "intensity"),
      drop = FALSE
    ],
    ms2_exp1$compound_id
  )
  
  sp2 <- split(
    ms2_exp2[
      ,
      c("mz", "intensity"),
      drop = FALSE
    ],
    ms2_exp2$compound_id
  )
  
  
  # Store accepted LCal rows
  keep <- logical(
    nrow(LCal)
  )
  
  
  # Compare candidate pairs
  for (i in seq_len(nrow(LCal))) {
    
    comp1 <- as.character(
      LCal[i, 1]
    )
    
    comp2 <- as.character(
      LCal[i, 7]
    )
    
    
    spec1 <- sp1[[comp1]]
    spec2 <- sp2[[comp2]]
    
    
    # Both compounds must have MS2 spectra
    if (is.null(spec1) ||
        is.null(spec2)) {
      next
    }
    
    
    # Remove missing values
    spec1 <- spec1[
      !is.na(spec1$mz) &
        !is.na(spec1$intensity),
      ,
      drop = FALSE
    ]
    
    spec2 <- spec2[
      !is.na(spec2$mz) &
        !is.na(spec2$intensity),
      ,
      drop = FALSE
    ]
    
    
    # Both spectra must contain at least 2 ions
    if (nrow(spec1) < 2 ||
        nrow(spec2) < 2) {
      next
    }
    
    
    # Find matching ions
    matches <- vapply(
      seq_len(nrow(spec1)),
      function(j) {
        
        idx <- which(
          abs(
            spec2$mz -
              spec1$mz[j]
          ) <= tol
        )
        
        
        # No matching ion
        if (length(idx) == 0) {
          return(
            NA_integer_
          )
        }
        
        
        # If several ions match, select closest m/z
        idx[
          which.min(
            abs(
              spec2$mz[idx] -
                spec1$mz[j]
            )
          )
        ]
      },
      integer(1)
    )
    
    
    # Fraction of spectrum 1 found in spectrum 2
    common <- !is.na(
      matches
    )
    
    ratio <- sum(common) /
      nrow(spec1)
    
    
    # Minimum fraction of common ions
    if (ratio < thr1) {
      next
    }
    
    
    # Intensities of matching ions
    int1 <- spec1$intensity[
      common
    ]
    
    int2 <- spec2$intensity[
      matches[common]
    ]
    
    
    # At least two common ions required
    if (length(int1) < 2) {
      next
    }
    
    
    # Calculate vector norms
    norm1 <- sqrt(
      sum(int1^2)
    )
    
    norm2 <- sqrt(
      sum(int2^2)
    )
    
    
    # Avoid division by zero
    if (!is.finite(norm1) ||
        !is.finite(norm2) ||
        norm1 == 0 ||
        norm2 == 0) {
      next
    }
    
    
    # Normalize intensities
    int1 <- int1 / norm1
    int2 <- int2 / norm2
    
    
    # Dot product
    dot <- sum(
      int1 * int2
    )
    
    
    # Keep candidate
    if (is.finite(dot) &&
        dot >= thr2) {
      
      keep[i] <- TRUE
    }
  }
  
  
  # Return unique accepted matches
  unique(
    LCal[
      keep,
      ,
      drop = FALSE
    ]
  )
}