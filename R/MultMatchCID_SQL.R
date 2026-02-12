MultMatchCID_SQL <- function(sql_path, dbkey, minIons = 4) {
  library(DBI)
  library(RSQLite)
  
  # Connect to SQLite database
  con <- dbConnect(SQLite(), sql_path)
  on.exit(dbDisconnect(con), add = TRUE)
  
  # Get target compound info
  compound <- dbGetQuery(con, sprintf(
    "SELECT compound_id, name, retention_time, expid
     FROM ms_compound
     WHERE compound_id = %d",
    dbkey
  ))
  
  if (nrow(compound) == 0) stop("Compound not found in ms_compound")
  
  # Get product ions of target compound
  unk <- dbGetQuery(con, sprintf(
    "SELECT DISTINCT p.Mz
     FROM msms_spectrum s
     JOIN msms_spectrum_peak p ON s.spectrum_id = p.spectrum_id
     WHERE s.compound_id = %d",
    dbkey
  ))$Mz
  
  if (length(unk) < minIons) stop("Not enough product ions for this compound")
  
  # Get all other compounds and their product ions
  all_compounds <- dbGetQuery(con, sprintf(
    "SELECT c.compound_id, c.name, c.retention_time, c.expid, p.Mz
     FROM ms_compound c
     JOIN msms_spectrum s ON c.compound_id = s.compound_id
     JOIN msms_spectrum_peak p ON s.spectrum_id = p.spectrum_id
     WHERE c.compound_id != %d",
    dbkey
  ))
  
 
  MS2list <- split(all_compounds$Mz, all_compounds$compound_id)
  

  msnSect <- MultMatch(MS2list, unk, minIons)
  

  ENTRIES <- subMultMatch(msnSect, unk, MS2list, minIons)
  

  fin <- data.frame(
    COMPID = integer(),
    mz = numeric(),
    FragIonCount = integer(),
    ComIons = integer(),
    DotIons = numeric(),
    ComLoss = numeric(),
    DotLoss = numeric(),
    COMPNAME = character(),
    tR = numeric(),
    EXPID = integer(),
    stringsAsFactors = FALSE
  )
  
  for (i in seq_along(ENTRIES)) {
    compid <- ENTRIES[i]
    out <- targMS2comp_SQL(dbkey1 = dbkey, dbkey2 = compid, MS2list = MS2list)
    
    comp_info <- dbGetQuery(con, sprintf(
      "SELECT name, retention_time, expid FROM ms_compound WHERE compound_id = %d",
      compid
    ))
    
    fin <- rbind(fin, data.frame(
      COMPID = compid,
      mz = NA,  
      FragIonCount = nrow(out),
      ComIons = out$ComIons,
      DotIons = out$DotIons,
      ComLoss = out$ComLoss,
      DotLoss = out$DotLoss,
      COMPNAME = comp_info$name,
      tR = comp_info$retention_time,
      EXPID = comp_info$expid,
      stringsAsFactors = FALSE
    ))
  }
  
  # Rank by combined score
  fin <- fin[order(fin$ComIons * fin$DotIons / fin$FragIonCount, decreasing = TRUE), ]
  
  writeLines(paste("The CID spectrum of compound", dbkey,
                   "matches to:"))
  
  return(fin)
}
