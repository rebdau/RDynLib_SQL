#' @title Local alignment between FTMS and QTOF datasets
#'
#' @description
#'
#' The `Aligning_General_SQL()` function performs local alignment between two
#' LC-MS datasets stored in SQL databases (FTMS and QTOF).  
#' The function:
#' - extracts compound data from the both SQL databases,
#' - performs local alignment using `LCalign()`,
#' - adjusts `mass_measured` values when the two experiments have different 
#'   polarities.
#' @param FT_con `character(1)` A DBI database connection object to the FTMS 
#'        SQL database.
#'
#' @param QTOF_con `character(1)` A DBI database connection object to the QTOF
#'        SQL database.
#'
#' @param FT_expnr `numeric(1)` Experiment number (expid) to load from the FTMS
#'         SQL database.
#'
#' @param QTOF_expnr `numeric(1)` Experiment number to load from the QTOF SQL
#'         database.
#'
#' @param err `numeric(1)` is the m/z tolerance window used to decide which QTOF
#'        peaks are potential isomeric matches for a given FT peak.
#'
#' @param t.ini `numeric(1)`The number of neighboring peaks
#'
#' @return A data.frame representing the locally aligned features
#' 
#' @importFrom DBI dbGetQuery
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
Aligning_General_SQL <- function(FT_con, QTOF_con, FT_expnr, QTOF_expnr, err,
                                 t.ini) {

  ft_sql <- dbGetQuery(
    FT_con,
    paste0(
      "SELECT c.*
     FROM ms_compound c
     WHERE c.expid = ", FT_expnr, "
       AND EXISTS (
         SELECT 1
         FROM msms_spectrum p
         WHERE p.compound_id = c.compound_id
       )"
    )
  )
  
              
    
  syn_sql <- dbGetQuery(
    QTOF_con,
    paste0(
      "SELECT c.*
     FROM ms_compound c
     WHERE c.expid = ", QTOF_expnr, "
       AND EXISTS (
         SELECT 1
         FROM msms_spectrum p
         WHERE p.compound_id = c.compound_id
       )"
    )
  )
  

    ft_mode  <- dbGetQuery(FT_con, "SELECT mode FROM experiment")$mode[1]
    syn_mode <- dbGetQuery(QTOF_con, "SELECT mode FROM experiment")$mode[1]

    ft.exp <- ft_sql[, c("retention_time", "mass_measured", "compound_id",
                         "expid", "name", "smiles", "nodename")]
    syn.exp <- syn_sql[, c("retention_time", "mass_measured", "compound_id",
                           "expid", "name", "smiles", "nodename")]
    ft.exp.o <- ft.exp[order(ft.exp$retention_time), ]

    reg <- nrow(ft.exp.o)
    if (reg %% 2 != 0) reg <- reg - 1
    COMPID <- ft.exp.o[reg/2, "compound_id"]

    ft.sh <- Make_loc_shortlist(COMPID, reg, ft.exp)

    if (ft_mode != syn_mode) {
        syn.exp$mass_measured <- syn.exp$mass_measured - 2.01456
    }

    ft.sh <- Remov_empties(err, syn.exp, ft.sh, reg)

    LCal <- LCalign(err, t.ini, syn.exp, ft.sh)

    if (ft_mode != syn_mode) {
        LCal[,5] <- as.numeric(LCal[,5]) + 2.01456
    }

    return(LCal)
}
