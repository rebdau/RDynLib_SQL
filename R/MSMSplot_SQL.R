#' @title plot MS2 spectra 
#'
#' @import DBI
#' @import RSQLite
#' 
#' @author Ahlam Mentag
#' 
#' @export
MSMSplot_SQL <- function(sql_path,
                         dbkey,
                         spectrum_id,
                         prdion,
                         neutloss,
                         err = 0.015,
                         minum = 2,
                         oldpar = NULL) {
  
  con <- dbConnect(SQLite(), sql_path)
  on.exit(dbDisconnect(con), add = TRUE)
  
  # precursor
  compound <- dbGetQuery(con, sprintf(
    "SELECT precursor_mz
     FROM msms_spectrum
     WHERE spectrum_id = %d
     AND ms_level = 2",
    spectrum_id
  ))
  
  if (nrow(compound) == 0)
    stop("Compound not found")
  
  precursor_mz <- compound$precursor_mz
  
  # peaks
  peaks <- dbGetQuery(con, sprintf(
    "SELECT p.mz, p.intensity
     FROM msms_spectrum s
     JOIN msms_spectrum_peak p ON s.spectrum_id = p.spectrum_id
     WHERE s.spectrum_id = %d
     ORDER BY p.mz",
    spectrum_id
  ))
  
  if (nrow(peaks) == 0)
    stop("No MS/MS peaks found")
  
  prod_ion <- round(peaks$mz, 2)
  intensity <- peaks$intensity
  
  sel <- which(intensity >= minum)
  
  # =========================
  # PRODUCT ION ANNOTATION
  # =========================
  cat("\nCandidate product ions:\n")
  prod_ann <- ProdIonMatch(prod_ion, prdion, err = err)
  
  # =========================
  # NEUTRAL LOSS (AUTO)
  # =========================
  neut_loss <- round(precursor_mz - prod_ion, 2)
  
  cat("\nCandidate neutral losses:\n")
  neut_ann <- NeutLossMatch(neut_loss, neutloss, err = err)
  
  # =========================
  # COMPLEMENTARY IONS
  # =========================
  cat("\nComplementary product ions:\n")
  ComplemIons_SQL(prod_ion, intensity, neut_loss, err)
  
  # =========================
  # PLOT
  # =========================
  par(mfrow = c(1,1))
  par(cex = 0.7)
  
  plot(prod_ion,
       intensity,
       type = "h",
       xlab = "m/z",
       ylab = "intensity",
       main = round(precursor_mz, 2))
  
  # labels: product ions
  text(prod_ion[sel],
       intensity[sel],
       labels = prod_ion[sel],
       pos = 3)
  
  # labels: neutral losses (automatic overlay)
  text(prod_ion[sel],
       intensity[sel],
       labels = neut_loss[sel],
       pos = 1,
       col = "red")
  
  if (!is.null(oldpar)) par(oldpar)
  
  return(list(
    product_ions = prod_ion,
    intensities = intensity,
    neutral_losses = neut_loss
  ))
}