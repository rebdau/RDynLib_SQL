#' @title Plot MSn spectra
#'
#' @description
#' Plots individual MSn spectra associated with a compound or a specific
#' individual spectrum identified by its spectrum_id in a DynLib SQLite
#' database. When multiple spectra are plotted, they are ordered by
#' increasing MS level and decreasing precursor intensity within each
#' MS level.
#'
#' Spectra belonging to the same original experimental fragmentation
#' tree are identified by a common MSntreeID, whereas spectra belonging
#' to the same representative or extended fragmentation tree are grouped
#' by a common MergedMSntreeID.
#'
#' For spectra whose original fragmentation pathway can be reconstructed,
#' the complete sequence of precursor m/z values is displayed in the
#' plot title. For complementary spectra incorporated during tree
#' enrichment whose original fragmentation pathway cannot be
#' reconstructed, only the MS level and immediate precursor m/z
#' are displayed.
#'
#' MS2 precursor m/z values are displayed with two decimal places,
#' while higher-level precursor m/z values are displayed at unit precision.
#'
#' @param sql_path Character. Path to the DynLib SQLite database.
#' @param compound_id Integer or character. Compound identifier in the
#'   DynLib database. Optional if spectrum_id is provided.
#' @param prcx Numeric. Character expansion factor for fragment ion
#'   labels. Defaults to 0.6.
#' @param thr_int Numeric. Minimum relative intensity, expressed as a
#'   percentage of the base peak intensity, for plotting fragment ions.
#'   Defaults to 0.
#' @param label_peaks Logical. Whether to display fragment ion m/z
#'   labels. Defaults to TRUE.
#' @param thr_label Numeric. Minimum relative intensity, expressed as
#'   a percentage of the base peak intensity, for labelling fragment
#'   ions. Defaults to 5.
#' @param digits Integer. Number of decimal places for fragment ion
#'   m/z labels. Defaults to 0 (unit precision).
#' @param max_cols Integer. Maximum number of columns in the graphical
#'   layout. Defaults to 2.
#' @param spectrum_id Integer or character. Optional identifier of a
#'   specific individual spectrum to plot. If NULL, all individual
#'   spectra associated with the specified compound are plotted.
#'   Defaults to NULL.
#'
#' @details
#' Relative intensity thresholds are calculated independently for each
#' spectrum using its base peak intensity before filtering.
#' Absolute ion intensities are retained on the y-axis.
#'
#' Only individual spectra are plotted, including those originating
#' from the representative MSn tree and complementary spectra
#' incorporated during tree enrichment using \code{enrich_msntree()}.
#' Pseudo-MS2 spectra are excluded.
#'
#' The MSntreeID preserves the association of each spectrum with its
#' original experimental fragmentation tree, while the MergedMSntreeID
#' groups spectra belonging to the same representative or extended tree.
#'
#' Fragmentation paths are reconstructed using MSntreeID, acquisition
#' scan numbers, precursor scan numbers and data origin. A common
#' MergedMSntreeID does not necessarily imply a common experimental
#' fragmentation pathway. Therefore, fragmentation paths are not
#' inferred solely from membership in the same extended tree.
#'
#' For example, an MS4 spectrum with a reconstructed experimental
#' fragmentation pathway may have the title "MS4: 625.14 -> 463 -> 301",
#' whereas a complementary MS4 spectrum for which the original
#' fragmentation pathway cannot be reconstructed is displayed with
#' the title "MS4: 301".
#' 
#' Individual spectra associated with the compound but lacking tree
#' identifiers are also plotted, with titles indicating only their
#' MS level and immediate precursor m/z.
#' 
#' When a specific spectrum is requested, all individual spectra
#' associated with its compound are considered for fragmentation
#' path reconstruction, but only the requested spectrum is plotted.
#' If both compound_id and spectrum_id are provided, the spectrum
#' must belong to the specified compound.
#' 
#' @return
#' No return value. The function generates MSn spectrum plots on the
#' active graphics device.
#'
#' @importFrom Spectra peaksData
#' @importFrom graphics par plot plot.new title text
#' @importFrom DBI dbConnect dbDisconnect dbGetQuery
#' @importFrom RSQLite SQLite
#'
#' @author Ahlam Mentag, Rebecca Dauwe
#'
#' @export
plotMSn <- function(sql_path,
                    compound_id = NULL,
                    prcx = 0.6,
                    thr_int = 0,
                    label_peaks = TRUE,
                    thr_label = 5,
                    digits = 0,
                    max_cols = 2,
                    spectrum_id = NULL) {
  
  # Check input parameters
  if (!file.exists(sql_path)) {
    stop("The specified database does not exist.")
  }
  
  if (is.null(compound_id) && is.null(spectrum_id)) {
    stop("Provide either compound_id or spectrum_id.")
  }
  
  stopifnot(
    is.null(compound_id) || length(compound_id) == 1,
    is.null(spectrum_id) || length(spectrum_id) == 1,
    thr_int >= 0, thr_int <= 100,
    thr_label >= 0, thr_label <= 100,
    digits >= 0, digits == as.integer(digits),
    max_cols >= 1, max_cols == as.integer(max_cols)
  )
  
  # Retrieve compound_id if only spectrum_id is provided
  if (is.null(compound_id)) {
    
    con <- DBI::dbConnect(RSQLite::SQLite(), sql_path)
    on.exit(DBI::dbDisconnect(con), add = TRUE)
    
    result <- DBI::dbGetQuery(
      con,
      "SELECT DISTINCT compound_id
     FROM msms_spectrum
     WHERE spectrum_id = ?",
      params = list(spectrum_id)
    )
    
    if (nrow(result) != 1) {
      stop("No unique compound found for this spectrum_id.")
    }
    
    compound_id <- result$compound_id[1]
  }
  
  # Retrieve all spectrum metadata for the compound
  sp_id <- get_spectrum_id(
    sql_path = sql_path,
    compound_id = compound_id
  )
  
  # Select individual spectra
  sp_selected_id <- subset(
    sp_id,
    spectrum_type == "individual"
  )
  
  if (nrow(sp_selected_id) == 0) {
    stop("No individual spectra found for this compound.")
  }
  
  # Sort by MS level and decreasing precursor intensity
  sp_selected_id <- sp_selected_id[
    order(
      sp_selected_id$ms_level,
      -sp_selected_id$precursorIntensity
    ),
  ]
  
  totalplots <- nrow(sp_selected_id)
  
  # Generate spectrum titles
  plot_titles <- character(totalplots)
  
  for (i in seq_len(totalplots)) {
    
    current <- i
    level <- sp_selected_id$ms_level[i]
    
    # Format immediate precursor m/z
    path <- formatC(
      sp_selected_id$precursor_mz[i],
      format = "f",
      digits = if (level == 2) 2 else 0
    )
    
    # Reconstruct fragmentation path for connected spectra
    if (!is.na(sp_selected_id$MSntreeID[i])) {
      
      while (sp_selected_id$ms_level[current] > 2) {
        
        # Identify parent spectrum
        parent <- which(
          sp_selected_id$MSntreeID ==
            sp_selected_id$MSntreeID[current] &
            sp_selected_id$dataOrigin ==
            sp_selected_id$dataOrigin[current] &
            sp_selected_id$acquisitionNum ==
            sp_selected_id$precScanNum[current] &
            sp_selected_id$ms_level ==
            sp_selected_id$ms_level[current] - 1
        )
        
        # Retain only the immediate precursor if the
        # complete fragmentation path cannot be reconstructed
        if (length(parent) != 1) {
          
          path <- formatC(
            sp_selected_id$precursor_mz[i],
            format = "f",
            digits = 0
          )
          
          break
        }
        
        # Format parent precursor m/z
        parent_mz <- formatC(
          sp_selected_id$precursor_mz[parent],
          format = "f",
          digits = if (sp_selected_id$ms_level[parent] == 2) 2 else 0
        )
        
        path <- c(parent_mz, path)
        current <- parent
      }
    }
    
    # Save spectrum title
    plot_titles[i] <- paste0(
      "MS", level, ": ",
      paste(path, collapse = " \u2192 ")
    )
  }
  
  sp_selected_id$plot_title <- plot_titles
  
  # Optionally select a specific spectrum
  if (!is.null(spectrum_id)) {
    
    sp_selected_id <- sp_selected_id[
      which(sp_selected_id$spectrum_id == spectrum_id),
      , drop = FALSE
    ]
    
    if (nrow(sp_selected_id) == 0) {
      stop("The specified individual spectrum was not found.")
    }
  }
  
  totalplots <- nrow(sp_selected_id)
  
  # Retrieve spectra
  sp <- sql_to_spectra(sql_path)
  
  sp_selected <- sp[
    sp$spectrum_id %in% sp_selected_id$spectrum_id
  ]
  
  # Reorder spectra according to selected metadata
  idx <- match(
    sp_selected_id$spectrum_id,
    sp_selected$spectrum_id
  )
  
  stopifnot(!anyNA(idx))
  
  sp_selected <- sp_selected[idx]
  
  stopifnot(
    all(sp_selected$spectrum_id == sp_selected_id$spectrum_id)
  )
  
  # Retrieve peak data
  peaks <- peaksData(sp_selected)
  
  # Define graphical layout
  col.nr <- min(max_cols, totalplots)
  row.nr <- ceiling(totalplots / col.nr)
  
  # Save and restore graphical parameters
  oldpar <- par(no.readonly = TRUE)
  on.exit(par(oldpar), add = TRUE)
  
  par(
    mfrow = c(row.nr, col.nr),
    mar = c(3, 4, 2, 1)
  )
  
  # Plot spectra
  for (i in seq_along(peaks)) {
    
    p <- peaks[[i]]
    
    # Retrieve spectrum title
    plot_title <- sp_selected_id$plot_title[i]
    
    # Handle empty spectra
    if (nrow(p) == 0) {
      plot.new()
      title(plot_title)
      text(0.5, 0.5, "No peaks")
      next
    }
    
    # Calculate base peak intensity before filtering
    base_peak <- max(p[, "intensity"])
    
    # Filter peaks by relative intensity
    p <- p[
      p[, "intensity"] >= base_peak * thr_int / 100,
      , drop = FALSE
    ]
    
    # Handle spectra with no remaining peaks
    if (nrow(p) == 0) {
      plot.new()
      title(plot_title)
      text(0.5, 0.5, "No peaks above threshold")
      next
    }
    
    prod_ion <- p[, "mz"]
    intens <- p[, "intensity"]
    
    # Plot spectrum
    plot(
      prod_ion,
      intens,
      type = "h",
      xlab = "m/z",
      ylab = "Ion intensity",
      main = plot_title,
      ylim = c(0, if (base_peak > 0) base_peak * 1.15 else 1)
    )
    
    # Optionally label peaks based on relative intensity
    if (label_peaks && base_peak > 0) {
      
      rel_int <- 100 * intens / base_peak
      sel <- rel_int >= thr_label
      
      text(
        prod_ion[sel],
        intens[sel],
        labels = round(prod_ion[sel], digits),
        pos = 3,
        cex = prcx
      )
    }
  }
  
  invisible(NULL)
}
