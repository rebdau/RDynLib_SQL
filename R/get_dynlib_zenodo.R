#' Retrieve a DynLib database from Zenodo
#'
#' Downloads a published DynLib database from Zenodo and stores the
#' original file in a BiocFileCache cache. A copy of the cached database
#' is then created in `destdir`, allowing the working database to be
#' modified without changing the cached original.
#'
#' @param doi `character(1)`. DOI of the Zenodo record, for example
#'   `"10.5281/zenodo.1234567"`.
#'
#' @param file `character(1)`. Name of the file to retrieve from the
#'   Zenodo record.
#'
#' @param destdir `character(1)`. Directory in which the working copy
#'   should be created.
#'
#' @param overwrite `logical(1)`. If `FALSE`, an existing working copy
#'   is kept. If `TRUE`, it is replaced by a fresh copy of the cached
#'   Zenodo file.
#'
#' @return The normalized local path to the working copy of the database.
#'
#' @importFrom BiocFileCache BiocFileCache bfcrpath
#' @importFrom jsonlite fromJSON
#'
#' @export
get_dynlib_zenodo <- function(doi,
                              file,
                              destdir = ".",
                              overwrite = FALSE) {
  
  ## Check inputs
  stopifnot(
    is.character(doi),
    length(doi) == 1L,
    !is.na(doi),
    nzchar(doi)
  )
  
  stopifnot(
    is.character(file),
    length(file) == 1L,
    !is.na(file),
    nzchar(file)
  )
  
  stopifnot(
    is.character(destdir),
    length(destdir) == 1L,
    !is.na(destdir),
    nzchar(destdir)
  )
  
  stopifnot(
    is.logical(overwrite),
    length(overwrite) == 1L,
    !is.na(overwrite)
  )
  
  
  ## Accept DOI URL as well as plain DOI
  doi <- sub(
    "^https?://(dx\\.)?doi\\.org/",
    "",
    doi
  )
  
  
  ## Extract Zenodo record ID
  record_id <- sub(
    "^10\\.5281/zenodo\\.",
    "",
    doi
  )
  
  if (
    identical(record_id, doi) ||
    !grepl("^[0-9]+$", record_id)
  ) {
    stop(
      "'doi' does not look like a Zenodo DOI: ",
      doi
    )
  }
  
  
  ## Retrieve Zenodo metadata
  api_url <- paste0(
    "https://zenodo.org/api/records/",
    record_id
  )
  
  record <- jsonlite::fromJSON(
    api_url,
    simplifyVector = FALSE
  )
  
  
  ## Get available filenames
  available_files <- vapply(
    record$files,
    function(x) x$key,
    character(1)
  )
  
  
  ## Check that requested file exists
  idx <- match(file, available_files)
  
  if (is.na(idx)) {
    stop(
      "File '", file,
      "' was not found in Zenodo record ",
      doi,
      ".\nAvailable files:\n",
      paste(
        paste0(" - ", available_files),
        collapse = "\n"
      )
    )
  }
  
  
  ## Extract the download URL
  zenodo_file <- record$files[[idx]]
  
  url <- zenodo_file$links$self
  
  if (
    is.null(url) ||
    !is.character(url) ||
    length(url) != 1L
  ) {
    stop(
      "Could not determine the download URL for file '",
      file,
      "'."
    )
  }
  
  
  ## Create/open BiocFileCache without interactive questions
  bfc <- BiocFileCache::BiocFileCache(
    ask = FALSE
  )
  
  
  ## Unique cache name
  cache_name <- paste(
    "RDynLib",
    doi,
    file,
    sep = "::"
  )
  
  
  ## Retrieve cached copy or download once
  cached_file <- BiocFileCache::bfcrpath(
    bfc,
    rnames = cache_name,
    fpath = url
  )
  
  
  ## Create destination directory if necessary
  if (!dir.exists(destdir)) {
    dir.create(
      destdir,
      recursive = TRUE
    )
  }
  
  
  ## Working copy
  destfile <- file.path(
    destdir,
    basename(file)
  )
  
  
  ## Keep existing working copy unless overwrite = TRUE
  if (file.exists(destfile) && !overwrite) {
    
    message(
      "Using existing working copy: ",
      destfile
    )
    
    return(
      normalizePath(
        destfile,
        mustWork = TRUE
      )
    )
  }
  
  
  ## Copy cached original to working directory
  copied <- file.copy(
    from = cached_file,
    to = destfile,
    overwrite = TRUE
  )
  
  if (!copied) {
    stop(
      "Unable to copy cached DynLib database to: ",
      destfile
    )
  }
  
  
  message(
    "Working copy created: ",
    destfile
  )
  
  
  normalizePath(
    destfile,
    mustWork = TRUE
  )
}