#' Annotate from Database
#'
#' Simple iterative database search to annotate spectra in a tibble
#' @param spectra Grouped (e.g. by \code{scan} tibble with columns \code{mz}, \code{intensity}.
#' @param dir_db Directory containing database spectrum files in an mzR-readable format (mzXML, mzML, etc.)
#' @param fname_regex Regex string that can be used to prefilter database files (generally the suffix.)
#' @keywords annot annotate db database
#' @export
#' @examples
#' # effectively add an annotation column to peaks_matched:
#' peaks_matched <- annot_from_db(peaks_matched, "~/msdb/", "methyl")

annot_from_db <- function(spectra, dir_db, bin, fname_regex = "*.mzXML", cores = 1){
  ## expect spectra to come in grouped
  if(!length(group_vars(spectra))){
    stop("query spectra must be grouped!")
  }

  if(cores > 1){
    clust <- makeCluster(cores, type="FORK")
    on.exit(stopCluster(clust)) # important!
  }else{
    clust <- NULL
  }

  #dir_db = "/Users/jwinnikoff/Documents/MBARI/Lipids/GCMSData/db/supel37"
  # system find won't work on Windoze!
  # there's probably another good nonrecursive solution
  files_db <- system2(c("find", dir_db, "-name", fname_regex, "-type", "f"), stdout = T)

  if(!missing(bin)){
    # if custom bin specified,
    # bin both the query and db spectra as requested
    bin_width_mz <- bin
    spectra <- spectra %>%
      bin_spectra(bin_width = bin)
  }else{
    # by default, bin db spectra to integer
    bin_width_mz <- 1
  }

  # again, as a disk I/O operation, aim to open and close each file only once
  ids <- pblapply(
    cl = clust,
    files_db,
    function(file_load){
      spectra %>%
        # generate a column of cosine scores for each file
        summarize(
          cos = cosine_spectra(
            # spectrum data from the present group
            cbind(mz, intensity) %>%
              as_tibble(),
            # spectrum data loaded from file
            file_load %>%
              read_tidymass() %>%
              bin_spectra(bin_width = bin_width_mz)
          )
        ) %>%
        select(cos) %>%
        # name the column for the file
        set_names(basename(file_load))
        # NTS: this would be better, but I can wait to make it work
        #summarize_at(
        #  vars(file_load) = cosine_spectra(
        #    # spectrum data from the present group
        #    cbind(mz, intensity) %>%
        #      as_tibble(),
        #    # spectrum data loaded from file
        #    file_load %>%
        #      read_tidymass() %>%
        #      bin_spectra(bin_width = bin_width_mz)
        #  )
        #)
    }
  ) %>%
    # bind all the columns together
    do.call(cbind, .) %>%
    as_tibble() %>%
    # label them with ROI
    # I sure hope these are in the right order!!
    bind_cols(group_keys(spectra)) %>%
    # melt it down
    gather(-roi, key = "file", value = "cos") %>%
    # and get the best match for each ROI
    group_by_at(group_vars(spectra)) %>%
    filter(cos == max(cos))

  return(ids)
}
