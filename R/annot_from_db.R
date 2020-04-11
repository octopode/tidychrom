#' Annotate from Database
#'
#' Simple iterative database search to annotate spectra in a tibble
#' @param peaks_matched Tibble with column 'spectrum' containing MSnBase Spectrum objects encapsulated in lists
#' @param dir_mgf Directory containing MS1 spectra (.mgf format, suffix) annotated using TITLE= tag
#' @param fname_regex Regex string that can be used to prefilter the TITLE= tags
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
  #spectra_list <- spectra %>%
  #  group_split()

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

annot_from_db_old <- function(peaks_matched, dir_mgf, fname_regex){
  #dir_mgf = "/Users/jwinnikoff/Documents/MBARI/Lipids/GCMSData/db/methyl "
  if(missing(fname_regex)){
    mgfs <- list.files(path = dir_mgf, pattern = ".mgf", full.names = T)
  }else{
    # use grep; won't work on Windoze
    mgfs <- system2(c("grep", "-lR", fname_regex, dir_mgf), stdout = T)
  }

  # previously, peaks_matched_all was filtered:
  #peaks_matched_stds <- peaks_matched_all %>%
  #  filter(samp == "Supel37_1_8") %>%
  peaks_matched <- peaks_matched %>%
    mutate(
      cos_db = 0,
      id_db = "",
      spectrum_db = c(),
      file_db = ""
    )
  # for each spectrum in the library
  for(i in seq(length(mgfs))){
    #spec_db <- readMgfData("/Users/jwinnikoff/Documents/MBARI/Lipids/GCMSData/MoNA/mgf/JP000010.mgf")
    spec_db <- readMgfData(mgfs[i])
    filename <- spec_db@phenoData$sampleNames %>%
      as.character() %>%
      basename()
    title <- spec_db@featureData@data$TITLE %>%
      as.character()
    #print(plot_spectrum(spec_db@assayData$X1))
    message(paste("(", i, "/", length(mgfs), ") checking ", filename, sep = "", collapse = ""))
    peaks_matched <- peaks_matched %>%
      mutate(
        # calculate cosine distance for each standard spectrum
        cos_new = compareSpectra(spec_db@assayData$X1, spectrum, fun="dotproduct"),
        # if greater than current cosine distance
        # replace hit identity
        id_db = ifelse(cos_new > cos_db, title, id_db),
        spectrum_db = ifelse(cos_new > cos_db, c(spec_db@assayData$X1), c(spectrum_db)),
        file_db = ifelse(cos_new > cos_db, filename, file_db),
        # and distance value
        cos_db = ifelse(cos_new > cos_db, cos_new, cos_db)
      )
  }
  peaks_matched <- peaks_matched %>%
    select(-cos_new)
  return(peaks_matched)
}
