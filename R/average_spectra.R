#' Average Spectra
#'
#' Take a list of MSnbase Spectrum objects and average the intensities, no scaling or anything.
#' Return another Spectrum. Use e.g. for averaging blank runs.
#' @param specs A list of spectra
#' @param bin Bin spectra to this m/z tolerance. If missing, do not bin at all.
#' @keywords subtract spectrum spectra baseline blank
#' @export
#' @examples
#' spec_avg <- average_spectrum(spec1, spec2)
#'
average_spectra <- function(specs, bin){
  # bin spectra if desired
  if(!missing(bin)){
    spectra <- lapply(spectra, function(x){bin(breaks = seq(0.5, max(x@mz), bin))})
  }
  # iteratively merge peaklists
  peaklists <- specs[[1]] %>%
    spec2peaklist()
  for(i in seq(2, length(specs))){
    peaklists <- peaklists %>%
      # keep all peaks present in sample spectrum
      full_join(
        specs[[i]] %>%
          spec2peaklist(),
        by = "mz"
      )
  }
  # named vector %>% list required for replace_na() function
  repl_na_list <- rep(0, ncol(peaklists) - 1)
  names(repl_na_list) <- names(peaklists)[2:ncol(peaklists)]
  spec_out <- peaklists %>%
    # replace NAs with 0s
    replace_na(repl_na_list %>% as.list()) %>%
    # average instensities
    mutate(intensity = select(., -mz) %>% rowMeans(.)) %>%
    # return as spectrum
    select(mz, intensity) %>%
    peaklist2spec()

  return(spec_out)
}
