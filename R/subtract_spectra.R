#' Subtract Spectra
#'
#' Take two MSnbase Spectrum objects and subtract one from the other, no scaling or anything.
#' Return another Spectrum. Use for subtracting a blank baseline.
#' @param spec_samp The spectrum with material in it
#' @param spec_blank The blank spectrum
#' @param bin Bin spectra to this m/z tolerance. If missing, do not bin at all.
#' @keywords subtract spectrum spectra baseline blank
#' @export
#' @examples
#' spec_corrected <- subtract_spectrum(specSamp, specBL)
#'
subtract_spectra <- function(spec_samp, spec_bl, bin){
  # bin both spectra if desired
  if(!missing(bin)){
    spec_samp <- spec_samp %>%
      bin(breaks = seq(0.5, max(.@mz), bin))
    spec_samp <- spec_bl %>%
      bin(breaks = seq(0.5, max(.@mz), bin))
  }
  spec_out <- spec_samp %>%
    spec2peaklist() %>%
    # keep all peaks present in sample spectrum
    left_join(
      spec_bl %>%
        spec2peaklist(),
      by = "mz"
      ) %>%
    replace_na(list(intensity.y = 0)) %>%
    mutate(intensity = intensity.x - intensity.y) %>%
    filter(intensity >= 0) %>%
    select(mz, intensity) %>%
    # pass along RT of the sample spectrum
    peaklist2spec(., spec_samp@rt)

  return(spec_out)
}
