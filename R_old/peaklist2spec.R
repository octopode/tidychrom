#' Convert (tbl) Peaklist to Spectrum
#'
#' Take tibble with mz and intensity columns (rt column optional).
#' Return a Spectrum object.
#' @param peaklist The tibble in question
#' @keywords convert encode spectrum peaklist
#' @export
#' @examples
#' my_spectrum <- peaklist2spec(my_peaklist)
#'
peaklist2spec <- function(peaklist){
  #NTS: should probably clean up this conditional and use ifelse() for rt, too -20200402
  if("rt" %in% names(peaklist)){
    # uses first value in the RT column. Hopefully they're all the same!
    spectrum <- new(
      "Spectrum1",
      mz = peaklist$mz,
      intensity = peaklist$intensity,
      rt = peaklist$rt[[1]],
      fromFile = ifelse("samp_num" %in% names(peaklist), peaklist$samp_num[[1]] %>% as.integer(), integer(0)),
      polarity = ifelse("polarity" %in% names(peaklist), peaklist$polarity[[1]] %>% as.integer(), NA %>% as.integer())
      )
  }else{
    spectrum <- new(
      "Spectrum1",
      mz = peaklist$mz,
      intensity = peaklist$intensity)
  }

  return(spectrum)
}
