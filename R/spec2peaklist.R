# take Spectrum, return tbl peaklist
spec2peaklist <- function(spectrum, bin){
  # bin the mz column if desired
  if(!missing(bin)){
    spectrum <- spectrum %>%
      bin(breaks = seq(0.5, max(.@mz), bin))
  }
  peaklist <- cbind(spectrum@mz, spectrum@intensity) %>%
  # builtin "as" method more versatile?
  # looks like there's no difference
  #peaklist <- spectrum %>%
    #as("data.frame") %>%
    as.tibble() %>%
    setNames(c("mz", "intensity"))
  return(peaklist)
}
