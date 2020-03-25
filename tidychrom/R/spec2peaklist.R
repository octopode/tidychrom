# take Spectrum, return tbl peaklist
spec2peaklist <- function(spectrum){
  peaklist <- cbind(spectrum@mz, spectrum@intensity) %>%
  # builtin "as" method more versatile?
  # looks like there's no difference
  #peaklist <- spectrum %>%
    #as("data.frame") %>%
    as.tibble() %>%
    setNames(c("mz", "intensity"))
  return(peaklist)
}
