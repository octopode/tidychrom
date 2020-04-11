# take Spectrum, return base m/z
get_basepeak <- function(spectrum){
  basepeak <- spec2peaklist(spectrum) %>%
    filter(intensity == max(intensity)) %>%
    select(mz) %>%
    unlist() %>% unname()
  return(basepeak)
}
