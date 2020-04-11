# take an MSnExp and a tibble of peaks
# return a list of XChromatogram objects bounded around the peaks
# note that the peaks tbl needs to have an "rt" column
divide_rois <- function(exp, peaks, lead = 1, lag = 1){
  chroms <- peaks %>%
    # add a margin to either side of the standard peak
    mutate(rtmin = rt - lead, rtmax = rt + lag) %>%
    select(rtmin, rtmax) %>%
    unname() %>%
    split(., seq(nrow(.))) %>%
    unname() %>%
    # pass the list of rtmin/max pairs
    lapply(., unlist) %>%
    # and make a chromatogram set for each
    lapply(., function(x){chromatogram(exp, rt = unlist(x), aggregationFun = "max")})
  return(chroms)
}
