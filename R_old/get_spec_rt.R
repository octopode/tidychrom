# take MSnExp/raw data, sample index, and rt
# return mass spectrum
# old version as of 20200330
# fast, but not robust against different run parameters (scan timing, etc)
# NTS 20200330: build in a conditional check to see if method for all runs is identical
get_spec_rt <- function(exp, index_samp, rt, rtime_tbl){
  # number of samples
  num_samps <- length(fileNames(exp))
  # number of spectra per sample; runs should be identical
  num_spectra <- length(exp) / num_samps
  rt_beg <- min(rtime(exp))
  rt_end <- max(rtime(exp))
  # which spectrum does rt correspond to?
  index_spec <- ceiling(num_spectra * (rt - rt_beg) / (rt_end - rt_beg))
  # put leading zeroes on indices
  index_samp <- str_pad(index_samp, width=nchar(num_samps), side="left", pad="0")
  index_spec <- str_pad(index_spec, width=nchar(num_spectra), side="left", pad="0")
  feature <- paste("F", index_samp, ".S", index_spec, sep = "", collapse = "")
  return(exp[[feature]])
}

## more robust algorithm
## slower?

#
#get_spec_rt <- function(exp, index_samp, rt, rtime_tbl){
#  # rtime_tbl can be passed into the call to save time.
#  # With a large experiment, the gather() call takes 5ever.
#  if(missing(rtime_tbl)){
#    rtime_tbl <- exprtime2tbl(exp)
#  }
#
#  feature <- rtime_tbl %>%
#    # first, filter experiment for requested sample
#    filter(samp_num == index_samp) %>%
#    # get nearest scan to query RT
#    mutate(delta = abs(rtime - rt)) %>%
#    # need to group for delta == min(delta) to work
#    group_by(1) %>%
#    # get the best-hit feature
#    filter(delta == min(delta)) %>%
#    pull(feature)
#
#  return(exp[[feature]])
#}
#
