# take MSnExp/raw data, sample index, and rt
# return mass spectrum
get_spec_rt <- function(exp, index_samp, rt){
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
