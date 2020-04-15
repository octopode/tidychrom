#' Convert OnDiskMSnExp to MSNExp
#'
#' Take an OnDiskMSNExp and load it into memory.
#' Specifically, this means taking the files in it and loading a new experiment object.
#' @param exp OnDiskMSNExp object. Be sure to subset as desired beforehand (e.g. filterFile())
#' or else returned object could take gigs of RAM.
#' @keywords experiment disk memory
#' @export
#' @examples
#' fname_master_blank <- average_runs(exp, samp_indices = c(1,2,3))
#'
exp_disk2mem <- function(exp_disk){
  # "works," but seems to never finish :(
  # per https://github.com/lgatto/MSnbase/issues/258, backend = "inMemory" only actually works for MS2 data
  exp_mem <- readMSData(
    files = exp_disk %>% fileNames(),
    pdata = exp_disk %>% .@phenoData,
    msLevel = msLevel(exp_disk) %>% unique()
    )

  return(exp_mem)
}
