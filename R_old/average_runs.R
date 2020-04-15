#' Average Chromatographic Runs
#'
#' Take an MSNExp and a vector of sample indices. Average those runs,
#' save the resultant CDF, and return that filename.
#' @param exp MSNExp object; all samples to be averaged need to be in the experiment.
#' @param samp_indices Vector of sample indices to be averaged.
#' @param bin Bin spectra to this m/z tolerance. If missing, do not bin at all.
#' @keywords average run experiment blank
#' @export
#' @examples
#' fname_master_blank <- average_runs(exp, samp_indices = c(1,2,3))
#'
average_runs <- function(exp, samp_indices, bin){
  # filterFile makes this quick
  exp_subset <- exp %>%
    filterFile(samp_indices)

  # load input to memory
  exp_subset <- exp_subset %>% as("MSnExp")

  # put ions in memory
  if(missing(bin)){
    ions <- exp_subset %>%
      exp2tbl()
  }else{
    ions <- exp_subset %>%
      exp2tbl(bin = bin)
  }

  # bin mz column if desired
  if(!missing(bin)){
    ions$mz <- bin_mz(ions$mz, bin)
  }

  # average the intensities
  nscans = max(ions$scan_num)
  message("averaging ion intensities")
  tbl_out <- ions %>%
    group_by(scan_num, mz) %>%
    summarise(
      rt = mean(rt),
      intensity = mean(intensity),
      ) %>%
    # important stuff to make the feature names appropriate
    group_by(scan_num) %>%
    mutate(
      samp_num = 1,
      feature = paste("F1.S", str_pad(scan_num, str_length(nscans), "0", side = "left"), sep = ""))

  # convert tbl to exp
  # but exp is for some reason not valid
  exp_new <- tbl_out %>%
    tbl2exp()

  # so for now, take contents of new exp and stuff them into old:
  exp_subset@assayData <- exp_new@assayData
  exp_subset@featureData <- exp_new@featureData
  exp_subset@processingData <- exp_new@processingData

  return(exp_subset)
}
