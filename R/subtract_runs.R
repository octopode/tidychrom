#' Subtract Chromatographic Runs
#'
#' Take an MSNExp and two sample indices. Subtract the second index from the first,
#' save the resultant CDF, and return that filename.
#' @param exp MSNExp object; both runs need to be in the experiment.
#' @param samp_index_A Index of first sample, where result = A - B.
#' @param samp_index_B Index of second sample (e.g. blank), where result = A - B.
#' @param bin Bin spectra to this m/z tolerance. If missing, do not bin at all.
#' @keywords subtract spectrum spectra baseline blank
#' @export
#' @examples
#' exp_blanked <- subtract_blank_exp(exp, blank_i = c(1,2))
#'
subtract_runs <- function(exp, samp_index_A, samp_index_B, bin){

  # get CDF filenames
  filenames_raw <- exp %>% fileNames()
  # sample, rt -> feature mapping tbl
  # this call takes awhile
  features <- exp %>% exprtime2tbl()

  # get indices of non-blank samples. Needed?
  samp_i <- setdiff(features$samp_num %>% unique(), blank_i)

  # get all the blank spectra into a tibble
  spectra_blank <- features %>%
    filter(samp_num %in% blank_i)# %>%
    # this extraction is slow and memory-hogging
    mutate(spectrum = list(exp[[feature]]))

  # average coincident blank spectra
  spectra_blank_avg <- spectra_blank %>%
    group_by(rtime) %>%
    summarise(
      spectrum = average_spectra(
        # to unpack each spectrum
        lapply(spectrum, unlist)
      )
    )

  return(spectra_blank_avg)
}
