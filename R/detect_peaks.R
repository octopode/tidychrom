#' Detect Peaks
#' Take chromatographic data, return summary of peaks above a hard threshold, SNR threshold, or biggest n peaks.
#' Data should be input as one or more groups of \code{rt}: \code{intensity}, without gaps in RT.
#' NOTE: improper grouping (e.g. by \code{rt}) will cause this function to return an empty tbl.
#'
#' @param chromdata chromatographic data, preferable single-ion or -wavelength
#' @param thres_snr intensity SNR threshold for peak detection
#' @param thres_intensity Hard intensity threshold for peak detection
#' @param int_min Minimum nonzero intensity, used to reasonably calculate SNR
#' @param n number of peaks to return
#' @return Peak rows for each group of the passed tibble, with a coords column containing the peak domain
#' @keywords detect peaks
#' @export
#' @examples
#' peaks <- chromdata %>%
#' get_peaks(thres_snr = 20)
#'
#'
detect_peaks <- function(chromdata, thres_snr = 0, thres_intensity = 0, n = Inf, int_min = 1, target_rt){

  # make sure data are in order
  chromdata <- chromdata %>%
    arrange(rt)

  # filter for peaks
  peaks <- chromdata %>%
    filter((intensity > lag(intensity)) & (intensity > lead(intensity))) %>%
    mutate(peak = TRUE)

  # filter for valleys
  valleys <- chromdata %>%
    filter(
      (intensity <= lag(intensity)) & (intensity <= lead(intensity)) |
        # put fenceposts at the edges
        (rt == first(rt)) | (rt == last(rt))
    ) %>%
    mutate(peak = FALSE) # i.e. it's a valley

  # calc SNR for each peak and filter using threshold
  peaks_out <- peaks %>%
    bind_rows(valleys) %>%
    arrange(rt) %>%
    mutate(
      # don't bother calculating for valleys
      snr = ifelse(peak == T, intensity / pmax(int_min, pmin(lead(intensity), lag(intensity))), NA),
      rt_min = ifelse(peak == T, lag(rt), NA),
      rt_max = ifelse(peak == T, lead(rt), NA)
    ) %>%
    filter(
      peak == T,
      snr >= thres_snr
    ) %>%
    select(-peak) %>%
    # filter by intensity threshold
    filter(intensity >= thres_intensity) %>%
    # get the top n
    arrange(desc(intensity)) %>%
    slice(1:min(n, n())) %>%
    # and put them in RT order again
    arrange(rt)

  return(peaks_out)
}
