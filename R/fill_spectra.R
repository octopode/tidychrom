#' Fill spectra (with 0s)
#'
#' Take tibble of full-scan chromatographic data grouped by scan,
#' fill every missing x-coord (m/z or wl) in each scan with 0.
#' This is important prior to single-ion or -wavelength peak detection and integration.
#' @param chromdata Tibble with columns \code{rt}, \code{mz}/\code{wl} and \code{intensity}, pre-grouped by peak (e.g. by ROI.)
#' @param x Name of spectrum x-coordinate (e.g. \code{mz} or \code{wl}.)
#' @return The input tibble with 0-intensity observations added.
#' @keywords fill spectra zero zeroes gaps
#' @export
#' @examples
#'
#' chromdata_filled <- chromdata %>%
#'  group_by(scan_rt) %>%
#'  fill_spectra

fill_spectra <- function(chromdata, x = "mz"){

  group_vars_in <- group_vars(chromdata)
  # this is a workaround to scoped nesting()
  chromdata$group_id <- group_indices(chromdata)

  chromdata_filled <- chromdata %>%
    ungroup() %>%
    complete_(c("group_id", x), fill = list(intensity = 0)) %>%
    # replace NAs in original grouping vars
    #group_by(group_id) %>%
    #mutate_at(c(group_vars_in) = )
    select(-group_id) %>%
    group_by_at(group_vars_in) %>%

  return(chromdata_filled)
}
