#' Bin Spectra
#'
#' Convenience function takes a list of m/z values and rounds them to the nearest integer
#' with specified precision. Breaks on the 0.5.
#' @param chromdata Tibble with columns identifying spectrum x-coord (e.g. mz or wl), y-coord (intensity, abs, etc.),
#' and scan number
#' @param x Name of x-coordinate column
#' @param y Name of y-coordinate column
#' @param scan Name of scan ID column
#' @param bin_width Width of x bins
#' @return \code{chromdata} tibble with the x-axis binned as specified
#' @keywords bin spectrum spectra mass wavelength
#' @export
#' @examples
#' ions_binned <- ions %>% bin_spectra(x = "mz", y = "intensity", scan = "scan")
#'

bin_spectra <- function(chromdata, x = "mz", y = "intensity", scan = "scan", bin_width = 1){

  # gotta get these right beforehand
  # maybe could use some work 20200404
  breaks <- seq(0.5, max(chromdata %>% pull(x)) + bin_width, bin_width)
  labels <- lapply(seq(length(breaks) - 1), function(x){mean(breaks[[x]], breaks[[x+1]]) + 0.5})

  chromdata <- chromdata %>%
    # run the whole x-axis column thru cut()
    pull(x) %>%
    cut(
      breaks = breaks,
      labels = labels
    ) %>%
    # convert the "number factor" to actual numeric
    # handy on why as.numeric() won't work
    # https://stackoverflow.com/questions/3418128/how-to-convert-a-factor-to-integer-numeric-without-loss-of-information
    as.numeric(levels(.))[.] %>%
    as_tibble() %>%
    # bind it back to the source data
    bind_cols(chromdata) %>%
    # drop the unbinned x-axis
    select(-one_of(x)) %>%
    # reassign the binned values
    rename("value" = x) %>%
    # now, sum values in the same bin
    # the _ats enable pasting of passed column names
    group_by_at(c(scan, "rt", x)) %>%
    summarise_at(
      .vars = y,
      .funs = sum
    )

  return(chromdata)
}
