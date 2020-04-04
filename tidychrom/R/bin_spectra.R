#' Bin Spectra
#'
#' Convenience function takes a list of m/z values and rounds them to the nearest integer
#' with specified precision. Breaks on the 0.5.
#' @param chromdata Tibble with columns identifying spectrum x-coord (e.g. mz or wl), y-coord (intensity, abs, etc.),
#' and scan number.
#' @param x Name of x-coordinate column.
#' @param y Name of y-coordinate column.
#' @param scan Name of scan ID column.
#' @param bin_width Width of x bins.
#' @keywords bin spectrum spectra mass wavelength
#' @export
#' @examples
#' ions_binned <- ions %>% bin_spectra(x = "mz", y = "intensity", scan = "scan")
#'

bin_spectra <- function(chromdata, x = "mz", y = "intensity", scan = "scan", bin_width = 1){

  chromdata <- chromdata %>%
    # run the whole x-axis column thru cut()
    pull(x) %>%
    cut(
      breaks = seq(0.5, max(.) + bin_width, bin_width),
      labels = seq(1, max(.), bin_width)
    ) %>%
    as_tibble() %>%
    # bind it back to the source data
    bind_cols(chromdata) %>%
    # drop the unbinned x-axis
    select(-one_of(x)) %>%
    # and reassign the binned values
    rename("value" = x) %>%
    # now, sum values in the same bin
    # the _ats enable pasting of passed column names
    group_by_at(c(scan, x)) %>%
    summarise_at(
      .vars = y,
      .funs = sum
    )

  return(chromdata)
}
