#' Separate Spectral Signals
#' Takes a grouped tibble of spectral data, returns a summary of the input
#' with an x-coordinate that is diagnostic of each of the input spectra.
#'
#' @param spectra Grouped tibble \code{mz}/\code{wl} and \code{intensity}. Should be grouped by scan.
#' @param thres_ortho Minimum orthogonality that a signal is allowed to have (1 is max)
#' @param thres_intensity Minimum relative intensity that a signal is allowed to have (0-1)
#' @param x Name of x-axis variable (\code{mz} or \code{wl}, etc.)
#' @return Summary of input tibble with diagnostic x-coordinate
#' @keywords cluster spectra cosine
#' @export
#' @examples
#' signal_ions <- conflicting_spectra %>%
#'   group_by(scan) %>%
#'   separate signals()
#'
separate_signals <- function(
  spectra,
  thres_ortho = 1,
  thres_intensity = 0.5,
  x = "mz",
  bin,
  cores = 1
  ){

  # save the original groupings for output
  group_vars_in <- group_vars(spectra)

  # bin for comparison if requested
  if(!missing(bin)){
    spectra <- spectra %>%
      bin_spectra(bin_width = bin)
  }

  signals <- spectra %>%
    ungroup() %>%
    # waiting for a better implementation of scoped complete()
    # may be able to use across() in dplyr 1.0
    complete_(c(group_vars_in, x), fill = list(intensity = 0)) %>%
    # now, estimate diagnostic power of each x [m/z]
    # 1 is a perfect score
    group_by_at(x) %>%
    mutate(
      # orthogonality
      ortho = intensity / sum(intensity),
      # "power"
      power = ortho * intensity
      ) %>%
    group_by_at(group_vars_in) %>%
    filter(ortho >= thres_ortho) %>%
    filter(intensity >= thres_intensity) %>%
    filter(power == max(power))

  if(nrow(signals) == 0){
    warning("Signals could not be separated according to criteria!")
  }

  return(signals)
}
