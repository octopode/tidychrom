#' Average Spectra
#' Takes a tibble of grouped NORMALIZED spectrum data (multiple spectra per group)
#' returns ion-wise averaged spectra with averaged continuous statistics and nested categorical identifiers.
#' Some reserved column names get special treatment; this should become more generic in the future.
#'
#' @param spectra Tibble featuring \code{mz}/\code{wl} and \code{intensity}. Should be grouped by spectrum.
#' @keywords average consensus spectra
#' @export
#' @examples
#' spectra_norm <- spectra %>%
#'   group_by(scan) %>%
#'   normalize_spectra()
#'
normalize_spectra <- function(spectra){
  spectra_normd <- spectra %>%
    mutate(intensity = intensity/max(intensity))
  return(spectra_normd)
}
