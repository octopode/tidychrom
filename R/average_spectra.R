#' Average Spectra
#' Takes a tibble of grouped NORMALIZED spectrum data (multiple spectra per group)
#' returns ion-wise averaged spectra with averaged continuous statistics and nested categorical identifiers.
#' Some reserved column names get special treatment; this should become more generic in the future.
#'
#' @param spectra Grouped tibble with \code{mz}/\code{wl} and \code{intensity}. Should be grouped by clusters to be averaged.
#' @param x Name of x-axis variable (\code{mz} or \code{wl}, etc.)
#' @keywords average consensus spectra
#' @export
#' @examples
#' spectra_cons <- spectra_clust %>%
#'   group_by(clust_id) %>%
#'   average_spectra()
#'
average_spectra <- function(spectra, x = "mz"){

  cluster_cols <- group_vars(spectra)
  spectra_avg <- spectra %>%
    # count spectra per cluster by counting basepeaks
    mutate(n_spectra = sum(intensity == max(intensity))) %>%
    # group along x axis, in addition to by cluster
    group_by_at(c(cluster_cols, x, "n_spectra")) %>%
    # manual calculation of mean saves time by avoiding expand() or complete()
    summarize(intensity = sum(intensity) / first(n_spectra)) %>%
    select(-n_spectra) %>%
    # restore original grouping
    group_by_at(cluster_cols)

  return(spectra_avg)
}
