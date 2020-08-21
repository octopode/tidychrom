#' Area Under the Curve (AUC)
#'
#' Summarize grouped chromatographic data (windowed chromatograms),
#' returning the RT window as well as overall and baseline-corrected integrals.
#' @param chromdata Tibble with columns \code{rt} and \code{intensity}, pre-grouped by peak (e.g. by ROI.)
#' @return A summary of the input tibble with additional columns:
#' \code{rt_min}, \code{rt_max},
#' \code{into} (integral overall),
#' \code{intb} (integral baseline-corrected)
#' @keywords area integral peak curve
#' @export
#' @examples
#'
#' areas_rois <- xics %>%
#'  group_by(roi) %>%
#'  auc()

auc <- function(chromdata){
  areas <- chromdata %>%
    arrange(rt) %>%
    summarise(
      rt_min = min(rt),
      rt_max = max(rt),
      into = mean(intensity) * (rt_max - rt_min),
      intb = into - (mean(dplyr::first(intensity), dplyr::last(intensity)) * (rt_max - rt_min)),
      intensity = max(intensity)
    ) %>%
    # join back to input grouping vars
    # these can hold metadata such as the mz/wl, cosine, peak rt and intensity, etc.
    left_join(
      chromdata %>%
        summarize(),
      by = dplyr::group_vars(chromdata)
      )

  return(areas)
}
