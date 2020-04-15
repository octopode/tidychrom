#' Compare Spectra by ROI and Scan
#'
#' Take candidate spectra, return a summary of the best matches against a set of master spectra.
#' A prepackaged function to identify a known compound within an ROI in targeted analysis.
#' @param spectra_candidate Tibble containing full scans and an \code{roi} column. Can be multiple scans/ROI.
#' @param spectra_master Tibble containing standard spectra, also identified with an roi number. Should be only 1 scan/ROI.
#' @return a summary of \code{spectra_candidate} with only the best match for each ROI and additional column \code{cos}.
#' @keywords cosine roi scan targeted
#' @export
#' @examples
#' scan_cosines <- cosine_by_roi_scan(spectra_candidate, spectra_master)
#'

# NTS 20200413 JRW: Should be refactored to work on a pre-grouped (e.g. by ROI) tibble,
# instead of having rois_scans$roi hardcoded
# the tidy solution at the bottom can be reinstated with column renaming,
# but the current version is nice because it can be parallelized.

cosine_by_roi_scan <- function(spectra_candidate, spectra_master, viz = F){

  # candidate ROI and scan combinations
  rois_scans <- spectra_candidate %>%
    select(roi, scan, rt) %>%
    distinct()

  # get cosine scores for each of those combinations
  cosines <- mapply(
      rois_scans$roi,
      rois_scans$scan,
      FUN = function(x, y){
        cosine_spectra(
          spectra_candidate %>%
            filter(
              roi == x,
              scan == y
            ),
          spectra_master %>%
            filter(roi == x)
        )
      })

  rois_scans$cos <- cosines

  return(rois_scans)

  # Would be nice if this tidy solution worked, but cos values all come out the same :/
  #spectra_candidate %>%
  #  group_by(roi, scan) %>%
  #  summarize(
  #    cos = match_spectrum(
  #      spectra_candidate %>%
  #        filter(
  #          roi == roi,
  #          scan == scan
  #        ),
  #      spectra_master %>%
  #        filter(roi == roi)
  #    )
  #  )
}
