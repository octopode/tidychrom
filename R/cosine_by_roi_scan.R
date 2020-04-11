#' Compare Spectra by ROI and Scan
#'
#' Take candidate spectra, return a summary of the best matches against a set of master spectra.
#' @param spectra_candidate Tibble containing full scans and an roi column. Can be multiple scans/ROI.
#' @param spectra_master Tibble containing standard spectra, also identified with an roi number. Should be only 1 scan/ROI.
#' @keywords spectra cosine dot product
#' @export
#' @examples
#' scan_cosines <- cosine_by_roi_scan(spectra_candidate, spectra_master)
#'

cosine_by_roi_scan <- function(spectra_candidate, spectra_master){

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

# Helper function: take two (binned) spectra as tibbles
# Join them and calculate the dot product.
cosine_spectra <- function(spec_1, spec_2, x = "mz"){
  specs <- inner_join(spec_1, spec_2, by = x)
  # cosine similarity function used here is lifted from
  # Fridolin Wild's LSA package (https://cran.r-project.org/web/packages/lsa/)
  return(cosine(specs$intensity.x, specs$intensity.y)[[1]])
}
