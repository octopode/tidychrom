#' Melt chromatographic coordinates
#'
#' Takes RT and intensity coordinates from a list of matching chromatograms (TIC, BPC, XIC etc.)
#' return a long tibble of RT, intensity, samp
#' @param chroms An XChromatograms object.
#' @keywords melt chrom chromatogram coord coordinates
#' @export
#'
melt_coords <- function(chroms){
  chroms_coords <- unlist(chroms) %>%
    # need to have samps in phenoData
    setNames(chroms$samps) %>%
    lapply(., chrom2coords) %>%
    as_tibble()
  # melt all the rt/intensity pairs
  chroms_coords.m <- gather_listoftbl(chroms_coords)
  return(chroms_coords.m)
}
