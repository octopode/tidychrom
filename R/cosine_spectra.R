#' Pairwise cosine similarity scores
#'
#' Takes two spectra in tibble form, calculates their cosine similarity score.
#' @param spec_1 Spectrum 1 in tibble form, with x-axis named as passed
#' @param spec_2 Spectrum 2 in tibble form, with x-axis named as passed
#' @param thres_pks Minimum number of matching peaks, else return 0
#' @param x Name of x-coordinate column (e.g. \code{mz} or \code{wl})
#' @return cosine similarity score of the two spectra (dbl)
#' @keywords cosine similarity score spectra
#' @export
#' @examples
#' cos <- cosine_spectra(ions %>% filter(scan == 114), ions %>% filter(scan == 114))
#' > 0.94
#'

# Helper function: take two (binned) spectra as tibbles
# Join them and calculate the dot product.
cosine_spectra <- function(spec_1, spec_2, thres_pks = 2, x = "mz"){
  specs <- inner_join(spec_1, spec_2, by = x)
  if(nrow(specs) < thres_pks){
    return(0)
  }else{
    # cosine similarity function used here is lifted from
    # Fridolin Wild's LSA package (https://cran.r-project.org/web/packages/lsa/)
    return(cosine(specs$intensity.x, specs$intensity.y)[[1]])
  }
}
