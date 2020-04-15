#' Cosine Similarity of Two Spectra
#'
#' Takes two spectra in tibble form, calculates their cosine similarity score.
#' In other words, performs a full join on the x-axis, fills with 0s, then gets cosine.
#' @param spec_1 Spectrum 1 in tibble form, with x-axis named as passed
#' @param spec_2 Spectrum 2 in tibble form, with x-axis named as passed
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
cosine_spectra <- function(spec_1, spec_2, x = "mz"){
  specs <- full_join(spec_1, spec_2, by = x) %>%
    # peaks present in 1 spectrum should be 0s in the other
    mutate(
      intensity.x = ifelse(is.na(intensity.x), 0, intensity.x),
      intensity.y = ifelse(is.na(intensity.y), 0, intensity.y)
    ) %>%
    # constrains comparison to the tighter mz window between the two spectra
    filter(
      # on the low end
      mz >= max(min(spec_1[[x]]), min(spec_2[[x]])),
      # on the high end
      # this one is problematic
      #mz <= min(max(spec_1[[x]]), max(spec_2[[x]])),
      )
  # cosine similarity function used here is lifted from
  # Fridolin Wild's LSA package (https://cran.r-project.org/web/packages/lsa/)
  return(cosine(specs$intensity.x, specs$intensity.y)[[1]])
}
