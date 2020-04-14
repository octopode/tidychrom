# Helper function: take two (binned) spectra as tibbles
# Join them and calculate the dot product.
cosine_spectra <- function(spec_1, spec_2, x = "mz"){
  specs <- full_join(spec_1, spec_2, by = x) %>%
    mutate(
      intensity.x = ifelse(is.na(intensity.x), 0, intensity.x),
      intensity.y = ifelse(is.na(intensity.y), 0, intensity.y)
    )
  # cosine similarity function used here is lifted from
  # Fridolin Wild's LSA package (https://cran.r-project.org/web/packages/lsa/)
  return(cosine(specs$intensity.x, specs$intensity.y)[[1]])
}
