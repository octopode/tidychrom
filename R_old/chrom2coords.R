# take a chromatogram, return a tibble with rt and intensity
chrom2coords <- function(chrom){
  rtime(chrom) %>%
    cbind(intensity(chrom)) %>%
    as_tibble() %>% setNames(c("rt", "intensity")) %>%
    return()
}
