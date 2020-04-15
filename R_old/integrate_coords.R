# integrate a table of peaks; wrapper for auc()
integrate_coords <- function(xic_peaks, xics_coords){
  xic_peaks <- xic_peaks %>%
    rowwise() %>%
    mutate(
      # use of logical subsetting here is intentional
      area = list(auc(xics_coords[xics_coords$samp == samp,], rt))
    ) %>%
    # turn the list column into separate columns
    bind_cols(do.call(rbind, .$area) %>% as_tibble()) %>%
    select(-area)
  return(xic_peaks)
}
