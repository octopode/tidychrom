## CHECK PEAK OVERLAPS
# for assessment of RT window and cosine threshold

# take a matched-peaks tibble containing gg objects and rt_std column
# plot integrals for the specified standard peak index (in RT order)
plot_std_peak_match <- function(peaks_matched_tbl, std_peak_index){

  chprs <- peaks_matched_tbl %>%
    arrange(rt_std) %>%
    mutate(rt_std = as_factor(rt_std)) %>%
    group_by(rt_std) %>%
    # select the index of the standard peak
    filter(rt_std == unique(.$rt_std)[std_peak_index]) %>%
    # pull the lists (pairs) of ggproto objects
    pull(gg)

  # plot all the above objects together
  if(length(chprs) > 1){
    pal = c(palette(rainbow(length(chprs)))[2:length(chprs)], "black")
  }else{
    pal = c("black")
  }

  gg_chrom <- ggplot() +
    chprs +
    #scale_color_brewer(palette = "Dark2") +
    #scale_fill_brewer(palette = "Dark2") +
    scale_color_manual(values = pal) +
    scale_fill_manual(values = pal) +
    theme_pubr() +
    theme(legend.position = "none")
  return(gg_chrom)
}
