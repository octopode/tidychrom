# take a tbl column of spectra and plot them all stacked up
# label the top label_pctile of peaks
plot_spectra <- function(spectra, master_index = 1,
                         # default palette is bad, but takes infinite spectra
                         pal = c("black", palette(rainbow(nrow(spectra)))[2:nrow(spectra)]), label_pctile = 5){
  peaklists <- spectra %>%
    # make sure all the spectrum objs are unpacked
    unlist() %>%
    lapply(., spec2peaklist)
  if("samp" %in% colnames(spectra)){
    names(peaklists) <- spectra$samp
  }
  peaklists <- peaklists %>%
    gather_listoftbl() %>%
    mutate(
      # put the master sample first
      samp = fct_relevel(samp, levels = c(unique(.$samp)[master_index], sort(unique(.$samp)[setdiff(seq(length(samp)), master_index)])))
    ) %>%
    group_by(samp) %>%
    mutate(
      # normalize intensity
      rel_abund = intensity / max(intensity),
      # and label the upper xtile of peaks
      mzlabel = ifelse(rel_abund > quantile(rel_abund, (1 - (label_pctile / 100))), mz, "")
    )
  fig <- peaklists %>% ggplot(aes(x = mz, y = rel_abund, fill = samp, label = mzlabel)) +
    geom_col() +
    geom_text(size = 2, color = "black", nudge_y = 0.05) +
    facet_grid(rows = vars(samp)) +
    scale_fill_manual(values = pal) +
    xlab("m/z") +
    ylab("relative abundance") +
    theme_pubr() +
    theme(legend.position = "none")
  return(fig)
}
