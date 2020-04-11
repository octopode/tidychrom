plot_spectrum <- function(spectrum, label_pctile = 5){
  gg <- spectrum %>%
    spec2peaklist() %>%
    mutate(
      rel_abund = intensity / max(intensity),
      mzlabel = ifelse(rel_abund > quantile(rel_abund, (1 - (label_pctile / 100))), mz, "")
    ) %>%
    ggplot(aes(x = mz, y = intensity, label = mzlabel)) +
    geom_col() +
    geom_text(size = 2) +
    theme_pubr()
  return(gg)
}
