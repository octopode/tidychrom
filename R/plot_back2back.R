#' Plot Spectra Back-to-Back
#'
#' Take candidate spectra, return a summary of the best matches against a set of master spectra.
#' @param ions Grouped tibble containing mz and intensity.
#' @param colors Character vector of colors in which to plot the spectra.
#' @param label_pctile Upper percentile of peaks to label with m/z or wl.
#' @param alpha_unmatched Alpha value to fade unmatched peaks.
#' @keywords plot spectra match cosine
#' @export
#' @examples
#' scan_cosines <- cosine_by_roi_scan(spectra_candidate, spectra_master)
#'

plot_back2back <- function(ions, pal = c("blue", "red"), x = "mz", label_pctile = 5, alpha_unmatched = 0.8){
  # group index column
  ions$gidx <- group_indices(ions)
  # make plot object
  data <- ions %>%
    # silently get just the first two groups
    filter(
      (gidx <= 2) &
        !is.na(intensity) &
        intensity != 0
    ) %>%
    mutate(
      intensity = (intensity / max(intensity)), # normalize intensity
      mzlabel = ifelse(intensity > quantile(intensity, (1 - (label_pctile / 100))), mz, ""),
      intensity = ifelse(gidx == 1, intensity, -1*intensity),
      ) %>%
    group_by_at(x) %>% # group by x to check matches
    mutate(alpha = ifelse(n() > 1, 1, alpha_unmatched))
  #print(data %>% pull(alpha))
  #print(data %>% select(mz, scan, alpha))

    # set up a list of gg objects for return
    gg <- list(
      geom_col(data = data,
        aes(
          x = get(x),
          y = intensity,
          #fill = as.factor(.[[group_vars(.)]])
          fill = as.factor(gidx),
          alpha = alpha
        ),
      position = "identity"
      ),
      geom_text(data = data,
        aes(
          x = get(x),
          y = intensity,
          label = mzlabel,
          alpha = alpha
          ),
        size = 2
        ),
      scale_fill_manual(values = pal),
      scale_y_continuous(breaks = seq(-1,1,0.25), labels = abs(seq(-1,1,0.25))),
      xlab("m/z"),
      ylab("relative abundance"),
      theme_pubr(),
      theme(legend.position = "none")
    )

  return(gg)

}
