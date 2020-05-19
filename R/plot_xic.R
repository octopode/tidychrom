#' Plot Extracted Ion Chromatogram
#'
#' Take a table of rt and intensity coordinates, return a list of gg objects
#' representing a chromatogram with baseline-corrected area filled.
#' @param coords rt and intensity coords
#' @keywords plot spectra match cosine
#' @export
#' @examples
#' xic_gg <- plot_xic(chromdata %>% filter(roi == 2))
#'
plot_xic <- function(coords, alpha = 0.4){
  coords <- coords %>%
    arrange(rt)
  gg = list(
    # plot the line
    geom_line(data = coords, aes(x = rt, y = intensity)),
    # and the baselined area
    geom_polygon(data = coords, aes(x = rt, y = intensity), alpha = alpha),
    #geom_area(data = coords, aes(x = rt, y = intensity), alpha = alpha),
    #geom_area(data = bind_rows(head(coords, 1), tail(coords, 1)), aes(x = rt, y = intensity), fill = "white"),
    theme_pubr()
  )
  return(gg)
}
