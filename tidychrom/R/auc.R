#' Area Under the Curve (AUC)
#'
#' take molten chrom coords and an RT
#' return info about the peak, including
#' i.e. integrate a single peak
#' @param coords Tibble with columns \code{rt} and \code{intensity}
#' @param rt_peak Retention time of the focal peak in units matching coords
#' @param plot_intb Logical specifying whether to return ggplot objects
#' @return a named list of:
#' \code{rtmin}, \code{rtmax},
#' \code{into} (integral overall),
#' \code{into}, (integral overall),
#' \code{intb}, (integral baseline-corrected),
#' \code{gg}, a vector of gg objects depicting the integrated peak
#' @keywords area integral peak curve
#' @export
#' @examples
#'
#' peaks_matched <- annot_from_db(peaks_matched, "~/msdb/", "methyl")

auc <- function(coords, rt_peak, plot_intb = T){
  #print(coords) #TEST the whole window
  rtmin <- coords %>%
    filter((rt < rt_peak) & ((intensity <= lag(intensity)) | is.na(lag(intensity)) | rt == min(rt)) & (intensity < lead(intensity))) %>%
    filter(rt == max(rt)) %>%
    pull(rt)
  rtmax <- coords %>%
    filter((rt > rt_peak) & (intensity < lag(intensity)) & ((intensity <= lead(intensity)) | is.na(lead(intensity)) | rt == max(rt))) %>%
    filter(rt == min(rt)) %>%
    pull(rt)
  # failsafe in case there is no peak (trace is monotonic within window)
  if(length(rtmin) == 0){
    rtmin <- rtmax
  }
  if(length(rtmax) == 0){
    rtmax <- rtmin
  }
  coords <- coords %>%
    filter((rtmin <= rt) & (rt <= rtmax)) %>%
    arrange(rt)
  #print(coords) #TEST the integrated points
  delta_rt <- (coords[["rt"]][nrow(coords)] - coords[["rt"]][1])
  into <- mean(coords[["intensity"]]) * delta_rt
  intb <- into - (mean(c(coords[["intensity"]][1], coords[["intensity"]][nrow(coords)])) * delta_rt)
  gg <- c()
  if (plot_intb == T){
    #gg <- coords %>%
    #  ggplot(aes(x = rt, y = intensity)) +
    #  geom_line(aes(color = samp)) +
    #  geom_area(data = . %>% filter((rt >= rtmin) & (rt <= rtmax)), aes(fill = samp), alpha = 0.4) +
    #  theme_pubr()
    gg <- c(
      geom_line(data = coords, aes(color = samp, x = rt, y = intensity)),
      geom_area(data = coords %>% filter((rt >= rtmin) & (rt <= rtmax)), aes(fill = samp, x = rt, y = intensity), alpha = 0.4)
      # should add a geom_area here that matches bkgd color, to illustrate BL correction
    )
  }
  return(list(rtmin = rtmin, rtmax = rtmax, into = into, intb = intb, gg = gg))
}
