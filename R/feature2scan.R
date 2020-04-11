#' Extract Sample Index from MSNbase Feature String
#'
#' @param feature string of form "Fxx.Sxxxx"
#' @keywords extract feature samp sample
#' @export
#' @examples
#' feature2scan("F01.S0004")
#' >[1] 4
#'
feature2scan <- function(feature){
  scan_num = strsplit(feature, "\\.") %>%
    .[[1]] %>% .[2] %>%
    substr(., 2, str_length(.)) %>%
    as.numeric()
  return(scan_num)
}
