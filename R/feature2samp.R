#' Extract Sample Index from MSNbase Feature String
#'
#' @param feature string of form "Fxx.Sxxxx"
#' @keywords extract feature samp sample
#' @export
#' @examples
#' feature2samp("F01.S0004")
#' >[1] 1
#'
feature2samp <- function(feature){
  samp_num = strsplit(feature, "\\.") %>%
    .[[1]] %>% .[1] %>%
    substr(., 2, str_length(.)) %>%
    as.numeric()
  return(samp_num)
}
