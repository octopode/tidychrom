#' Read Tidy Mass Chromatogram from File
#'
#' @param file Path to an mzR-compatible MS1 file
#' @return A tibble with columns \code{scan}, \code{rt}, \code{mz}, and \code{intensity}
#' @keywords write tibble mz
#' @export
#' @examples
#' ions <- read_tidymass(file = "mychrom.mzxml")
#'
# NTS 20200415 JRW: add MS level filtering

read_tidymass <- function(file){

  mzr <- openMSfile(file)
  rts <- header(mzr)$retentionTime

  # different routine for single- than for multiscan file
  # could be cleaner but it works for now
  if(length(rts) > 1){
    ions <- mzr %>%
      spectra() %>%
      mapply(., rts, FUN=function(spec, rt_val){
        spec %>%
          as_tibble() %>%
          mutate(rt = as.numeric(rt_val))
      }, SIMPLIFY = FALSE) %>%
      do.call(rbind, .) %>%
      #dplyr::rename(mz = V1, intensity = V2) %>%
      arrange(rt) %>%
      group_by(rt)
    ions$scan <- group_indices(ions)
    # set irrespective of original names
    colnames(ions) <- c("mz", "intensity", "rt", "scan")
  }else{
    # implementation of the above this way (with group_indices) might speed it up
    # would need to use a rep() call or something to replicate the header over ions
    ions <- mzr %>%
      spectra() %>%
      cbind(., header(mzr)$retentionTime) %>%
      as_tibble() %>%
      set_names(c("mz", "intensity", "rt")) %>%
      mutate(scan = group_indices(., rt))
  }

  return(ions)
}
