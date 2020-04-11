#' Read Tidy Mass Chromatogram from File
#'
#' @param file A tibble with columns scan or scan_num, mz, and intensity
#' @keywords write tibble mz
#' @export
#' @examples
#' ions %>% tbl2exp(file = "mychrom.mzxml")
#'
read_tidymass <- function(file){

  mzr <- openMSfile(file)

  # different routines for single, multiscan files
  # could be cleaner but it works for now
  if(length(header(mzr)$retentionTime) > 1){
    ions <- mzr %>%
      spectra() %>%
      setNames(., header(mzr)$retentionTime) %>%
      gather_listoftbl(., key = "rt", index = "scan") %>%
      mutate(rt = as.numeric(rt)) %>%
      # when loading from some filetypes, mz and intensity columns are "V1", "V2"
      set_names(c("mz", "intensity", names(.)[3:length(names(.))]))
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
