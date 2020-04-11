#' Write Tidy Mass Chromatogram to File
#'
#' @param ions A tibble with columns scan or scan_num, mz, and intensity
#' @keywords write tibble mz
#' @export
#' @examples
#' ions %>% tbl2exp(file = "mychrom.mzxml")
#'
write_tidymass <- function(ions, file){
  #
  if("scan_num" %in% names(ions)){
    ions <- ions %>%
      mutate(scan = scan_num)
  }
  # convert spectra to matrix format
  spectra <- ions %>%
    group_by(scan) %>%
    group_split() %>%
    lapply(., function(x){
      x %>%
        select(mz, intensity) %>%
        as.matrix()
    })
  # write the headers
  headers <- ions %>%
    #select(scan, rt) %>%
    distinct(scan, .keep_all = T) %>%
    # scan numbers must start with 1 and be consecutive for save to work
    ungroup() %>%
    mutate(scan = seq(length(.$scan))) %>%
    dplyr::rename(
      retentionTime = rt
    ) %>%
    mutate(
      seqNum = scan,
      acquisitionNum = scan,
      msLevel = ifelse("mslevel" %in% names(ions), mslevel, 1),
      spectrumId = paste("scan=", scan, sep=""),
      peaksCount = 1,
      totIonCurrent = 1E7, #placeholder
      basePeakMZ = -1,
      basePeakIntensity = -1,
      collisionEnergy = -1,
      ionisationEnergy = -1,
      lowMZ = -1,
      highMZ = -1,
      precursorScanNum = -1,
      precursorMZ = -1,
      precursorCharge = -1,
      precursorIntensity = -1,
      mergedScan = -1,
      mergedResultScanNum = -1,
      mergedResultStartScanNum = -1,
      mergedResultEndScanNum = -1,
      filterString = "",
      injectionTime = -1,
      centroided = NA,
      ionMobilityDriftTime = -1,
      isolationWindowTargetMZ = 0,
      isolationWindowLowerOffset = 0,
      isolationWindowUpperOffset = 0,
      scanWindowLowerLimit = 0,
      scanWindowUpperLimit = 0,
      polarity = -1,
    ) %>%
    select(-scan) %>%
    as.data.frame()

  #this totally works! The header data is bunk, but can fix that later.
  mzR::writeMSData(spectra, headers, file=file)
}
