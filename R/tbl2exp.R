#' Convert Tibble to MSnExp
#'
#' Pack a tibble into an MSNExp object.
#' @param ions A tibble with columns feature, mz, and intensity
#' @keywords experiment convert
#' @export
#' @examples
#' expt <- tbl2exp(ions)
#'
tbl2exp <- function(ions){
  # group ions by the scan identifier
  ions <- ions %>%
    group_by(feature)# %>%
    # trim the fat
    #select(mz, intensity, rt, feature)

  # step 1: put ions into spectra
  message("consolidating spectra")
  # break tibble into scans
  peaklists <- ions %>%
    group_split(keep = F)
  names(peaklists) <- ions %>%
    group_keys() %>%
    unlist()

  # Calculate the number of cores
  num_cores <- detectCores()
  # Initiate cluster with current environment
  clust <- makeCluster(num_cores, type="FORK")
  # failsafe
  on.exit(stopCluster(clust))
  # takes maybe 5 s/sample
  spectra <- pblapply(cl = clust, peaklists, peaklist2spec)
  stopCluster(clust)
  # step 2: put spectra into experiment
  # environment containg spectra
  spectra_env <- list2env(spectra)
  # AnnotatedDataFrame with the right rownames
  features <- ions %>%
    select(feature, samp_num, scan_num, rt) %>%
    distinct(feature, .keep_all = T) %>%
    dplyr::rename(
      fileIdx = samp_num,
      spIdx = scan_num,
      retentionTime = rt
      ) %>%
    rowid_to_column(var = "spectrum") %>%
    mutate(
      seqNum = spIdx,
      acquisitionNum = spIdx,
      msLevel = ifelse("mslevel" %in% names(ions), mslevel, 1),
      spectrumId = paste("scan=", spIdx, sep=""),
      ## seems this is all so much garbage
      #originalPeaksCount = 1,
      #totIonCurrent = 1E7, #placeholder
      #basePeakMZ = -1,
      #basePeakIntensity = -1,
      #collisionEnergy = -1,
      #ionisationEnergy = -1,
      #highMZ = -1,
      #precursorScanNum = -1,
      #precursorMZ = -1,
      #precursorCharge = -1,
      #precursorIntensity = -1,
      #mergedScan = -1,
      #mergedResultScanNum = -1,
      #mergedResultStartScanNum = -1,
      #mergedResultEndScanNum = -1,
      #injectionTime = -1,
      #centroided = NA,
      #ionMobilityDriftTime = -1,
      #isolationWindowTargetMZ = NA,
      #isolationWindowLowerOffset = NA,
      #isolationWindowUpperOffset = NA,
      #scanWindowLowerLimit = NA,
      #scanWindowUpperLimit = NA,
      #polarity = NA,
    ) %>%
    #select(feature) %>%
    column_to_rownames(var = "feature") %>%
    as(., "AnnotatedDataFrame")
  # MSnProcess object
  # get filenames, or placeholders for them, as efficiently as possible
  if(F & "samp" %in% names(ions)){
    files <- unique(ions$samp) # this creates a mismatch error when writing, so disabled for now
  }else{
    if("samp_num" %in% names(ions)){
      files <- unique(ions$samp_num) %>%
        as.character()
    }else{
      files <- ions %>%
        select(feature) %>%
        distinct() %>%
        transmute(samp_num = feature2samp(feature)) %>%
        pull(samp_num) %>%
        as.character() %>%
        unique()
    }
  }
  process <- new(
    "MSnProcess",
    # length needs to match the number of indices after "F"
    files = files
  )
  # dumb placeholder
  # may be unnecessary
  pheno <- files %>%
    as_tibble() %>%
    as.data.frame() %>%
    as("AnnotatedDataFrame")
  exp <- new(
    "MSnExp",
    assayData = spectra_env,
    featureData = features,
    processingData = process,
    phenoData = pheno
             )
  return(exp)
}
