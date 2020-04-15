#' Convert MSnExp to Tibble
#'
#' Take an MSNExp and convert it losslessly (but exclusive of metadata)
#' into a tbl with columns rt, mz, and intensity (the tidychrom-ms format.)
#' @param exp MSNExp object
#' @param bin Bin spectra to this m/z tolerance. If missing, do not bin at all.
#' @keywords experiment convert
#' @export
#' @examples
#' ions <- exp2tbl(exp)
#'

exp2tbl <- function(exp, bin){

  # Third and Best Way 20200402
  # using rtime(), mz(), and intensity() functions
  # these appear to harness MSnbase's parallel I/O backend
  ions <- mapply(mz(exp), intensity(exp), FUN=cbind) %>%
    mapply(., rtime(exp), FUN=cbind) %>%
    lapply(., as_tibble) %>%
    # would be nice to have a more standard solution for this
    gather_listoftbl(key = "feature", parallel = T) %>%
    transmute(mz = V1, intensity = V2, rt = V3, feature = feature) %>%
    group_by(feature) %>% # need for the below mutates to work efficiently
    mutate(
      samp_num = feature2samp(feature),
      scan_num = feature2scan(feature)
    )

  # bin mz column if desired
  if(!missing(bin)){
    ions$mz <- bin_mz(ions$mz, bin)
  }

  # Second Way ~5 min/run = 0.08 s/scan
  ## get feature names
  ## m/z-bin experiment if desired
  ## mz(.) %>% unlist() %>% max() ALSO TOO SLOW, doubles overall exec time
  #if(!missing(bin)){
  #  exp <- exp %>% bin(bin)
  #}

  #features <- exprtime2tbl(exp)
#
  ## disk read ops are parallelized using library(pbapply)
  ## Slow and power-hungry, but probably as fast as things will go on an existing R backend.
  ## At least there's a progress bar to keep you sane.
  ## Calculate the number of cores
  #num_cores <- detectCores()
  ## Initiate cluster with current environment
  #clust <- makeCluster(num_cores, type="FORK")
  #message("loading spectra into memory")
  #peaklists <- pblapply(cl = clust, features$feature, function(x){exp_subset[[x]] %>% spec2peaklist()})
  #stopCluster(clust)
  #names(peaklists) <- features$feature
#
  #tbl_out <- gather_listoftbl(peaklists, key = "feature", parallel = T) %>%
  #  left_join(features, by = "feature")

  # NOPE, TOO SLOW (~2 s/scan)
  ## using a loop and df indexing here to protect memory
  #tbl_out <- NULL
  #for(i in seq(nrow(features))){
  #  message(features[i,] %>% pull(feature))
  #  tbl_out <- tbl_out %>%
  #    bind_rows(
  #      exp[[features[i,] %>% pull(feature)]] %>%
  #        spec2peaklist() %>%
  #        mutate(
  #          samp_num = features[i,] %>% pull(samp_num),
  #          rt = features[i,] %>% pull(rtime)
  #        )
  #      )
  #}

  return(ions)
}
