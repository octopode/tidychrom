#' Filter Chromatogram data to Regions of Interest (ROIs)
#'
#' @param chromdata A tibble of chromatography data
#' @param peaks A tibble of peak RTs on which to center the ROIs. Should also contain an x-axis column (e.g. mz or wl)
#' @param roi_width Width of ROIs to generate, in seconds
#' @param x Name of the x-axis column, like "mz" or "wl"
#' @param cores Number of cores to use for parallel ops on the input list. Only beneficial in huge tables (1E7s of rows)
#' @keywords roi region interest
#' @export
#' @examples
#' rois_sample <- chromdata_sample %>% filter_rois(peaks_master, roi_width = 4)
#'
filter_rois <- function(chromdata, peaks, roi_width, x = "mz", cores = 1){

  message("extracting regions of interest")

  if(cores > 1){
    clust <- makeCluster(cores, type="FORK")
    on.exit(stopCluster(clust)) # important!
  }else{
    clust <- NULL
  }

  chromdata <- chromdata %>%
    # the grouping scheme here appears to have no effect on speed - 20200406 JRW
    #group_by_at(c("rt", x))
    rowwise()
  # for some reason, mapply won't do the job here, but lapply with the index will!
  rts_peaks <- peaks %>% pull(rt)
  xs_peaks <- peaks %>% pull(x)
  chromdata_rois <- pblapply(
    cl = clust,
    seq(nrow(peaks)),
    function(i){
      chromdata %>%
        filter(
          (rt > rts_peaks[i] - roi_width/2) &&
            (rt < rts_peaks[i] + roi_width/2) &&
            (mz == xs_peaks[i])
        )
    }
  ) %>%
    setNames(seq(nrow(peaks))) %>%
    gather_listoftbl(index = "roi", key = "temp") %>%
    select(-temp)

  return(chromdata_rois)
}

filter_rois2 <- function(chromdata, peaks, roi_width, x = "mz", cores = 1){

  # NTS 20200405: FIGURE OUT HOW TO LET USER SPECIFY x NAME!

  if(cores == 1){
    chromdata_rois <- chromdata %>%
      rowwise() %>%
      mutate(
        roi = list(mapply(
          # feed ROI rts, mzs, indices in parallel
          peaks_master$rt,
          peaks_master$mz,
          seq(nrow(peaks_master)),
          FUN = function(x, y, z){
            # should be the fastest-eliminating check order
            if(
              (x > rt - roi_width/2) &&
              (x < rt + roi_width/2) &&
              (y == mz)
            ){
              return(z)
            }else{
              #return(0)
            }
          }
        ))
      ) %>%
      # filter out rows not belonging to a ROI
      filter(roi != 0)
  }else{
    # if parallel, use pbapply
    # I think this is a bad idea; it parallelizes over ROIs, not over chromdata.
    clust <- makeCluster(cores, type="FORK")
    on.exit(stopCluster(clust)) # failsafe
    message("filtering for ROIs")
    chromdata_rois <- chromdata %>%
      rowwise() %>%
      filter(any(pbmapply(
        # feed ROI rts, mzs into filter in parallel
        peaks_master$rt,
        peaks_master$mz,
        FUN = function(x, y){
          (x > rt - roi_width/2) &&
            (x < rt + roi_width/2) &&
            (y == mz)
        }
      )))
    stopCluster(clust)
  }

  return(chromdata_rois)
}
