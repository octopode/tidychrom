#' Calculate Cosine Scores Pairwise
#' Takes a grouped tibble of scan data, calculates pairwise cosine similarities
#'
#' @param spectra Grouped tibble \code{mz}/\code{wl} and \code{intensity}. Should be grouped by scan.
#' @param x Name of x-axis variable (\code{mz} or \code{wl}, etc.)
#' @param thres_cos Min allowable cosine score for scans to cluster
#' @param thres_rt Max allowable RT difference for scans to cluster
#' @return Input tibble but with neighbor and cos columns added
#' @keywords cluster spectra
#' @export
#' @examples
#' spectra_clustered <- spectra %>%
#'   group_by(scan, file) %>%
#'   cluster_spectra()
#'
cosine_pairwise <- function(spectra, x = "mz", thres_cos = 0, thres_rt = Inf, cores = 1, bin){

  # bin for comparison if requested
  if(missing(bin) == F){
    spectra <- spectra %>%
      bin_spectra(bin_width = bin)
  }

  # give scans unique IDs
  spectra$spec_id <- spectra %>% group_indices()

  # set up parallel
  if(cores > 1){
    clust <- makeCluster(cores, type="FORK")
    on.exit(stopCluster(clust)) # important!
  }else{
    clust <- NULL
  }

  # split scans into a list of tbls for this op
  scans <- spectra %>%
    # ensures index matches spec_id
    arrange(spec_id) %>%
    group_split()

  # make pairwise similarity matrix
  spec_ids <- scans %>% length() %>% seq()
  pairs <- bind_cols(x = spec_ids, y = spec_ids) %>%
    complete(x, y) %>%
    # no diagonal
    filter(x != y)

  # calc all pairwise cosines in parallel
  # for some unknown reason, pbmapply has no parallel functionality
  pairs$cos <- pairs %>% pbapply(.,
    cl = clust,
    MARGIN = 1,
    FUN = function(row){
      cosine_spectra(
        scans[[row[[1]]]],
        scans[[row[[2]]]]
      )
    }
  )

  pairs_clust <- pairs %>%
    filter(cos >= thres_cos) %>%
    dplyr::rename(
      spec_id.x = x,
      spec_id.y = y
    ) %>%
    left_join(spectra %>% dplyr::rename(spec_id.x = spec_id), by = "spec_id.x") %>%
    left_join(spectra %>% dplyr::rename(spec_id.y = spec_id), by = "spec_id.y")

  return(pairs_clust)
}
