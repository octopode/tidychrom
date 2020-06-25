#' Cluster Spectra
#' Takes a tibble of scan data and some clustering thresholds,
#' returns consensus spectra with cluster stats
#'
#' @param spectra Grouped tibble \code{mz}/\code{wl} and \code{intensity}. Should be grouped by scan.
#' @param x Name of x-axis variable (\code{mz} or \code{wl}, etc.)
#' @param thres_cos Min allowable cosine score for scans to cluster
#' @param thres_pks Min allowable number of matched peaks
#' @param thres_drt Max allowable RT difference for scans to cluster
#' @return Input tibble but with neighbor and cos columns added
#' @keywords cluster spectra
#' @export
#' @examples
#' spectra_clustered <- spectra %>%
#'   group_by(scan, file) %>%
#'   cluster_spectra()
#'
cluster_spectra <- function(
  spectra, bin, x = "mz",
  thres_cos = 0,
  thres_pks = 2,
  thres_drt = Inf,
  cores = 1
  ){

  # bin for comparison if requested
  if(!missing(bin)){
    spectra <- spectra %>%
      bin_spectra(bin_width = bin)
  }

  # give scans unique IDs
  # these are used to retrieve metadata later
  spectra$spec_id <- spectra %>% group_indices()

  # group spectra additionally by spec_id
  spectra <- spectra %>%
    group_by_at(c(group_vars(.), "spec_id"))

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

  # init pairwise similarity matrix
  pairs_all <- tibble(
    x = scans %>% length() %>% seq(),
    y = x
    ) %>%
    tidyr::expand(x, y) %>%
    filter(x < y)

  if(thres_drt != Inf){
    # calc pairwise delta RTs
    message("calculating pairwise delta RTs")
    rts <- spectra %>%
      distinct(rt) %>%
      pull(rt)
    # filter pairs by RT tolerance
    pairs_rt <- pairs_all %>%
      mutate(drt = abs(rts[x] - rts[y])) %>%
      filter(drt <= thres_drt)
  }else{
    pairs_rt <- pairs_all
  }

  if(thres_cos != 0){
    ## calc pairwise matched peaks
    ## this filter is run first to optimize parallel processing
    #message("calculating pairwise matched peaks")
    #peaks <- pairs_rt %>%
    #  pbapply(.,
    #          cl = clust,
    #          MARGIN = 1,
    #          FUN = function(row){
    #            inner_join(
    #              scans[[row[[1]]]],
    #              scans[[row[[2]]]]
    #            ) %>%
    #              nrow()
    #          }
    #  )
    #pairs_pks <- pairs_rt %>%
    #  mutate(pks = peaks) %>%
    #  filter(pks >= thres_pks)

    # calc pairwise cosines in parallel
    message("calculating pairwise cosine scores")
    #cosines <- pairs_pks %>%
    cosines <- pairs_rt %>%
      pbapply(.,
              cl = clust,
              MARGIN = 1,
              FUN = function(row){
                cosine_spectra(
                  scans[[row[[1]]]],
                  scans[[row[[2]]]],
                  thres_pks = thres_pks
                )
              }
      )
    # filter pairs by cosine threshold
    pairs_cos <- pairs_rt %>%
      mutate(cos = cosines) %>%
      filter(cos >= thres_cos)
  }else{
    pairs_cos <- pairs_rt
  }

  # use igraph to condense the spectrum pairs into clusters (disjoint sets)
  message("clustering...")
  graph <- graph_from_data_frame(pairs_cos, directed = FALSE)

  # just the input spectra summarized by cluster
  clusters <- tibble(spec_id = split(V(graph)$name %>% as.integer(), clusters(graph)$membership)) %>%
    mutate(clust_id = row_number()) %>%
    unnest(spec_id)

  scans_clustered <- spectra %>%
    summarize() %>%
    left_join(clusters, by = "spec_id") %>%
    # spec_id is internal; remove it
    select(-spec_id)

  return(scans_clustered)
}
