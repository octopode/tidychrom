#' Cluster Spectra
#' Takes a grouped tibble of scan data, clusters scans based on a cosine threshold,
#' returns the same tibble with a cluster index.
#'
#' @param spectra Grouped tibble \code{mz}/\code{wl} and \code{intensity}. Should be grouped by scan.
#' @param x Name of x-axis variable (\code{mz} or \code{wl}, etc.)
#' @param thres_cos Min allowable cosine score for scans to cluster
#' @param thres_pks Min allowable number of matched peaks
#' @return Input tibble with \code{clust_id} column added
#' @keywords cluster spectra cosine
#' @export
#' @examples
#' spectra_clustered <- spectra %>%
#'   group_by(scan, file) %>%
#'   cluster_spectra(thres_cos = 0.95)
#'
cluster_spectra <- function(
  spectra,
  thres_cos,
  thres_pks = 2,
  x = "mz",
  bin,
  cores = 1
  ){

  # save the original groupings for output
  group_vars_in <- group_vars(spectra)

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
  pairs <- tibble(
    x = scans %>% length() %>% seq(),
    y = x
    ) %>%
    tidyr::expand(x, y) %>%
    filter(x < y)

  # calc pairwise cosines in parallel
  message("calculating pairwise cosine scores")
  cosines <- pairs %>%
    pbapply(.,
            cl = clust,
            MARGIN = 1,
            FUN = function(pair){
              cosine_spectra(
                scans[[pair[[1]]]],
                scans[[pair[[2]]]],
                thres_pks = thres_pks
              )
            }
    )
  # filter pairs by cosine threshold
  pairs_cos <- pairs %>%
    mutate(cos = cosines) %>%
    filter(cos >= thres_cos)

  # use igraph to condense the spectrum pairs into clusters (disjoint sets)
  message("clustering...")
  graph <- graph_from_data_frame(pairs_cos, directed = FALSE)

  # just the input spectra summarized by cluster
  clusters <- tibble(spec_id = split(V(graph)$name %>% as.integer(), clusters(graph)$membership)) %>%
    mutate(clust_id = row_number()) %>%
    unnest(spec_id)

  spectra_clustered <- spectra %>%
    left_join(clusters, by = "spec_id") %>%
    # regroup like input
    group_by_at(group_vars_in) %>%
    # spec_id is internal; remove it
    select(-spec_id)

  return(spectra_clustered)
}
