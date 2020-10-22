#' Cluster scans
#' Takes a summary table of scans with \code{rt_min}/\code{rt_max} columns,
#' clusters the scans based on RT proximity,
#' returns the same tibble with cluster index.
#'
#' @param scans Tibble with \code{rt_min}/\code{rt_max}. If there is >1 row per scan, they should be grouped.
#' @param thres_drt Maximum difference in RT allowed for scans to cluster
#' @return Input tibble but with \code{clust_id} column added
#' @keywords cluster scans
#' @export
#' @examples
#' scans_clustered <- scans %>%
#'   cluster_scans(thres_lap = 0.8)
#'
cluster_scans <- function(
  scans,
  thres_drt
  ){

  # save the original groupings for output
  group_vars_in <- group_vars(scans)

  # give scans unique IDs
  # these are used to retrieve metadata later
  scans$scan_id <- scans %>% group_indices()

  # group scans additionally by scan_id
  scans <- scans %>%
    group_by_at(c(group_vars(.), "scan_id", "rt"))

  # init pairwise similarity matrix
  pairs <- tibble(
    x = scans %>% pull(scan_id) %>% max() %>% seq(),
    y = x
    ) %>%
    tidyr::expand(x, y) %>%
    filter(x < y)

  # calc pairwise delta RTs
  message("calculating pairwise delta RTs")
  rts <- scans %>%
    summarize() %>%
    pull(rt)
  # filter pairs by delta RT threshold
  pairs_drt <- pairs %>%
    mutate(drt = abs(rts[x] - rts[y])) %>%
    filter(drt <= thres_drt)

  # use igraph to condense the spectrum pairs into clusters (disjoint sets)
  message("clustering...")
  graph <- graph_from_data_frame(pairs_drt, directed = FALSE)

  # just the input scans summarized by cluster
  clusters <- tibble(scan_id = split(V(graph)$name %>% as.integer(), clusters(graph)$membership)) %>%
    mutate(clust_id = row_number()) %>%
    unnest(scan_id)

  scans_clustered_nas <- scans %>%
    left_join(clusters, by = "scan_id") %>%
    ungroup()

  # give the singletons cluster IDs
  scans_clustered <- bind_rows(
    scans_clustered_nas %>% filter(!is.na(clust_id)),
    scans_clustered_nas %>% filter(is.na(clust_id)) %>%
      mutate(clust_id = (scan_id %>% as.factor() %>% as.integer()) + (scans_clustered_nas %>% pull(clust_id) %>% max(c(., 0), na.rm=TRUE)))
  ) %>%
  # regroup like input
  group_by_at(group_vars_in) %>%
  # scan_id is internal; remove it
  select(-scan_id)

  return(scans_clustered)
}
