#' Cluster peaks
#' Takes a summary table of peaks with \code{rt_min}/\code{rt_max} columns,
#' clusters the peaks based on fraction overlap,
#' returns the same tibble with cluster index.
#'
#' @param peaks Tibble with \code{rt_min}/\code{rt_max}. If there is >1 row per peak, they should be grouped
#' @param thres_lap Fraction of RT overlap required for peaks to cluster (calculated as a fraction of the narrower peak)
#' @return Input tibble but with \code{clust_id} column added
#' @keywords cluster peaks
#' @export
#' @examples
#' peaks_clustered <- peaks %>%
#'   cluster_peaks(thres_lap = 0.8)
#'
cluster_peaks <- function(
  peaks,
  thres_lap
  ){

  # save the original groupings for output
  group_vars_in <- group_vars(peaks)

  # give peaks unique IDs
  # these are used to retrieve metadata later
  peaks$peak_id <- peaks %>% group_indices()

  # group peaks additionally by peak_id
  # doubles as a check for presence of rt_min, rt_max
  peaks <- peaks %>%
    group_by_at(c(group_vars(.), "peak_id", "rt_min", "rt_max"))

  # init pairwise similarity matrix
  pairs <- tibble(
    x = peaks %>% pull(peak_id) %>% max() %>% seq(),
    y = x
    ) %>%
    tidyr::expand(x, y) %>%
    filter(x < y)

  # calc pairwise RT overlaps
  message("calculating pairwise RT overlaps")
  peaks_summary <- peaks %>%
    summarize()
  rt_mins <- peaks_summary %>%
    pull(rt_min)
  rt_maxs <- peaks_summary %>%
    pull(rt_max)
  # filter pairs by RT tolerance
  pairs_lap <- pairs %>%
    mutate(
      # divisor for the overlap value (width of narrower peak)
      div = min((rt_maxs[x]-rt_mins[x]), (rt_maxs[y]-rt_mins[y])),
      lap = ifelse(
        # peaks overlap perfectly
        (rt_mins[x] == rt_mins[y]) & (rt_maxs[x] == rt_maxs[y]),
        # return 1
        1,
        ifelse(
          # peak x partially leads y
          ((rt_mins[x] > rt_mins[y]) & (rt_mins[x] < rt_maxs[y])),
          (rt_maxs[y] - rt_mins[x]) / div,
          ifelse(
            # peak y partially leads x
            ((rt_mins[y] > rt_mins[x]) & (rt_mins[y] < rt_maxs[x])),
            (rt_maxs[x] - rt_mins[y]) / div,
            # if no overlap, return 0
            0
          )
        )
      )
    ) %>%
    filter(lap > thres_lap)

  # use igraph to condense the spectrum pairs into clusters (disjoint sets)
  message("clustering...")
  graph <- graph_from_data_frame(pairs_lap, directed = FALSE)

  # just the input peaks summarized by cluster
  clusters <- tibble(peak_id = split(V(graph)$name %>% as.integer(), clusters(graph)$membership)) %>%
    mutate(clust_id = row_number()) %>%
    unnest(peak_id)

  peaks_clustered <- peaks %>%
    left_join(clusters, by = "peak_id") %>%
    # regroup like input
    group_by_at(group_vars_in) %>%
    # peak_id is internal; remove it
    select(-peak_id)

  return(peaks_clustered)
}
