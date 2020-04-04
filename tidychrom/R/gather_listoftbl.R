#' Gather List of Tibbles (or other dataframe-like type)
#' Names and index of the list are saved to columns
#'
#' @param list_in A named list of tibbles or another dataframe-like type. Must all have matching columns.
#' @param key Name for the tbl column that original list names will end up in.
#' @param index Name for the tbl column that original list indices will end up in.
#' @param cores Number of cores to use for parallel ops on the input list. Only beneficial in huge tables (1E7s of rows.)
#' @keywords gather list tibble named
#' @export
#' @examples
#' ions <- ions <- mzr %>%
#'    spectra() %>%
#'    setNames(., header(mzr)$retentionTime) %>%
#'    gather_listoftbl(., key = "rt", index = "scan")
#'

# 20200404: new, faster approach using arrays and do.call
gather_listoftbl <- function(list_in, key = "key", index = "index", cores = 1){

  # presently set up for a named list only
  keys <- names(list_in)

  if(cores == 1){
    tbl_out <- lapply(seq(length(keys)), function(x){
      list_in %>%
        .[[x]] %>%
        as_tibble() %>%
        mutate(
          key = keys[x],
          index = x
        )
    })
  }else{
    # if parallel, use pbapply
    # though this doesn't seem to buy much time
    clust <- makeCluster(cores, type="FORK")
    on.exit(stopCluster(clust)) # failsafe
    tbl_out <- pblapply(cl = clust, seq(length(keys)), function(x){
      list_in %>%
        .[[x]] %>%
        as_tibble() %>%
        mutate(
          key = keys[x],
          index = x
        )
    })
    stopCluster(clust)
  }

  # fastest per https://stackoverflow.com/questions/28339514/convert-list-of-matrices-of-the-same-order-to-a-array
  tbl_out <- tbl_out %>%
    do.call(rbind, .) %>%
    as_tibble() %>%
    rename(
      key = key,
      index = index
    )

  return(tbl_out)
}

## utility function from https://stackoverflow.com/questions/3318333/split-a-vector-into-chunks-in-r
#chunk <- function(x,n) split(x, cut(seq_along(x), n, labels = FALSE))
#
#gather_listoftbl <- function(data, key = "samp", parallel = F){
#  if(!parallel){
#    # run the core routine
#    data.m <- NULL
#    # going off the rails and using a for loop, but it's pretty fast
#    for(name in names(data)){
#      ## if it's a list of some other df-like type, convert it
#      #if(!is_tibble(data[[name]])){
#      #  data[[name]] %>% data[[name]] %>%
#      #    as_tibble()
#      #}
#      data.m <- data[[name]] %>%
#        mutate(., !!key := rep(name, nrow(data[[name]]))) %>%
#        bind_rows(data.m, .)
#    }
#    return(data.m)
#  }else{
#    #build a cluster and run the routine recursively
#    # Calculate the number of cores
#    num_cores <- detectCores()
#    # split up the passed named list
#    data_map <- data %>%
#      chunk(num_cores)
#    # Initiate cluster with current environment
#    clust <- makeCluster(num_cores, type="FORK")
#    message("concatenating peaktable")
#    out_map <- pblapply(cl = clust, data_map, function(x){gather_listoftbl(x, key = key, parallel = F)})
#    stopCluster(clust)
#    return(bind_rows(out_map))
#  }
#}
