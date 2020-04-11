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
