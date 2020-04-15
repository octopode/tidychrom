# helper function extracts a tbl of feature, sample index, rtime
exprtime2tbl <- function(exp){
  rtime_tbl <- rtime(exp) %>%
    bind_rows() %>%
    gather(key = feature, value = rt) %>%
    # not sure why this is needed, but it is
    rowwise() %>%
    mutate(
      samp_num = strsplit(feature, "\\.") %>%
        .[[1]] %>% .[1] %>%
        substr(., 2, str_length(.)) %>%
        as.numeric()
    )
  return(rtime_tbl)
}
