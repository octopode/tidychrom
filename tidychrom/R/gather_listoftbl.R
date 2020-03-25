# utility func that takes a named list of tibbles,
# rbinds them, and adds a column corresponding to element names (samp)
gather_listoftbl <- function(data){
  data.m <- NULL
  # going off the rails and using a for loop, but it's pretty fast
  for(name in names(data)){
    data.m <- data[[name]] %>%
      add_column(., `samp` = rep(name, nrow(data[[name]]))) %>%
      bind_rows(data.m, .)
  }
  return(data.m)
}
