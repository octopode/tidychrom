#' Author's utility functions, to be loaded with \code{tidychrom} namespace

# using devtools and roxygen2, reload the tidychrom namespace
reload_tidychrom <- function(){
  document()
  detach("package:tidychrom")
  library(tidychrom)
}
