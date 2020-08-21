#' Author's utility functions, to be loaded with \code{tidychrom} namespace

# using devtools and roxygen2, reload the tidychrom namespace
reload_tidychrom <- function(){
  library(roxygen2)
  library(devtools)
  document()
  detach("package:tidychrom")
  library(tidychrom)
}
