#' Upper Limit of Linearity (LOL)
#'
#' Take lists of x and y, return the upper limit of linearity.
#' @param x independent variable, e.g. concentration or dilution
#' @param y dependent variable, e.g. into or intb
#' @param rsq minimum R^2 at which signal is considered linear
#' @param min_dils minimum number of unique x values (dilutions) that must be present for a compound
#' @param force_origin Force regression through the origin? Set T for pre-blanked data, otherwise F.
#' @return a named list:
#' \code{$y_max}: maximum \code{y} for which R^2 >= \code{rsq}
#' \code{$n_dils}: number of dilutions used in the linreg
#' \code{$rr}: R^2 of the linreg
#' # NTS: can put slope and intercept in here too, if need be
#' @keywords upper limit linearity
#' @export
#' @examples
#' # use inside a dplyr::summarise() call
#' ldrs_roi <- areas_all %>%
#'   group_by(roi) %>%
#'   summarise(intb_max = lol(intb))
#'

# simple version returns only y_max
lol <- function(x, y, rsq = "0.95", min_dils = 3, force_origin = F){
  coords <- cbind(x, y) %>%
    as_tibble()
  for(x_max in sort(unique(x), decreasing = T)){
    coords <- coords %>%
      filter(x <= x_max)
    # if there are enough dilutions left
    n_dils <- length(unique(coords$x)) %>%
      as.integer()
    if(n_dils >= min_dils){
      # get R^2
      rr <- r_squared(x, y, force_origin = force_origin)
      if(rr >= rsq){
        return(max(coords$y))
      }
    }else{
      # if too few dilutions, just
      return(NA)
    }
  }
  # failsafe? In case someone sets min_dils = 0...
  return(NA)
}

# helper function, does what it looks like
r_squared <- function(x, y, force_origin = F){
  if(force_origin){
    linreg <- lm(y~x+0)
  }else{
    linreg <- lm(y~x)
  }
  # get the R^2
  rr <- linreg %>%
    summary() %>%
    .$r.squared
  return(rr)
}
