# old LDR-related functions

lol_fancy <- function(x, y, rsq = "0.95", min_dils = 3, force_origin = F){
  coords <- cbind(x, y) %>%
    as_tibble()
  for(x_max in sort(unique(x), decreasing = T)){
    coords <- coords %>%
      filter(x <= x_max)
    # if there are enough dilutions left
    n_dils <- length(unique(coords$x)) %>%
      as.integer()
    if(n_dils >= min_dils){
      # do a linear regression
      # with or without origin forced
      if(force_origin){
        linreg <- lm(y~x+0)
      }else{
        linreg <- lm(y~x)
      }
      # get the R^2
      rr <- linreg %>%
        summary() %>%
        .$r.squared
      if(rr >= rsq){
        return(list(
          y_max = max(coords$y),
          rr = rr,
          n_dils = n_dils
        ))
      }
    }else{
      # if too few dilutions, just
      return(list(
        y_max = NA,
        rr = NA,
        n_dils = NA
      ))
    }
  }
  # failsafe? In case someone sets min_dils = 0...
  return(list(
    y_max = NA,
    rr = NA,
    n_dils = NA
  ))
}

#' Calculate LDR Upper Limit (\code{intb_max})
#'
#' @param peak_group Tibble with columns \code{dil}, and \code{intb},
#' subsetted to a single compound, e.g. by \code{group_map()}
#' @param key compound identifier, passed as a one-element tibble, as by \code{group_map()}
#' @param rsq minimum R^2 at which signal is considered linear
#' @param min_dils minimum number of dilutions that must be present for a compound
#' @return tibble with max. linear peak area, number of data points, R^2, and compound id
#' @keywords ldr linear dynamic range
#' @export
#' @examples
#' use this function through \code{calc_ldrs}

# get intbmax and R^2 for a given compound
#calc_intb_max <- function(dil_col, intb_col, cpd_col, cpd, rsq){
calc_intb_max <- function(peak_group, key, rsq = 0.95, min_dils = 3){
  rsq <- unlist(rsq)
  rr <- 0
  dils <- unique(peak_group$dil) %>% sort()
  i <- length(dils)
  # if r^2 is too low, keep removing the upper data points
  # until it is acceptable
  while((rr < rsq) & (i >= min_dils - 1)){
    good_peaks <- peak_group %>%
      filter(dil %in% dils[1:i])
    rr <- good_peaks %>%
      lm(dil~intb, .) %>%
      summary() %>%
      .$r.squared
    i <- i-1
  }
  # return max linear peak area, number of dils, r-squared
  if(i + 1 >= min_dils){
    return(tibble(intb_max = c(max(good_peaks$intb)), npts = c(i+1), rr = c(rr), !!colnames(key) := key[[1]]))
  }else{
    # if there is not an adequate number of concs, return NAs
    return(tibble(intb_max = c(NA), npts = c(NA), rr = c(NA), !!colnames(key) := key[[1]]))
  }
}

#' Calculate LDRs for all peaks
#' A convenience wrapper for \code{calc_intb_max()}
#'
#' @param peaks_matched Tibble with columns \code{dil}, and \code{intb}
#' @param cpd_id String naming the column to be used as a unique compound ID.
#' Conservatively, \code{rt_std}. If you trust your dereplication, can be \code{id_db}.
#' @param rsq minimum R^2 at which signal is considered linear
#' @param min_dils minimum number of dilutions that must be present for a compound
#' @return tibble with max. linear peak area, number of data points, R^2, and compound id
#' @keywords ldr linear dynamic range
#' @export
#' @examples
#'

# convenience function to apply LDR values to the main tibble
calc_ldrs <- function(peaks_matched, cpd_id = "rt_std", rsq = 0.98, min_dils = 3){
  ldrs <- peaks_matched %>%
    group_by_(cpd_id) %>%
    group_map(~calc_intb_max(
      peak_group = .x,
      key = .y,
      rsq = rsq,
      min_dils = min_dils)
    ) %>% bind_rows()
  print(ldrs)
  peaks_matched <- peaks_matched %>%
    left_join(ldrs, cpd_id)
  return(peaks_matched)
}
