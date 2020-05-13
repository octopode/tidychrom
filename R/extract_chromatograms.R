#' Extract Chromatograms
#'
#' @param peaks Tibble with columns \code{scan} and \code{mz}/\code{wl}, and other identifiers
#' @param chromdata Full-spectrum chromatographic data
#' @param x Name of x-axis variable (\code{mz} or \code{wl}, etc.)
#' @param scans_gap Number of consecutive missing scans allowed in a peak
#' @return A tibble of XICs or single-wavelength chromatograms for the full peaks
#' centered on \code{scan} (ready for integration!) Peak ID columns like roi are preserved.
#' @keywords extract chromatogram
#' @export
#' @examples
#' xics_matched <- peaks_matched %>%
#' extract_chromatograms(chromdata)
#'
extract_chromatograms <- function(peaks, chromdata, x = "mz", scans_gap = 3, cores = 1){

  # differentiate XICs x3
  chromdata_diff <- chromdata %>%
    ungroup() %>%
    group_by(mz) %>%
    arrange(rt) %>%
    mutate(
      # first derivative (interpolated)
      d1I_dt1 =   ((intensity - lag(intensity)) / (rt - lag(rt)) + (lead(intensity) - intensity) / (lead(rt) - rt))/2,
      # second "
      d2I_dt2 = (lead(d1I_dt1) - lag(d1I_dt1)) / (lead(rt) - lag(rt)),
      # third "
      d3I_dt3 = (lead(d2I_dt2) - lag(d2I_dt2)) / (lead(rt) - lag(rt))
    )

  ## get minimum time between scans
  scan_time <- chromdata %>%
    select(rt) %>%
    mutate(delta = lead(rt) - rt) %>%
    filter(delta != 0) %>%
    pull(delta) %>%
    .[2:length(.)] %>%
    min()

  x_peak <- paste(x, "peak", sep="_")

  # add rt_min, max columns to peak table
  peaks <- peaks %>%
    group_by(roi) %>%
    # so this symbol comes unambiguously from chromdata
    select(-intensity) %>%
    rename(
      #scan = "scan_peak", # this breaks the function(??) So does removing it :/
      rt = "rt_peak",
      mz = "mz_peak" #NTS: NEED TO WORK OUT X-SUBSTITUTION HERE AND IN BELOW FILTERS!!
      # maybe best to just rename to `x` and `x_peak` at top of the function
      ) %>%
    mutate(
      rt_min = chromdata_diff %>%
        filter(
          (mz == mz_peak) & # extract signal ion
          (rt < rt_peak) # comes before the peak
          ) %>%
        # peak start criteria
        filter(lead(rt) <= rt + scans_gap*scan_time) %>%  # candidate point must be followed by a point <= n scans away
        filter(
          (
            # conditions for a fully resolved peak
            (d3I_dt3 >= 0) & (lead(d3I_dt3) < 0) # next point is a concavity maximum
          ) | # or
          (
            # for a partially resolved peak
            (lead(intensity) > intensity) & (lag(intensity) >= intensity) # this point is a crotch
          )
        ) %>%
        # last value before the peak
        pull(rt) %>%
        max() %>%
        # if there is no valid peak limit, nix the peak
        ifelse(. != -Inf, ., NA),
      rt_max = chromdata_diff %>%
        filter(
          (mz == mz_peak) &
          (rt > rt_peak) # comes after the peak
        ) %>%
        # peak end criteria
        filter(lag(rt) >= rt - scans_gap*scan_time) %>%  # candidate point must be followed by a point <= n scans away
        filter(
          (
            # conditions for a fully resolved peak
            (lag(d3I_dt3) > 0) & (d3I_dt3 <= 0) # previous point is a concavity maximum
          ) | # or
          (
            # for a partially resolved peak
            (lead(intensity) >= intensity) & (lag(intensity) > intensity) # this point is a crotch
          )
        ) %>%
        # first value after the peak
        pull(rt) %>%
        min() %>%
        # if there is no valid peak limit, nix the peak
        ifelse(. != Inf, ., NA),
    )

  # need to iterate over ROIs rather than filtering chromdata,
  # because theoretically, ROIs *could* overlap (even in mz)
  # best in parallel
  if(cores > 1){
    clust <- makeCluster(cores, type="FORK")
    on.exit(stopCluster(clust)) # important!
  }else{
    clust <- NULL
  }
  message("extracting windowed chromatograms")
  chromdata_out <- pblapply(
    cl = clust,
    seq(nrow(peaks)),
    function(i){
      row <- peaks[i,]
      chromdata_xtract <- chromdata %>%
        rowwise() %>%
        # so far as filter() is concerned,
        # NA (i.e. when peak DNE and rt_min or rt_max = NA) is the same as FALSE
        filter(
          (mz == row$mz_peak) &&
            (rt > row$rt_min) &&
            (rt < row$rt_max)
        )
      # add row to output tbl
      chromdata_xtract <- chromdata_xtract %>%
        bind_cols(
          row %>%
            # drop mapping columns
            select(
            -mz_peak,
            -rt_min,
            -rt_max,
            -scan,
            -rt_peak
            ) %>%
            slice(rep(1, each = nrow(chromdata_xtract)))
        )
      return(chromdata_xtract)
    }
  ) %>%
    # quick rbind the extracted chromdata
    do.call(rbind, .) %>%
    as_tibble()

  return(chromdata_out)
}
