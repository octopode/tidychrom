#' Extract Chromatograms
#'
#' @param peaks Tibble with columns \code{scan} and \code{mz}/\code{wl}, and other identifiers
#' @param chromdata Full-spectrum chromatographic data
#' @param x Name of x-axis variable (\code{mz} or \code{wl}, etc.)
#' @param threshold Hard intensity threshold for the XIC
#' @param scans_gap Number of consecutive missing scans allowed in a peak
#' @return A tibble of XICs or single-wavelength chromatograms for the full peaks
#' centered on \code{scan} (ready for integration!) Peak ID columns like roi are preserved.
#' @keywords extract chromatogram
#' @export
#' @examples
#' xics_matched <- peaks_matched %>%
#' extract_chromatograms(chromdata)
#'
extract_chromatograms <- function(peaks, chromdata, x = "mz", scans_gap = 3, threshold = 0, cores = 1){

  ## differentiate XICs x3
  chromdata <- chromdata %>%
    ungroup() %>%
    filter(intensity >= threshold) %>%
    group_by(mz) %>%
    arrange(rt)# %>%
  #  mutate(
  #    # first derivative (interpolated)
  #    d1I_dt1 =   ((intensity - lag(intensity)) / (rt - lag(rt)) + (lead(intensity) - intensity) / (lead(rt) - rt))/2,
  #    # second "
  #    d2I_dt2 = (lead(d1I_dt1) - lag(d1I_dt1)) / (lead(rt) - lag(rt)),
  #    # third "
  #    d3I_dt3 = (lead(d2I_dt2) - lag(d2I_dt2)) / (lead(rt) - lag(rt))
  #  )

  ## get minimum time between scans
  scan_time <- chromdata %>%
    ungroup() %>%
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
      rt_min = chromdata %>%
        filter(
          (mz == mz_peak) & # extract signal ion
          (rt < rt_peak) # comes before the peak
        ) %>%
        # peak start criteria
        filter(
          # is a local minimum
          (
            (lead(intensity) > intensity) &
            (lag(intensity) >= intensity)
          ) | # or
          # follows a gap
          (lag(rt) <= (rt - scans_gap*scan_time)) | # or
          # is the first data point of the XIC
          (rt == dplyr::first(rt))
        ) %>%
        # last value before the peak
        pull(rt) %>%
        max(),# %>%
        # if there is no valid peak limit, nix the peak
        #ifelse(. != -Inf, ., NA),
      rt_max = chromdata %>%
        filter(
          (mz == mz_peak) &
          (rt > rt_peak) # comes after the peak
        ) %>%
        # peak end criteria
        filter(
          # is a local minimum
          (
            (lead(intensity) >= intensity) &
            (lag(intensity) > intensity)
          ) | # or
          # precedes a gap
          (lead(rt) >= (rt + scans_gap*scan_time)) | # or
          # is the last data point of the XIC
          (rt == dplyr::last(rt))
        ) %>%
        # first value after the peak
        pull(rt) %>%
        min()# %>%
        # if there is no valid peak limit, nix the peak
        #ifelse(. != Inf, ., NA),
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
