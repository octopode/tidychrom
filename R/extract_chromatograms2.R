#' Extract Chromatograms 2
#' Same interface as \code{extract_chromatograms()}, but more robust to peak shape, saturation, etc.
#'
#' @param peaks Tibble with columns \code{scan} and \code{mz}/\code{wl}, and other identifiers
#' @param chromdata Full-spectrum chromatographic data
#' @param x Name of x-axis variable (\code{mz} or \code{wl}, etc.)
#' @param thres_hard Hard intensity threshold for the XIC

#' @return A tibble of XICs or single-wavelength chromatograms for the full peaks
#' centered on \code{scan} (ready for integration!) Peak ID columns like roi are preserved.
#' @keywords extract chromatogram
#' @export
#' @examples
#' xics_matched <- peaks_matched %>%
#' extract_chromatograms(chromdata)
#'
extract_chromatograms2 <- function(chromdata, x = "mz", thres_hard = 0, thres_snr = 100, cores = 1){

  ## get chromatographic data in line
  chromdata_complete <- chromdata %>%
    ungroup() %>%
    filter(intensity >= thres_hard) %>%
    # fill gaps in XICs!
    complete(rt, mz, fill = list(intensity = 0)) %>%
    group_by(mz) %>%
    arrange(rt)

  x_peak <- paste(x, "peak", sep="_")

  # add rt_min, max columns to peak table
  # peaks should be either grouped by ROI, or else rowwise
  peaks <- peaks %>%
    # so this symbol comes unambiguously from chromdata
    rename(
      rt = "rt_peak",
      mz = "mz_peak", #NTS: NEED TO WORK OUT X-SUBSTITUTION HERE AND IN BELOW FILTERS!!
      intensity = "intensity_peak"
      # maybe best to just rename to `x` and `x_peak` at top of the function
      ) %>%
    mutate(
      rt_min = chromdata_complete %>%
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
          ) &
          # at least 1/SNR below the peak "center"
          (intensity < intensity_peak / thres_snr) | # or
          # is the first data point of the XIC
          (rt == dplyr::first(rt))
        ) %>%
        # last value before the peak
        pull(rt) %>%
        max(),
      rt_max = chromdata_complete %>%
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
          ) &
          # at least 1/SNR below the peak "center"
          (intensity < intensity_peak / thres_snr) | # or
          # is the last data point of the XIC
          (rt == dplyr::last(rt))
        ) %>%
        # first value after the peak
        pull(rt) %>%
        min()
    ) %>%
    group_by(rt_min, rt_max) %>%
    summarize(mz_peak = fmode(mz_peak))

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
            #select(
            #-mz_peak,
            #-rt_min,
            #-rt_max,
            #-scan,
            #-rt_peak
            #) %>%
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
