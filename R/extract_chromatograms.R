#' Extract Chromatograms
#'
#' @param peaks Tibble with columns \code{scan} and \code{mz}/\code{wl}, and other identifiers
#' @param chromdata Full-spectrum chromatographic data
#' @param x Name of x-axis variable (\code{mz} or \code{wl}, etc.)
#' @return A tibble of XICs or single-wavelength chromatograms for the full peaks
#' centered on \code{scan} (ready for integration!) Peak ID columns like roi are preserved.
#' @keywords extract chromatogram
#' @export
#' @examples
#' xics_matched <- peaks_matched %>%
#' extract_chromatograms(chromdata)
#'
extract_chromatograms <- function(peaks, chromdata, x = "mz", cores = 1){

  chromdata <- chromdata %>%
    ungroup()
  x_peak <- paste(x, "peak", sep="_")

  # add rt_min, max columns to peak table
  peaks <- peaks %>%
    group_by(roi) %>%
    # so this symbol comes unambiguously from chromdata
    select(-intensity) %>%
    rename(
      rt = "rt_peak",
      mz = "mz_peak" #NTS: NEED TO WORK OUT X-SUBSTITUTION HERE AND IN BELOW FILTERS!!
      # maybe best to just rename to `x` and `x_peak` at top of the function
      ) %>%
    mutate(
      rt_min = chromdata %>%
        filter(mz == mz_peak) %>%
        # peak start criteria
        filter(
          (rt < rt_peak) & # comes before the peak
          (scan == lead(scan) - 1) & # next scan contains the same ion/wavelength
          (intensity < lead(intensity)) & # < the next value
          (
            (intensity <= lag(intensity)) | # <= the previous value (partially resolved peaks)
            is.na(lag(intensity)) | # or previous value DNE
            rt == min(rt) # or it's the very first point in the dataset
          )
        ) %>%
        # last value before the peak
        filter(rt == max(rt)) %>%
        pull(rt),
      rt_max = chromdata %>%
        filter(mz == mz_peak) %>%
        # peak end criteria
        filter(
          (rt > rt_peak) & # comes after the peak
          (scan == lag(scan) + 1) & # previous scan contains the same ion/wavelength
          (intensity < lag(intensity)) & # < the last value
          (
            (intensity <= lead(intensity)) | # <= the next value (partially resolved peaks)
            is.na(lead(intensity)) | # or previous value DNE
            rt == max(rt) # or it's the very last point in the dataset
          )
        ) %>%
        # first value after the peak
        filter(rt == min(rt)) %>%
        pull(rt)
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
        filter(
          (mz == row$mz_peak) &&
            (rt > row$rt_min) &&
            (rt < row$rt_max)
        )
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
