# take chromatogram coordinates
# return a tibble of either the tallest n peaks
# or all peaks with intensity >= int_min
get_peak_scans <- function(chromdata, n, int_min){
  if(missing(n) & missing(int_min)){
    stop("need to specify n or int_min!")
  }
  peaks <- chromdata %>%
    ungroup() %>%
    arrange(scan) %>%
    filter(
      # make sure it's within a contiguous ROI in case working with filtered (rt-windowed) data
      (scan == lag(scan) + 1) &&
        (scan == lead(scan) - 1) &&
        # and that it's a peak
        (intensity > lag(intensity)) &&
        (intensity > lead(intensity))
      ) %>%
    arrange(desc(intensity))
  # top-n-peaks option takes priority
  if(!missing(n) & (nrow(peaks) >= n)){
    peaks <- peaks %>%
      .[1:n,]
  }
  if(!missing(int_min)){
    peaks <- peaks %>%
      filter(intensity >= int_min)
  }
  return(peaks %>% arrange(rt))
}

# this version is supposed to take multiple samples in parallel
get_peaks2 <- function(chromdata, n, int_min){
  if(missing(n) & missing(int_min)){
    stop("need to specify n or int_min!")
  }
  peaks <- chromdata %>%
    # arrange by samp if column exists ("ask for permission")
    arrange(samp, rt) %>%
    filter(
      (intensity > lag(intensity)) &
        (intensity > lead(intensity)) &
        (samp == lag(samp)) &
        (samp == lead(samp))
        ) %>%
    arrange(desc(intensity))
  # top-n-peaks option takes priority
  if(!missing(n) & (nrow(peaks) >= n)){
    peaks <- peaks %>%
      # group by samp if column exists ("ask for permission")
      group_by(samp) %>%
      filter(row_number() <= n)
  }
  if(!missing(int_min)){
    peaks <- peaks %>%
      filter(intensity >= int_min)
  }
  return(peaks %>% arrange(samp, rt))
}
