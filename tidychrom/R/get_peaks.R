# take chromatogram coordinates
# return a tibble of either the tallest n peaks
# or all peaks with intensity >= int_min
get_peaks <- function(coords, n, int_min){
  if(missing(n) & missing(int_min)){
    stop("need to specify n or int_min!")
  }
  all_peaks_coords <- coords %>%
    filter((intensity > lag(intensity)) & (intensity > lead(intensity))) %>%
    arrange(desc(intensity))
  # top-n-peaks option takes priority
  if(!missing(n) & (nrow(all_peaks_coords) >= n)){
    all_peaks_coords <- all_peaks_coords %>%
      .[1:n,]
  }
  if(!missing(int_min)){
    all_peaks_coords <- all_peaks_coords %>%
      filter(intensity >= int_min)
  }
  return(all_peaks_coords)
}

# this version is supposed to take multiple samples in parallel
get_peaks2 <- function(coords, n, int_min){
  if(missing(n) & missing(int_min)){
    stop("need to specify n or int_min!")
  }
  all_peaks_coords <- coords %>%
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
  if(!missing(n) & (nrow(all_peaks_coords) >= n)){
    all_peaks_coords <- all_peaks_coords %>%
      # group by samp if column exists ("ask for permission")
      group_by(samp) %>%
      filter(row_number() <= n)
  }
  if(!missing(int_min)){
    all_peaks_coords <- all_peaks_coords %>%
      filter(intensity >= int_min)
  }
  return(all_peaks_coords %>% arrange(samp, rt))
}
