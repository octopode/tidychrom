## the heavy lifting
# take a raw data object, an index for the reference sample, target RT, and cos cutoff
# find either the tallest corresponding, or best-corresponding peaks wrt standard
# (biggest or best) return respective RTs, base m/z
match_peaks <- function(
  exp, std_index, std_rt,
  lead = 4, lag = 4, mztol = 1,
  cosine = 0.9, method = "biggest", integrate = T){
  # get spectrum from the passed std_rt (center of std peak)
  spec_std <- get_spec_rt(exp, std_index, std_rt) %>%
  # binning probably not needed if all data collected with
  # same instrument, method
  bin(binSize = mztol)
  # get the basepeak to peak-pick all samples with
  basepeak <- get_basepeak(spec_std)
  # extract std-referenced, constrained BPCs from raw data
  xics <- chromatogram(
    exp,
    rt = c(std_rt - lag, std_rt + lead),
    mz = c(basepeak - mztol, basepeak + mztol)
  )
  # melt them down (for use now and later)
  xics_coords <- xics %>%
    melt_coords() %>%
    group_by(samp)
  # and get the peaks and spectra from them
  xic_peaks <- xics_coords %>%
    group_map(~ get_peaks(.x, n = 5)) %>%
    set_names(exp@phenoData$samp) %>%
    gather_listoftbl() %>%
    # then retrieve the peak spectra
    rowwise() %>%
    mutate(
      samp_index = match(samp, exp@phenoData$samp),
      # need to encapsulate spectrum in a vector to put in column
      spectrum = c(get_spec_rt(exp, samp_index, rt) %>% bin(binSize = mztol)),
      mz = get_basepeak(spectrum),
      # calculate each peak's distance from standard!
      cos = compareSpectra(spec_std, spectrum, fun="dotproduct")
    ) %>%
    # drop inadequate matches
    filter(cos >= cosine) %>%
    group_by(samp) %>%
    arrange(desc(cos))
  if(method == "best"){
    xic_peaks <- xic_peaks %>%
      filter(cos == max(cos))
  }
  if(method == "biggest"){
    xic_peaks <- xic_peaks %>%
      filter(intensity == max(intensity))
  }
  # finally, use the tabular trace data to integrate each peak
  if(integrate){
    xic_peaks <- integrate_coords(xic_peaks, xics_coords)
  }
  return(xic_peaks)
}

match_peaks2 <- function(
  exp, std_index, std_rt,
  lead = 4, lag = 4, mztol = 1,
  cosine = 0.9, method = "biggest", integrate = T){
  # get spectrum from the passed std_rt (center of std peak)
  spec_std <- get_spec_rt(exp, std_index, std_rt) %>%
    # binning probably not needed if all data collected with
    # same instrument, method
    bin(binSize = mztol)
  # get the basepeak to peak-pick all samples with
  basepeak <- get_basepeak(spec_std)
  # extract std-referenced, constrained BPCs from raw data
  xics <- chromatogram(
    exp,
    rt = c(std_rt - lag, std_rt + lead),
    mz = c(basepeak - mztol, basepeak + mztol)
  )
  # melt them down (for use now and later)
  xics_coords <- xics %>%
    melt_coords() %>%
    group_by(samp)
  # and get the peaks and spectra from them
  xic_peaks <- xics_coords %>%
    group_map(~ get_peaks(.x, n = 5)) %>%
    set_names(exp@phenoData$samp) %>%
    gather_listoftbl() %>%
    # then retrieve the peak spectra
    rowwise() %>%
    mutate(
      samp_index = match(samp, exp@phenoData$samp),
      # need to encapsulate spectrum in a vector to put in column
      spectrum = c(get_spec_rt(exp, samp_index, rt) %>% bin(binSize = mztol)),
      mz = get_basepeak(spectrum),
      # calculate each peak's distance from standard!
      cos = compareSpectra(spec_std, spectrum, fun="dotproduct")
    ) %>%
    # drop inadequate matches
    filter(cos >= cosine) %>%
    group_by(samp) %>%
    arrange(desc(cos))
  if(method == "best"){
    xic_peaks <- xic_peaks %>%
      filter(cos == max(cos))
  }
  if(method == "biggest"){
    xic_peaks <- xic_peaks %>%
      filter(intensity == max(intensity))
  }
  # finally, use the tabular trace data to integrate each peak
  if(integrate){
    xic_peaks <- integrate_coords(xic_peaks, xics_coords)
  }
  return(xic_peaks)
}
