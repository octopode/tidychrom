## the heavy lifting
# take a raw data object, an index for the reference sample, target RT, and cos cutoff
# find either the tallest corresponding, or best-corresponding peaks wrt standard
# (biggest or best) return respective RTs, base m/z
match_peaks2 <- function(
  exp, std_index, std_rt,
  lead = 4, lag = 4, bin = F, mztol = 1,
  cosine = 0.9, method = "biggest", integrate = T,
  rtime_tbl # saves time
  ){
  if(missing(rtime_tbl)){
    # store all RT data from expt in a tbl (saves time)
    rtime_tbl <- exprtime2tbl(exp)
  }
  # get spectrum from the passed std_rt (center of std peak)
  # Not gonna trust the user to pass a perfectly on-peak RT, though.
  # Consequently, scotch the following and use the standard spectrum already identified within
  # the tibble group. - 20200330 JRW
  spec_std <- get_spec_rt(exp, std_index, std_rt, rtime_tbl)
  # binning probably not needed if all data collected with same instrument, method
  if(bin){
    spec_std <- spec_std %>%
      bin(breaks = seq(0.5, max(.@mz), mztol))
  }
  # get the basepeak to peak-pick all samples with
  # might wanna refactor this - 20200330 JRW
  basepeak <- get_basepeak(spec_std)
  # extract std-referenced, constrained BPCs from raw data
  xics <- chromatogram(
    exp,
    rt = c(std_rt - lag, std_rt + lead),
    mz = c(basepeak - mztol, basepeak + mztol),
    aggregationFun = "max"
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
      spectrum = ifelse(
        bin,
        get_spec_rt(exp, samp_index, rt, rtime_tbl) %>%
          # note bins broken on the 0.5
          bin(breaks = seq(0.5, max(.@mz), mztol)) %>%
          c(),
        get_spec_rt(exp, samp_index, rt) %>%
          c()
      ),
      mz = get_basepeak(spectrum)
    )
  # get the *real* master standard spectrum
  # see above about not trusting the caller
  spec_std <- xic_peaks %>%
    filter(index_samp == std_index) %>%
    group_by(1) %>%
    filter(intensity == max(intensity)) %>%
    pull(spectrum) %>%
    .[[1]]
  #print(spec_std) #TEST
  # then calculate each spectrum's distance from standard!
  xic_peaks <- xic_peaks %>%
      mutate(
        cos = compareSpectra(
        spec_std,
        spectrum,
        fun="dotproduct"
        )
    ) %>%
    # drop inadequate matches
    filter(cos >= cosine) %>%
    group_by(samp) %>%
    arrange(desc(cos))
  print(xic_peaks %>% select(-intensity, -spectrum)) #TEST
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

# sandbox 20200328
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
