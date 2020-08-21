## STEP 1: BLANKING

# Preprocessing script to read in blank chromatographic runs, average them scanwise,
# then subtract the average from a set of experimental runs

# file selection criteria
# location of all the session subdirectories (or symlinks)
base_path = "~/Documents/MBARI/Lipids/GCMSData/cdf/"
# filename hash pattern to identify blanks
pattern_blank = "blank"

# in each output file,
# replace this pattern
# also used to identify input patterns. Dot is required.
suffix_in  = "\\.cdf"
# with this one:
suffix_out = "_blanked.mzxml"

# hard threshold for denoising chromatographic data
# full scans are kept if the base peak meets threshold
# setting to zero disables thresholding
threshold_bp <- 0
threshold <- 0

cores <- detectCores()
if(cores > 1){
  clust <- makeCluster(cores, type="FORK")
  on.exit(stopCluster(clust)) # important!
}else{
  clust <- NULL
}

for(subdir in c(
  # session subdirectories
   #"20200212"
   #"20200213",
   #"20200214",
   #"20200215",
   #"20200608",
   #"20200609"
   #"20200610",
   #"20200611",
   #"20200612",
   "20200613"
)){
  dir_data = paste(base_path,subdir,sep="")
  # filter for file type
  cdfs <- list.files(path = dir_data,
                     pattern = suffix_in %>%
                       strsplit(., "\\.") %>%
                       unlist %>% .[[2]] %>%
                       paste(".", ., sep=""),
                     full.names = T)
  # blank runs
  files_blanks <- cdfs[which(str_detect(cdfs, pattern_blank))]
  # experimental runs
  files_exptal <- setdiff(cdfs, files_blanks)

  message(paste("preparing blank(s) in", basename(dir_data)))
  print(files_blanks)
  # load the blank runs
  chromdata_blanks <- files_blanks %>%
    # first into a list of tibbles
    pblapply(., read_tidymass, cl = clust) %>%
    # with m/z binned to integer
    pblapply(., function(x){bin_spectra(x, bin_width = 1)}, cl = clust) %>%
    # then condense to a single tibble
    setNames(basename(files_blanks)) %>%
    gather_listoftbl(key = "file", index = "blank_num")

  # average them scanwise
  chromdata_blank_avg <- chromdata_blanks %>%
    group_by(scan, mz) %>%
    summarise(
      rt = mean(rt),
      intensity = mean(intensity)
    )

  # plot BPC of the blanks and the average
  gg_bpcs <- chromdata_blank_avg %>%
    mutate(file = "average") %>%
    bind_rows(chromdata_blanks) %>%
    mutate(file = fct_rev(file)) %>%
    group_by(file, scan) %>%
    # gets base peak intensity
    filter(intensity == max(intensity)) %>%
    ggplot() +
      geom_line(aes(x = rt, y = intensity, color = file)) +
      scale_color_brewer(palette = "Dark2") +
      theme_pubr() +
      theme(legend.position = "right") +
      ggtitle(paste("BPCs: blank(s) in session", basename(dir_data)))

  print(gg_bpcs)

  # then, looping thru experimental files,
  for(file_in in files_exptal){
    message(paste("subtracting blank from", basename(file_in)))
    chromdata_exptal <- file_in %>%
      read_tidymass() %>%
      bin_spectra(bin_width = 1) %>%
      # subtract the blank scanwise
      left_join(chromdata_blank_avg, by = c("scan", "mz")) %>%
      replace_na(list(intensity.y = 0)) %>%
      group_by(scan, mz) %>%
      mutate(
        rt = rt.x,
        intensity = max(intensity.x - intensity.y, 0)
      ) %>%
      select(scan, rt, mz, intensity) %>%
      group_by(scan) %>%
      filter(max(intensity) > threshold_bp) %>%
      ungroup() %>%
      filter(intensity > threshold)

    # and save
    file_out <- str_replace(file_in, suffix_in, suffix_out)
    chromdata_exptal %>% write_tidymass(file_out)
  }
}
