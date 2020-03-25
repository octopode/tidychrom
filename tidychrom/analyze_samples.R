library(xcms)
library(tidyverse)
library(RColorBrewer)
library(ggpubr)

## MAIN ANALYSIS BODY ##

for(dir_data in c(
  "/Users/jwinnikoff/Documents/MBARI/Lipids/GCMSData/cdf/20200213",
  "/Users/jwinnikoff/Documents/MBARI/Lipids/GCMSData/cdf/20200214",
  "/Users/jwinnikoff/Documents/MBARI/Lipids/GCMSData/cdf/20200215"
)){
  ## load data
  # starting with one dir of sample files and a standard mix
  #dir_data    <- "/Users/jwinnikoff/Documents/MBARI/Lipids/GCMSData/cdf/20200214"
  file_blank  <- "blank.cdf"
  file_stds   <- "Supel37_1_8.cdf"

  # read data file names
  cdfs <- list.files(path = dir_data, pattern = ".cdf", full.names = T)
  samps <- unlist(lapply(cdfs, function(x){sub("^([^.]*).*", "\\1", basename(x))}))
  blank_index <- match(file_blank, lapply(cdfs, basename))
  stds_index <- match(file_stds, lapply(cdfs, basename))

  # load raw MS data
  pd <- as(as.data.frame(samps), "AnnotatedDataFrame")
  raw_data <- readMSData(files = cdfs,
                         pdata = pd,
                         msLevel = 1,
                         mode = "onDisk"
  )

  # BPCs
  # use aggregationFun = sum" for TIC
  master <- chromatogram(raw_data, aggregationFun = "max")
  nscans <- length(rtime(master[,1]))

  # extract the traces of all master
  # can pull out like master_coords$JWL12 or master_coords$JWL12$intensity
  master_coords.m <- melt_coords(master)

  # detect peaks on the master of the standards to get general ROIs
  num_stds <- 37
  peaks_stds <- master_coords.m %>%
    filter(samp == "Supel37_1_8") %>%
    get_peaks(n = num_stds)

  # ggplot the detected master peaks
  master_coords.m %>%
    filter(samp == "Supel37_1_8") %>%
    ggplot() +
    geom_line(aes(x = rt/60, y = intensity, color = samp)) +
    geom_point(data = peaks_stds, aes(x = rt/60, y = intensity), shape = 25) +
    theme_pubr()

  # split raw data out into list of ROI chroms; takes a minute
  # non-essential; fine for plotting and such
  #roi_chroms <- divide_rois(raw_data, peaks_stds, lead = 4, lag = 4)

  # use spectrum distance (cosine) to match the peaks in each ROI
  # In principle, this could work if I gave R more memory. Would it be faster?
  #test <- peaks_stds %>%
  #  pull(rt) %>%
  #  lapply(., function(x){match_peaks(raw_data, std_index = stds_index, std_rt = x, cosine = 0.9)})

  # using loop instead of lapply here to prevent overflow
  # obviously the disk reads take awhile
  peaks_matched <- NULL
  rt_stds <- peaks_stds %>% arrange(rt) %>% pull(rt)
  for (i in seq(length(rt_stds))){
    peaks_matched <- match_peaks(raw_data, std_index = stds_index, std_rt = rt_stds[i], cosine = 0.9) %>%
      add_column(rt_std = rt_stds[i]) %>%
      bind_rows(peaks_matched, .)
    message(paste("(", i, "/", nrow(peaks_stds), ") matched peaks around ", rt_stds[i], " s", sep = "", collapse = ""))
  }

  save(peaks_matched, file = paste(dir_data, "/areas.save", sep = "", collapse = ""))
}
