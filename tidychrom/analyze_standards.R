library(xcms)
library(tidyverse)
library(RColorBrewer)
library(ggpubr)

## CALIBRATION CURVES ##
## Read in serial dilutions of the standard to determine the LDR of area for each peak.
## These values can later be used for QC, to identify samples that are saturating.

# read data file names
dir_data    <- "/Users/jwinnikoff/Documents/MBARI/Lipids/GCMSData/cdf/20200212/stds"
cdfs <- list.files(path = dir_data, pattern = ".cdf", full.names = T)
samps <- unlist(lapply(cdfs, function(x){sub("^([^.]*).*", "\\1", basename(x))}))
file_stds   <- "Supel37_1_8.cdf"
stds_index <- match(file_stds, lapply(cdfs, basename))

# read data file names
cdfs <- list.files(path = dir_data, pattern = ".cdf", full.names = T)

# load raw MS data
pd <- as(as.data.frame(samps), "AnnotatedDataFrame")
raw_data <- readMSData(files = cdfs,
                       pdata = pd,
                       msLevel = 1,
                       mode = "onDisk"
)

# master chromatogram used to get ROIs
# use aggregationFun = "sum" for TIC
master <- chromatogram(raw_data, aggregationFun = "max")
nscans <- length(rtime(master[,1]))

# extract the traces of all master
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
  theme_pubr() +
  theme(legend.position = "right")

# match peaks by RT window and cosine score
# instead of just RT, since integration will be done on base peak
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

## IDENTIFY STANDARDS ##

peaks_matched <- annot_from_db(peaks_matched,
                               "/Users/jwinnikoff/Documents/MBARI/Lipids/GCMSData/db/supel37/")

# can visually inspect spectrum comparisons
peaks_master <- peaks_matched %>%
  filter(samp == "Supel37_1_8")
peaknum = 20
peaks_master %>%
  select(spectrum, spectrum_db) %>%
  .[peaknum,] %>%
  # put the spectrum and spectrum_db columns into the same column
  gather() %>% select(value) %>%
  plot_spectra() +
    ggtitle(paste("peak", as.character(peaknum), "\n",
                  peaks_master %>%
                    .[peaknum,] %>%
                    select(id_db)))

## DETERMINE LDRs ##

# named list to assign dilutions based on file naming convention
sname2dil <- c(
  Supel37_1_8     = 1/8,
  Supel37_1_8_2   = 1/8,
  Supel37_1_40    = 1/40,
  Supel37_1_40_2  = 1/40,
  Supel37_1_200   = 1/200,
  Supel37_1_200_2 = 1/200,
  Supel37_1_1k    = 1/1000,
  Supel37_1_1k_2  = 1/1000,
  Supel37_1_5k_2  = 1/5000,
  Supel37_1_5k    = 1/5000
)

# plot standard curves for each peak
peaks_matched %>%
  mutate(dil = sname2dil[samp], intb = unlist(intb)) %>%
  ggplot(aes(x = dil, y = intb)) +
    facet_wrap(facets = vars(rt_std), nrow = 6, ncol = 7) +
    geom_point() +
    geom_smooth(method = "lm") +
    stat_cor(size = 3, aes(label = ..rr.label..)) +
    theme_pubr() +
    scale_y_continuous(labels = function(x) format(x, scientific = TRUE)) +
    ggtitle("peak area standard curves by RT (s)") +
    xlab("dilution factor") +
    ylab("baseline-adjusted area")
