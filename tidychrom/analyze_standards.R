library(tidychrom)
library(RColorBrewer)

## Step 1: Blanking
## STEP 2: STANDARD CALIBRATION CURVES

## Read in serial dilutions of the standard to determine the LDR of area for each peak.
## The resultant LDR limits are important for Step 3: building a standard spectrum library,
## and Step 4: Experimental Peak Matching, where they are used for QC,
## to identify samples that are saturating.

# USER PARAMETERS
# location of blanked data files
dir_data    <- "/Users/jwinnikoff/Documents/MBARI/Lipids/GCMSData/cdf/20200212/stds"
# search patter for blanked data files
mzxmls <- list.files(path = dir_data, pattern = "blanked.mzxml", full.names = T)

# resolution of the mass analyzer
bin_width_mz <- 1
# index in mzxmls of the "master run" used to identify ROIs for the standards
index_master <- 10
# number of peaks to identify in the standard mix
n_stds <- 35
# width in sec of regions of interest (ROIs) used for matching standard peaks
roi_width <- 4
# minimum acceptable cosine similarity to consider two spectra matching
cos_min <- 0.9

# load raw data for the master
chromdata_master <- mzxmls[index_master] %>%
  read_tidymass() %>%
  bin_spectra(bin_width = bin_width_mz)
# then get the BPC
bpc_master <- chromdata_master %>%
  group_by(scan) %>%
  filter(intensity == max(intensity))
# then get the biggest n peaks in the BPC
peaks_master <- bpc_master %>%
  # filter for peaks
  ungroup() %>%
  filter((intensity > lag(intensity)) & (intensity > lead(intensity))) %>%
  # get the 35 tallest
  arrange(desc(intensity)) %>%
  .[1:n_stds,] %>%
  # and put them in order again
  arrange(rt) %>%
  # peak matching requires a roi column
  rownames_to_column(var = "roi") %>%
  mutate(roi = as.integer(roi))

# plot the master BPC with detected peaks: a visual check for n_stds
bpc_master %>%
  ggplot() +
    geom_line(aes(x = rt, y = intensity)) +
    geom_point(data = peaks_master, aes(x = rt, y = intensity), color = "red", shape = 25) +
    theme_pubr() +
    theme(legend.position = "none") +
    ggtitle(paste("BPC:", basename(mzxmls[index_master])))

# get peak spectra from the master sample, for peak matching
spectra_master <- chromdata_master %>%
  right_join(
    peaks_master,
    by = c("scan", "rt")) %>%
  rename(
    "mz.x" = "mz",
    "intensity.x" = "intensity"
  ) %>%
  select(-mz.y, -intensity.y) %>%
  ungroup() %>%
  arrange(scan, rt, mz)

# match and integrate peaks from standard dilutions
# this is done iteratively in a loop to avoid memory overflow
# (raw data from a single file can exceed 50 MB)

areas_all <- NULL
for(file in mzxmls){
  message(paste("(", match(file, mzxmls), "/", length(mzxmls), ")", sep="", collapese = ""))
  # load a run
  message(paste("loading", basename(file)))
  chromdata <- file %>%
    read_tidymass() %>%
    bin_spectra(bin_width = bin_width_mz)
  # filter to ROIs
  # multithreading this is a 2.5x speed boost
  chromdata_rois <- chromdata %>%
    filter_rois(peaks_master, roi_width, cores = detectCores())
  # identify peaks in ROIs
  message("detecting peaks of interest")
  peaks_candidate <- chromdata_rois %>%
    ungroup() %>%
    arrange(roi, scan) %>%
    filter(
      # cannot be edge of an ROI
      (roi == lag(roi)) &
        (roi == lead(roi)) &
        # and has to be a peak
        (intensity > lag(intensity)) &
        (intensity > lead(intensity))
      )

  # get spectra from candidate peaks
  message("extracting candidate spectra")
  spectra_candidate <- chromdata %>%
    right_join(peaks_candidate, by = c("scan", "rt")) %>%
    rename(
      "mz.x" = "mz",
      "intensity.x" = "intensity"
      ) %>%
    select(-mz.y, -intensity.y)

  # get the candidate scan numbers that best match the master
  message("matching spectra")
  peaks_matched <- spectra_candidate %>%
    cosine_by_roi_scan(spectra_master) %>%
    # get only the one best match for each ROI
    group_by(roi) %>%
    filter(
      cos > cos_min,
      cos == max(cos)
      ) %>%
    # and assign the master base peak as mz
    left_join(
      peaks_master %>%
        select(
          -scan,
          -rt
        ),
      by = "roi"
    )

  # get the XICs for those peaks (extracting on the master base peak)
  # console message is inside the function
  xics_matched <- peaks_matched %>%
    # parallelizing speeds it up
    extract_chromatograms(chromdata, cores = detectCores())

  # integrate those XICs
  message("integrating XICs")
  peak_areas <- xics_matched %>%
    group_by(roi) %>%
    auc() %>%
    # and append column identifying the file
    mutate(file = basename(file))

  areas_all <- areas_all %>%
    bind_rows(peak_areas)
}


## DETERMINE LDRs ##

# named list to assign dilutions based on file naming convention
filename2dil <- c(
  Supel37_1_1k_2_blanked.mzxml  = 1/1000,
  Supel37_1_1k_blanked.mzxml    = 1/1000,
  Supel37_1_200_2_blanked.mzxml = 1/200,
  Supel37_1_200_blanked.mzxml   = 1/200,
  Supel37_1_40_2_blanked.mzxml  = 1/40,
  Supel37_1_40_blanked.mzxml    = 1/40,
  Supel37_1_5k_2_blanked.mzxml  = 1/5000,
  Supel37_1_5k_blanked.mzxml    = 1/5000,
  Supel37_1_8_2_blanked.mzxml   = 1/8,
  Supel37_1_8_blanked.mzxml     = 1/8
)

areas_all <- areas_all %>%
  mutate(dil = filename2dil[file], intb = unlist(intb))

ldrs <- calc_ldrs(peaks_matched)

# plot standard curves for each ROI
areas_all %>%
  #mutate(dil = sname2dil[samp], intb = unlist(intb)) %>%
  ggplot(aes(x = dil, y = intb)) +
  facet_wrap(facets = vars(roi), nrow = 6, ncol = 7) +
  geom_point() +
  geom_smooth(method = "lm") +
  stat_cor(size = 3, aes(label = ..rr.label..)) +
  theme_pubr() +
  scale_y_continuous(labels = function(x) format(x, scientific = TRUE)) +
  ggtitle("peak area standard curves by ROI#") +
  xlab("dilution factor") +
  ylab("baseline-adjusted area")

# get the most robust, but sub-saturated spectrum from each ROI

# try to identify these standard spectra using a local EI-MS database
ids_roi <- spectra_master %>%
  group_by(roi) %>%
  annot_from_db("/Users/jwinnikoff/Documents/MBARI/Lipids/GCMSData/db/supel37/")

# check the identifications yourself and make any needed corrections!

# finally, save standard spectra to a new database of authentic standards

## old stuff as of 20200406


# match peaks by RT window and cosine score
# instead of just RT, since integration will be done on base peak
# using loop instead of lapply here to prevent overflow
# obviously the disk reads take awhile

peaks_matched <- NULL
rt_stds <- peaks_stds %>% arrange(rt) %>% pull(rt)
for (i in seq(length(rt_stds))){
  peaks_matched <- match_peaks(
    raw_data,
    std_index = stds_index,
    std_rt = rt_stds[i],
    cosine = 0.9,
    method = "best",
    bin = T,
    rtime_tbl = rtime_tbl
  ) %>%
    bind_rows(peaks_matched, .)
  message(paste("(", i, "/", nrow(peaks_stds), ") matched peaks around ", rt_stds[i], " s", sep = "", collapse = ""))
}

# ggplot the detected *matched* peaks
# NTS: singletons (detected in one dilution only) are being dropped from the matched tibble.
master_coords.m %>%
  filter(samp %in% c("Supel37_1_8", "Supel37_1_8_2")) %>%
  ggplot() +
  geom_line(aes(x = rt, y = intensity, color = samp)) +
  geom_point(data = peaks_matched, aes(x = rt, y = intensity, color = samp), shape = 4) +
  #geom_point(data = peaks_matched %>% select(rt_std) %>% unique(), aes(x = rt_std, y = rep(0, 37)), shape = 25) +
  theme_pubr() +
  theme(legend.position = "right")

## IDENTIFY STANDARDS FROM DATABASE ##
# in principle, these identifications should be 1:1 and unique
# as of 20200323, correct identification from spectra still requires some work
peaks_matched <- annot_from_db(peaks_matched,
                               "/Users/jwinnikoff/Documents/MBARI/Lipids/GCMSData/db/supel37/")

## WRITE SPECTRA FROM MASTER SAMPLE##
peaks_matched %>%
  filter(samp == "Supel37_1_8") %>%
  arrange(rt_std) %>%
  pull(spectrum) %>%
  lapply(function(x){writeMgfData(unlist(x))})

writeMgfData(itraqdata,file="itraqdata.mgf",COM="MSnbase itraqdata")

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

