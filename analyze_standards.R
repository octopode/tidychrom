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
dir_data  <-  "/Users/jwinnikoff/Documents/MBARI/Lipids/GCMSData/cdf/20200212/stds"
# search pattern for blanked data files
mzxmls <- list.files(path = dir_data, pattern = "blanked.mzxml", full.names = T)
# location of EI-MS database
#dir_db <-     "/Users/jwinnikoff/Documents/MBARI/Lipids/GCMSData/db/supel37/MoNA" # downloaded spectra
dir_db <-     "/Users/jwinnikoff/Documents/MBARI/Lipids/GCMSData/db/supel37/20200414_JRW" # my own library
# location of standard mix datasheet (as TSV)
file_coa <-   "example_data/Supel37_DB23.tsv"
# where to save ROI summary data
file_scans_best <- file.path(dir_data, "scans_best.RData")

# resolution of the mass analyzer
bin_width_mz <- 1
# index in mzxmls of the "master run" used to identify ROIs for the standards
index_master <- 6
# number of peaks to identify in the standard mix
n_stds <- 35
# width in sec of regions of interest (ROIs) used for matching standard peaks
roi_width <- 4
# minimum acceptable cosine similarity to consider two spectra matching
cos_min <- 0.8
# minimum acceptable R^2 value for calibration curves
rsq_min <- 0.95

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
    geom_text(data = peaks_master, aes(x = rt, y = intensity, label = roi), size = 2, position = position_nudge(y = 1E5)) +
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
    read_tidymass()# %>%
    # not required for pre-binned data
    #bin_spectra(bin_width = bin_width_mz)
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

# add standard dilution values to the area table
areas_all <- areas_all %>%
  mutate(dil = filename2dil[file], intb = unlist(intb))

# plot standard curves with all data for each ROI
# overlay R^2 on the plot
areas_all %>%
  # calculate and join R-squareds
  left_join(
    areas_all %>%
      group_by(roi) %>%
      summarize(rsq = r_squared(dil, intb, force_origin = T)),
    by = "roi"
  ) %>%
  ggplot(aes(x = dil, y = intb)) +
  facet_wrap(facets = vars(roi), nrow = 5, ncol = 7) +
  geom_point() +
  # note forcing thru origin
  geom_smooth(method = "lm", formula = y~x+0) +
  geom_text(size = 3, aes(x = 0.05, y = 1E7, label = round(rsq, 4))) +
  theme_pubr() +
  scale_y_continuous(labels = function(x) format(x, scientific = TRUE)) +
  ggtitle("peak area standard curves by ROI: all data") +
  xlab("dilution factor") +
  ylab("baseline-adjusted area")

# calculate upper Limit of Linearity (LoL) for each ROI
areas_all <- areas_all %>%
  group_by(roi) %>%
  # force origin because data are pre-blanked
  summarise(intb_max = lol(dil, intb, rsq = rsq_min, force_origin = T)) %>%
  # and join it back to the area table so we can filter
  right_join(areas_all, by = "roi")

# plot standard curves with in-range data for each ROI
areas_all %>%
  rowwise() %>%
  # only data for which there _is_ a LoL and it's bigger than the area
  # using a mutate() call keeps out-of-bounds observations as NAs,
  # so the facet_grid will still be full and we can see which ROIs are not quantifiable
  # but ggplot does throw warnings as a result!
  mutate(intb = ifelse(!is.na(intb_max) && (intb <= intb_max), intb, NA)) %>%
  # add in the R^2 values for valid curves
  left_join(
    areas_all %>%
      filter(!is.na(intb_max) & (intb <= intb_max)) %>%
      group_by(roi) %>%
      summarize(rsq = r_squared(dil, intb, force_origin = T)),
    by = "roi"
  ) %>%
  ggplot(aes(x = dil, y = intb)) +
  facet_wrap(facets = vars(roi), nrow = 5, ncol = 7) +
  geom_point() +
  # note forcing thru origin
  geom_smooth(method = "lm", formula=y~x+0) +
  geom_text(size = 3, aes(x = 0.05, y = 1E7, label = round(rsq, 4))) +
  theme_pubr() +
  scale_y_continuous(labels = function(x) format(x, scientific = TRUE)) +
  ggtitle(
    paste(
      "peak area standard curves by ROI: linear range only (rsq >= ",
      rsq_min,
      ")",
      sep = "", collapse = ""
    )
  ) +
  xlab("dilution factor") +
  ylab("baseline-adjusted area")

# get the most robust, but sub-saturated spectrum from each ROI
# get those files and scan numbers
scans_best <- areas_all %>%
  group_by(roi) %>%
  arrange(roi) %>%
  # only scans within LDR, comment out if your standard curves suck
  #filter(intb <= intb_max) %>%
  # only the most robust of those
  filter(intb == max(intb))

# now read in the actual spectra
# *opening only 1 file at a time*
clust <- NULL
spectra_best <- pblapply(
  cl = clust,
  unique(scans_best$file),
  function(file_load){
    spectra_best <- scans_best %>%
      filter(file == file_load) %>%
      # using inner join to keep roi, file metadata
      inner_join(
        # read from data file
        #file.path(dir_data, file) %>%
        file %>%
          read_tidymass(),
        by = "scan"
      ) %>%
      select(-rt.x, -mz.x, -intensity.x) %>%
      rename(rt.y = "rt", mz.y = "mz", intensity.y = "intensity")
  }) %>%
  do.call(rbind, .) %>%
  as_tibble()

# try to identify these standard spectra using a local EI-MS database
# and append the id columns to the scans_best tibble
scans_best <- spectra_best %>%
  group_by(roi) %>%
  annot_from_db(dir_db, cores = detectCores()) %>%
  rename(file = "file_db", cos = "cos_db") %>%
  right_join(scans_best, by = "roi")

# to undo the above step (remove DB mappings):
#scans_best <- scans_best %>% select(-file_db, -cos_db)

# check the identifications yourself and make any needed corrections!
# in this case, I am using the first two database IDs to determine that C4:0 and C6:0 were not detected,
# and thereby assigning ROIs to standards data loaded from this file:
stds_data <- read_tsv(file_coa) #%>%
  #mutate(roi = ifelse(order_elute > 2, order_elute - 2, NA))

scans_best <- scans_best %>%
  left_join(stds_data, by = "roi")

# these data can be saved for later sample analysis, in case the session is cleared
# The important mapping in this dataframe (for QC) is roi:intb_max.
#save(scans_best, file = file_scans_best)
#
# finally, save measured standard spectra to the database of authentic standards
#pbmapply(
#  scans_best$scan,
#  scans_best$file,
#  scans_best$id,
#  FUN = function(scan_pull, file_pull, id){
#    file_out <- file.path(dir_db, paste(id, ".mzXML", sep = ""))
#    spectra_best %>%
#      filter((scan == scan_pull) & (file == file_pull)) %>%
#      write_tidymass(file = file_out)
#    return(file_out)
#  }
#)

# Congrats, you determined the LoL for your standards and created a spectrum database!
# Next, on to Step 3: Relative Quantitation of Samples (analyze_samples.R)

# Generate back-to-back spectrum comparisons of standards vs. library.
scans_best <- scans_best %>%
  mutate(
    b2b = spectra_best %>%
      rename(
        roi = "roi_pull",
        scan = "scan_pull"
      ) %>%
      filter((roi_pull == roi) & (scan_pull == scan)) %>%
      rename(
        roi_pull = "roi",
        scan_pull = "scan"
      ) %>%
      bind_rows(
        file.path(dir_db, file_db) %>%
          read_tidymass()
      ) %>%
      # the loaded file always has scan=1,
      # so this makes sure it's on the bottom
      group_by(as.factor(-1*scan)) %>%
      plot_back2back() %>%
      # overlay some explanatory elements
      c(., list(
        ggtitle(paste(
          "ROI",
          roi,
          "\t\t\t\t",
          file_db,
          "\ncos =",
          round(cos_db, 2),
          "\t",
          id
          ))
      )) %>%
      list()
  )

# to plot individual b2bs from this array:
#ggplot() + scans_best$b2b[4] # e.g. for ROI #4

# to plot out a grid
# extract and plotify
b2b <- lapply(scans_best$b2b, function(x){ggplot() + x})
pdf("20200414_cosineMatches_1_40_newCoA.pdf", width = 50, height = 25)
do.call("grid.arrange", c(b2b, nrow = 5, ncol = 7))
dev.off()
