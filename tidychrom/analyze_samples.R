## Step 1: Blanking
## Step 2: Standard Calibration Curves
## STEP 3: RELATIVE QUANTITATION OF SAMPLES

## Read in GCMS runs from subdirectories corresponding to instrument sessions.
## Map all peaks to a master standard mix (run in that session,)
## integrate them, and flag samples with peaks outside LoL.
## Normalize peak areas using the master standard mix, and return a final
## table of molar percentages for identified compounds.

# USER PARAMETERS
# name of standard file in all the subdirectories
file_stds   <- "Supel37_1_8_blanked.mzXML"
# where to save matched peak area data
# autosaves after reading every data file
file_areas_all <- "/Users/jwinnikoff/Documents/MBARI/Lipids/CtenoLipids2020/tidychrom/example_data/areas_all.save"
# where ROI summary (ID, LoL) data is stored
file_scans_best <- "/Users/jwinnikoff/Documents/MBARI/Lipids/CtenoLipids2020/tidychrom/example_data/scans_best.save"

# resolution of the mass analyzer
bin_width_mz <- 1
# number of peaks to identify in the standard mix
n_stds <- 35
# width in sec of regions of interest (ROIs) used for matching standard peaks
roi_width <- 4
# minimum acceptable cosine similarity to consider two spectra matching
cos_min <- 0.9
# minimum acceptable area of all matched peaks (lower limit QC)
area_matched_min = 2E5

# init areas table
areas_all <- NULL
# init QC table
qc <- tibble(
  file = character(),
  status = character(),
  reason = character(),
)
files_too_low <- NULL
files_too_high <- NULL
for(dir_data in c(
  # session subdirectories
  "/Users/jwinnikoff/Documents/MBARI/Lipids/GCMSData/cdf/20200213",
  "/Users/jwinnikoff/Documents/MBARI/Lipids/GCMSData/cdf/20200214",
  "/Users/jwinnikoff/Documents/MBARI/Lipids/GCMSData/cdf/20200215"
)){
  ## list data files in directory
  mzxmls <- list.files(path = dir_data, pattern = "blanked.mzxml", full.names = T)
  # load raw data for the master
  chromdata_master <- file.path(dir_data, file_stds) %>%
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
  plot_bpc_master <- bpc_master %>%
    ggplot() +
      geom_line(aes(x = rt, y = intensity)) +
      geom_point(data = peaks_master, aes(x = rt, y = intensity), color = "red", shape = 25) +
      theme_pubr() +
      theme(legend.position = "none") +
      ggtitle(paste("BPC:", file.path(dir_data, file_stds)))
  print(plot_bpc_master)

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
        cos >= cos_min,
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
    # but only if there _are_ any matched peaks
    # console message is inside the function
    if(nrow(peaks_matched)){
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

      # save the areas dataframe as we go
      areas_all %>%
        save(file = file_areas_all)
    }else{
      message(paste(basename(file), "matched no peaks to standard!"))
      qc <- qc %>%
        add_case(c(
          "file" = file,
          "status" = "lo",
          "reason" = "no matches to standard",
          ))
    }
  } # end subdirectory
} # end all subdirectories

## Post-analysis: QC and area normalization
# get LoL data, if it's not already loaded
load(file_scans_best)

# Quality Control:
# Lower limit - total matched peak area threshold:
qc <- areas_all %>%
  group_by(file) %>%
  summarise(intb_tot = sum(intb)) %>%
  filter(intb_tot < area_matched_min) %>%
  transmute(
    file = file,
    status = "lo",
    reason = paste("matched area: ", intb_tot,"<", area_matched_min, sep="")
  ) %>%
  bind_rows(qc)

# Upper limit - all peaks must be within limit of linearity
qc <- areas_all %>%
  left_join(scans_best %>% select(-intb, -file), by = "roi") %>%
  filter(intb > intb_max) %>%
  group_by(file) %>%
  summarise(
    status = "hi",
    reason = paste("ROIs", paste(roi, sep=", "),"are saturated")
  ) %>%
  bind_rows(qc)

# filter out samples that fail QC
areas_all_qc <- areas_all %>%
  filter(!(file %in% qc$file))

# normalize peak areas to molar percentages
# NTS: to make more rigorous, need to calculate molar.per.area using each session's standard
areas_all_qc <- areas_all_qc %>%
  left_join(
    scans_best %>%
      mutate(molar.per.area = conc_mol.per.L / intb) %>%
      select(roi, id, molar.per.area),
    by = "roi"
  ) %>%
  group_by(file) %>%
  mutate(
    molar = molar.per.area * intb,
    molar.per.cent = 100 * molar / sum(molar)
      )
