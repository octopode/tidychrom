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
file_areas_all <- "example_data/areas_all.RData"
# where ROI summary (ID, LoL) data is stored
#file_scans_best <- "example_data/scans_best.RData"
file_scans_best <- "/Users/jwinnikoff/Documents/MBARI/Lipids/GCMSData/cdf/20200212/stds/scans_best.RData"
# where to save an output file of pct composition
file_tsv_out <- "20200518_pct_wide.tsv"

# resolution of the mass analyzer
bin_width_mz <- 1
# number of peaks to identify in the standard mix
n_stds <- 35
# width in sec of regions of interest (ROIs) used for matching standard peaks
roi_width <- 4
# minimum acceptable cosine similarity to consider two spectra matching
cos_min <- 0.5
# hard intensity threshold for XICs
threshold_quan <- 1000
# minimum acceptable area of all matched peaks (lower limit QC)
area_matched_min = 1E5

# init areas table
areas_all <- NULL
# init QC table
qc <- tibble(
  file = character(),
  status = character(),
  reason = character(),
)

base_path = "~/Documents/MBARI/Lipids/GCMSData/"

for(subdir in c(
  # session subdirectories
   "cdf/20200213",
   "cdf/20200214",
   "cdf/20200215"
)){
  dir_data = paste(base_path,subdir,sep="")
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
        extract_chromatograms(chromdata, threshold = threshold_quan, cores = detectCores())

      # integrate those XICs
      message("integrating XICs")
      peak_areas <- xics_matched %>%
        group_by(roi) %>%
        auc() %>%
        # and append column identifying the file
        mutate(file = basename(file))

      ## multithread plotting
      #message("generating plots")
      #peak_areas <- peak_areas %>%
      #  group_by(roi) %>%
      #  bind_cols(b2b_parallel(.)) %>%
      #  bind_cols(xic_parallel(.))

      message("generating plots")
      peak_areas <- peak_areas %>%
        group_by(roi) %>%
        mutate(
          # plot the spectral matchups
          b2b = chromdata %>%
            # get sample spectrum
            rename(
              scan = "scan_pull"
            ) %>%
            filter(scan_pull == scan) %>%
            rename(
              scan_pull = "scan"
            ) %>%
            bind_rows(
              # get standard spectrum
              spectra_master %>%
                rename(
                  roi = "roi_pull"
                ) %>%
                filter(roi_pull == roi) %>%
                rename(
                  roi_pull = "roi"
                ) %>%
                # and isolate the scan
                mutate(scan = 0)
            ) %>%
            # the loaded file always has scan=0,
            # so this makes sure it's on the bottom
            group_by(as.factor(-1*scan)) %>%
            plot_back2back() %>%
            # overlay some explanatory elements
            c(., list(
              ggtitle(paste(
                "ROI ",
                roi,
                "\tcos = ",
                round(cos, 2),
                "\n",
                basename(file),
                sep = "", collapse = ""))
            )) %>%
            list(),
          # plot the xics
          xic = xics_matched %>%
            # get the current ROI
            rename(
              roi = "roi_pull"
            ) %>%
            filter(roi_pull == roi) %>%
            rename(
              roi_pull = "roi"
            ) %>%
            plot_xic() %>%
            # overlay some explanatory elements
            c(., list(
              ggtitle(paste(
                "ROI ",
                roi,
                "\tcos = ",
                round(cos, 2),
                "\n",
                basename(file),
                sep = "", collapse = ""))
            )) %>%
            list()
        )

      areas_all <- areas_all %>%
        bind_rows(peak_areas %>% ungroup())

      # save the areas dataframe as we go
        save(areas_all, file = file_areas_all)
    }else{
      message(paste(basename(file), "matched no peaks to standard!"))
      qc <- qc %>%
        add_case(
          "file" = basename(file),
          "status" = "lo",
          "reason" = "no matches to standard"
          )
    }
  } # end subdirectory
} # end all subdirectories
# add column with EID
areas_all <- areas_all %>%
  group_by(file) %>%
  mutate(samp = filename2samp(file))

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
    reason = paste("matched area: ", round(intb_tot),"<", area_matched_min, sep="")
  ) %>%
  bind_rows(qc)

# Upper limit - all peaks must be within limit of linearity
qc <- areas_all %>%
  left_join(scans_best_id %>% select(-intb, -file), by = "roi") %>%
  filter(intb > intb_max) %>%
  group_by(file) %>%
  summarise(
    status = "hi",
    reason = paste("ROIs", paste(roi, sep="", collapse = ", "),"are saturated")
  ) %>%
  bind_rows(qc)

# filter out samples that fail QC
areas_all_qc <- areas_all %>%
  filter(!(file %in% qc$file))

# clean up the integrated peak data post hoc
# tune these parameters using windowed XIC plot at bottom
areas_qc_clean <- areas_all_qc %>%
  filter(intb > 0) %>%
  # e.g. cut out inordinately broad peaks
  filter((rt_max - rt_min) <= 12) %>%
  # and remove outlier RTs
  group_by(roi) %>%
  # could also weight by intensity here instead of integral
  filter(abs(rt - fmean(x=rt, w=intb)) <= 4*fsd(x=rt, w=intb))

# here's what you just threw out
areas_qc_trash <- areas_all_qc %>%
  anti_join(areas_qc_clean, by=c("roi", "file"))

# normalize peak areas to molar percentages
# NTS: to make more rigorous, need to calculate molar.per.area using each session's standard
areas_quan <- areas_qc_clean %>%
  left_join(
    scans_best_id %>%
      mutate(molar.per.area = conc_molar * dil / intb) %>%
      select(roi, id, molar.per.area),
    by = "roi"
  ) %>%
  group_by(file) %>%
  mutate(
    molar = molar.per.area * intb,
    molar.per.cent = 100 * molar / sum(molar)
  )

# data export
# save tsv in wide format
areas_quan %>%
  filter(samp != "Supel37") %>%
  ungroup() %>%
  select(samp, id, molar.per.cent) %>%
  spread(key = id, value = molar.per.cent) %>%
  # replace NA -> 0
  mutate_all(function(x){ifelse(is.na(x), 0, x)}) %>%
  rename(samp = "eid") %>%
  write_tsv(file_tsv_out)

## Visualization:
# Spectral matchups and integrated XICs for every single peak are stored
# in areas_all and areas_all_qc. They can be accessed individually like so:
samp_viz <- c("JWL0051")
cpd_viz <- "C20:5"
roi_viz <- 35

ggplot() + areas_quan %>%
  filter(
    samp == samp_viz &
      id == cpd_viz
    ) %>%
  pull(b2b)

# 20200507
ggplot() + areas_quan %>%
  filter(
    samp %in% samp_viz &
    id == cpd_viz
  ) %>%
  pull(xic)

# series of xics can also be overlaid:
# e.g. to show all the ROIs integrated in a given sample
ggplot() +
  areas_quan %>%
  filter(
    samp %in% samp_viz
  ) %>%
  pull(xic) +
  ggtitle(paste("all ROIs: sample", samp_viz))

# 20200507
ggplot() +
  # overlay a BPC
  #geom_line(
  #  data = chromdata %>%
  #    group_by(scan) %>%
  #    filter(rt >= 600) %>%
  #    filter(intensity == max(intensity)),
  #  aes(x = rt, y = intensity)
  #) +
  areas_quan %>%
  filter(
    samp == samp_viz
  ) %>%
  pull(xic) +
  ggtitle(paste("all ROIs: sample", samp_viz))

# or to show all the samples found in a given ROI
ggplot() + areas_qc_clean %>%
  filter(
    id == cpd_viz,
    #samp %in% samp_viz
  ) %>%
  pull(xic) +
  ggtitle(paste(cpd_viz, "in", samp_viz))
# the stored plots are already titled, so the top-level
# title should be overwritten with your own

# this can be used to check peak bounding conditions in extract_chromatograms()
# to plot out a grid
# extract and plotify
#frames <- lapply(areas_all_qc$xic, function(x){ggplot() + x})
frames <- lapply(
  areas_qc_clean %>%
  #areas_qc_trash %>%
    group_by(roi) %>%
    summarize(xic = list(xic)) %>%
    pull(xic),
  function(x){ggplot() + x})
pdf("20200518_xics_qc_clean_4wsd.pdf", width = 50, height = 25)
do.call("grid.arrange", c(frames, nrow = 5, ncol = 7))
dev.off()
