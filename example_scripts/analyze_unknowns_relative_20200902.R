library(tidychrom)
#library(ggrepel)

## Step 1: Blanking
## Step 2: Calculate Relative Ionization Coefficients
## Step 3: RELATIVE COMPOSITION OF UNKNOWN SAMPLES

## Read in unknown samples,
## extract single-ion chromatograms per the guide table,
## detect and integrate the largest peak within each XIC,
## use relative ionization coefficients to calculate composition of each sample.

# USER PARAMETERS

# cores to use for parallel ops
threads = 8L
# file selection criteria
# location of blanked data files
dir_data  <-  "/Users/jwinnikoff/Documents/MBARI/Lipids/20200901_PressureYeast/mzxml"
# filename patterns and dilutions
# suffix used to identify proper file type
suffix_in <- "blanked.mzxml"

# mass spec parameters
# resolution of the mass analyzer
bin_width_mz <- 1
# local SNR threshold for peak detection
# below this SNR won't get integrated
thres_snr <- 5
# RT tolerance for peak collision
# peaks closer than this RT margin will be assigned to the signal ion with greater intensity
tol_rt <- 0.345 # corresponds to 1 maximal scan period

# END USER PARAMETERS

# search for data files
files_unks <- list.files(path = dir_data, recursive = T, pattern = suffix_in, full.names = T) %>%
  enframe(name = NULL) %>%
  dplyr::rename(file_data = value) %>%
  rowwise()

# load BPCs from all files, for visual inspection later
# files only loaded into memory one at a time (per thread), so takes awhile
bpcs <- files_unks %>%
  group_split() %>%
  pblapply(
    .,
    function(row){
      message(paste("loading", basename(row %>% pull(file_data))))
      row %>%
        pull(file_data) %>%
        read_tidymass() %>%
        group_by(scan) %>%
        filter(intensity == max(intensity)) %>%
        mutate(file_data = row %>% pull(file_data))
    },
    cl = threads
  ) %>%
  do.call(rbind, .)

## EXTRACT all ROI XICs
# files only loaded into memory one at a time (per thread), so takes awhile
xics_unks <- files_unks %>%
  #slice(1:5) %>% #TEST
  pull(file_data) %>%
  pblapply(
    .,
    function(file_data){
      # load a run of an unknown sample
      message(paste("loading", basename(file_data)))
      chromdata <- file_data %>%
        read_tidymass() %>%
        ungroup() %>%
        # make 0s explicit for peak integration
        complete(nesting(scan, rt), mz, fill = list(intensity = 0)) %>%
        group_by(scan, rt)

      # extract all the ROIs determined above
      xics <- chromdata %>%
        inner_join(
          quant_guide %>% select(-rt),
          by = "mz"
        ) %>%
        filter((rt >= rt_min) & (rt <= rt_max)) %>%
        # and label with the file they came from
        mutate(file_data = file_data)

      return(xics)
    },
    cl = threads
  ) %>%
  do.call(rbind, .) %>%
  group_by(file_data, stds_mix, roi)

## DETECT *all* peaks in the ROIs
#NTS 20200907 faster on a single thread
peaks_unks_all <- xics_unks %>%
  # due to internal summarize(), this allows all metadata to be retained in output
  group_by_at(vars(-scan, -rt, -intensity)) %>%
  group_split() %>%
  pblapply(
    .,
    function(xic){
      xic %>% detect_peaks(thres_snr = thres_snr)
    },
    cl = NULL
  ) %>%
  do.call(rbind, .) %>%
  group_by_at(vars(-scan, -rt, -intensity)) %>%
  arrange(desc(intensity))

# ASSIGN one peak to each compound
peaks_unks_one <- peaks_unks_all %>%
  group_by(file_data) %>%
  select(-clust_id) %>%
  group_split() %>%
  #.[1] %>% #TEST
  pblapply(
    .,
    function(scans){
      scans %>%
        # cluster_scans requires each scan to be its own group
        rowwise() %>%
        cluster_scans(thres_drt = tol_rt)
    },
    cl = threads
  ) %>%
  do.call(rbind, .) %>%
  # if the *same peak* (same rt) is called twice
  #rowwise() %>%
  #mutate(drt = abs(rt - mean(rt_min, rt_max))) %>%
  #group_by(file_data, clust_id) %>%
  ## keep the one nearest the center of its ROI
  #filter(drt == min(drt)) %>%
  #select(-drt) %>%
  group_by(file_data, clust_id) %>%
  # then just get the strongest remaining peak in each cluster
  filter(intensity == max(intensity)) %>%
  # finally, regroup by file and compound
  group_by(file_data, stds_mix, roi) %>%
  # select the strongest remaining peak in each ROI
  filter(intensity == max(intensity)) %>%
  select(-clust_id) # clust_id doesn't mean what it used to!

# rewindow XICs to these peaks
xics_unks_one <- peaks_unks_one %>%
  select(-scan, -rt, -mz, -intensity) %>%
  right_join(
    xics_unks %>%
      select(file_data, stds_mix, roi, scan, mz, rt, intensity),
    by = c("file_data", "stds_mix", "roi")) %>%
  filter(rt >= rt_min & rt <= rt_max)

## INTEGRATE rewindowed XICs
areas_unks <- xics_unks_one %>%
  group_by_at(vars(-scan, -rt, -rt_min, -rt_max, -intensity)) %>%
  auc() %>%
  group_by(file_data)

## QC ##
files_bad <- areas_unks  %>%
  # flag saturated peaks
  group_by(file_data) %>%
  mutate(qc = ifelse(any(intb > intb_max), "hi", "ok")) %>%
  # flag low-conc files
  mutate(qc = ifelse(sum(intb) < 2.5E5, "lo", qc)) %>%
  group_by(file_data, qc) %>%
  summarize(intb = sum(intb)) %>%
  filter(qc != "ok") %>%
  select(file_data, qc)

## FILTER MOLAR % COMPOSITION ##
comp_unks <- areas_unks %>%
  # filter out bad files
  #anti_join(files_bad, by = "file_data") %>%
  # manual, in this case:
  filter(!any(
    basename(file_data) == "F000B_blanked.mzxml",
    basename(file_data) == "F500B_blanked.mzxml",
    basename(file_data) == "W500B_blanked.mzxml"
  )) %>%
  # also filter out standards
  filter(!str_detect(file_data, "Cayman1") %>% unlist()) %>%
  rowwise() %>%
  mutate(conc_molar = max(0, intb * ion_coeff)) %>%
  group_by(file_data) %>%
  mutate(frac_molar = conc_molar/sum(conc_molar, na.rm = TRUE))

## INSPECT integrations and assignments
# get all BPCs, XICs, and peaks into one dataframe
chroms_all <- xics_unks %>%
  ungroup() %>%
  bind_rows(bpcs %>% mutate(stds_mix = "BPC")) %>%
  bind_rows(xics_unks_one %>% mutate(integrated = TRUE))

# scale to max non-contam intensity within each file
chroms_all_normd <- chroms_all %>%
  group_by(file_data) %>%
  left_join(
    chroms_all %>%
      group_by(file_data) %>%
      summarize(intensity_max = max(intensity))
  ) %>%
  rowwise() %>%
  mutate(
    intensity = intensity/intensity_max,
    intensity = min(intensity, 1)
  )

# plot the chromatograms (all files)
chroms_all_normd %>%
  ggplot() +
  facet_wrap(~paste(file_data %>% as.character() %>% basename(), intensity_max %>% format(scientific = TRUE, digits = 2)), ncol = 1) +
  geom_polygon(
    aes(
      x = rt,
      y = intensity,
      group = paste(integrated, stds_mix, roi),
      fill = ifelse(integrated, stds_mix, "BPC")
    ),
    alpha = 0.5
  ) +
  geom_line(
    aes(
      x = rt,
      y = intensity,
      group = paste(stds_mix, roi, id), # fix stacking order later
      color = factor(stds_mix, levels = c("Cayman1", "BPC")),
      alpha = stds_mix,
      size = stds_mix
    )
  ) +
  scale_color_manual(guide = FALSE, values = c("BPC" = "black", "Cayman1" = "red")) +
  scale_fill_manual(name = "mix", values = c("BPC" = NA, "Cayman1" = "red")) +
  scale_alpha_manual(guide = FALSE, values = c("BPC" = 0.4, "Cayman1" = 1)) +
  scale_size_manual(guide = FALSE, values = c("BPC" = 0.25, "Cayman1" = 0.5)) +
  theme_pubr() +
  # save it to a GIANT PDF
  ggsave(file="/Users/jwinnikoff/Documents/MBARI/Lipids/20200901_PressureYeast/20200907_unks_inspect.pdf", width=10, height=2 * nrow(files_unks) + 2, units="in", limitsize=FALSE)

