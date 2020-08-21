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
n_cores <- detectCores()
# file selection criteria
# location of blanked data files
dir_data  <-  "/Users/jwinnikoff/Documents/MBARI/Lipids/GCMSData/cdf"
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

# 20200815 minor correction to quant guide to rectify misassignment of C18:1(Z) to C18:1(E)
quant_guide <- quant_guide %>%
  mutate(
    ion_coeff = ifelse(id == "C18:1(E)-ME", ion_coeff*3, ion_coeff)
  )

# END USER PARAMETERS

# search for data files
files_unks <- list.files(path = dir_data, recursive = T, pattern = suffix_in, full.names = T) %>%
  enframe(name = NULL) %>%
  dplyr::rename(file_data = value) %>%
  # leave out contaminant files
  filter(str_detect(file_data, "BHT") == FALSE) %>%
  rowwise()# %>%
  #slice(109:120) #TEST

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
    cl = 8L
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
    cl = 8L
  ) %>%
  do.call(rbind, .) %>%
  group_by(file_data, stds_mix, roi)

## DETECT *all* peaks in the ROIs
#NTS 20200815 may not be necessary to multithread!
peaks_unks_all <- xics_unks %>%
  # due to internal summarize(), this allows all metadata to be retained in output
  group_by_at(vars(-scan, -rt, -intensity)) %>%
  group_split() %>%
  pblapply(
    .,
    function(xic){
      xic %>% detect_peaks(thres_snr = thres_snr)
    },
    cl = 8L
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
    cl = 8L
  ) %>%
  do.call(rbind, .) %>%
  # if the *same peak* (same mz and intensity) is called twice
  rowwise() %>%
  mutate(drt = abs(rt - mean(rt_min, rt_max))) %>%
  group_by(file_data, clust_id, intensity) %>%
  # keep the one nearest the center of its ROI
  filter(drt == min(drt)) %>%
  select(-drt) %>%
  group_by(file_data, clust_id) %>%
  # then just get the strongest remaining peak in each cluster
  filter(intensity == max(intensity)) %>%
  # finally, regroup by file and compound
  group_by(file_data, stds_mix, roi) %>%
  # to select the strongest remaining peak in each ROI
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

# convert JRW filename to sample code
filename2samp <- function(filename){
  prefix <- substr(filename, 0, 3)
  samp <- str_split(filename, "[_-]") %>%
    .[[1]] %>% .[1] %>%
    substr(4, str_length(.)) %>%
    str_pad(., 4, "0", side = "left") %>%
    paste(prefix, ., sep="")
  return(samp)
}

## INTEGRATE rewindowed XICs
areas_unks <- xics_unks_one %>%
  group_by_at(vars(-scan, -rt, -rt_min, -rt_max, -intensity)) %>%
  auc() %>%
  mutate(eid = filename2samp(file_data %>% basename())) %>%
  group_by(eid) %>%
  mutate(count = n())

## QC ##
files_bad <- areas_unks  %>%
  filter(stds_mix != "BHT") %>%
  # flag saturated peaks
  group_by(file_data) %>%
  mutate(qc = ifelse(any(intb > intb_max), "hi", "ok")) %>%
  # flag low-conc files
  # maybe should be 2.5E5. There's quite a lot of gross s**t down there...
  mutate(qc = ifelse(sum(intb) < 1E5, "lo", qc)) %>%
  group_by(file_data, qc) %>%
  summarize(intb = sum(intb)) %>%
  mutate(eid = filename2samp(file_data %>% basename())) %>%
  group_by(eid) %>%
  # an optimistic bit added 20200816
  # for possibly saturated samples, keep the one with lowest total area
  # comment to be conservative!
  mutate(qc = ifelse(qc == "hi" & intb == min(intb), "ok", qc)) %>%
  group_by(file_data, qc) %>%
  filter(qc != "ok") %>%
  select(file_data, qc)

## FILTER MOLAR % COMPOSITION ##
comp_unks <- areas_unks %>%
  # filter out bad files
  anti_join(files_bad, by = "file_data") %>%
  rowwise() %>%
  mutate(conc_molar = max(0, intb * ion_coeff)) %>%
  group_by(file_data) %>%
  mutate(frac_molar = conc_molar/sum(conc_molar, na.rm = TRUE))

## INSPECT integrations and assignments
# get all BPCs, XICs, and peaks into one dataframe
chroms_all <- xics_unks %>%
  ungroup() %>%
  bind_rows(bpcs %>% mutate(stds_mix = "BPC")) %>%
  bind_rows(xics_unks_one %>% mutate(integrated = TRUE)) %>%
  # reorder the stds_mix factor!
  mutate(
    integrated = ifelse(!is.na(integrated), TRUE, FALSE),
    file_data = factor(file_data),
    stds_mix = factor(stds_mix, levels = c("Supel37", "xPUFAs", "fatAlcs", "BHT", "BPC"))
  ) %>%
  #filter((as.numeric(file_data) > 5) & (as.numeric(file_data) <= 10)) %>%  #TEST
  mutate(stds_mix = paste(stds_mix))

# scale to max non-contam intensity within each file
chroms_all_normd <- chroms_all %>%
  group_by(file_data) %>%
  left_join(
    chroms_all %>%
      group_by(file_data) %>%
      filter(!(stds_mix %in% c("BPC", "BHT"))) %>%
      summarize(intensity_max = max(intensity))
  ) %>%
  rowwise() %>%
  mutate(
    intensity = intensity/intensity_max,
    intensity = min(intensity, 1)
  )

# color mapping
chroma <- c(
  "BPC"     = "black",
  "Supel37" = "#1B9E77",
  "xPUFAs"  = "#D95F02",
  "fatAlcs" = "#7570B3",
  "BHT"     = "#E7298A"
)

# plot the chromatograms (all files)
chroms_all_normd %>%
  ggplot() +
  facet_wrap(~paste(file_data, intensity_max %>% format(scientific = TRUE, digits = 2)), ncol = 1) +
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
      color = factor(stds_mix, levels = c("BHT", "Supel37", "xPUFAs", "fatAlcs", "BPC")),
      alpha = stds_mix,
      size = stds_mix
    )
  ) +
  scale_color_manual(guide = FALSE, values = chroma) +
  scale_fill_manual(name = "mix", values = c("BPC" = NA, chroma[2:length(chroma)])) +
  scale_alpha_manual(guide = FALSE, values = c(rep(1, 4), 0.4) %>% setNames(c("Supel37", "xPUFAs", "fatAlcs", "BHT", "BPC"))) +
  scale_size_manual(guide = FALSE, values = c(rep(0.5, 4), 0.25) %>% setNames(c("Supel37", "xPUFAs", "fatAlcs", "BHT", "BPC"))) +
  theme_pubr() +
  # save it to a GIANT PDF
  ggsave(file="~/Downloads/20200815_unks_inspect.pdf", width=10, height=2 * nrow(files_unks) + 2, units="in", limitsize=FALSE)

# plot the chromatograms (QC pass only)
chroms_all_normd %>%
  anti_join(files_bad, by = "file_data") %>%
  ggplot() +
  facet_wrap(~paste(file_data, intensity_max %>% format(scientific = TRUE, digits = 2)), ncol = 1) +
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
      color = factor(stds_mix, levels = c("BHT", "Supel37", "xPUFAs", "fatAlcs", "BPC")),
      alpha = stds_mix,
      size = stds_mix
    )
  ) +
  scale_color_manual(guide = FALSE, values = chroma) +
  scale_fill_manual(name = "mix", values = c("BPC" = NA, chroma[2:length(chroma)])) +
  scale_alpha_manual(guide = FALSE, values = c(rep(1, 4), 0.4) %>% setNames(c("Supel37", "xPUFAs", "fatAlcs", "BHT", "BPC"))) +
  scale_size_manual(guide = FALSE, values = c(rep(0.5, 4), 0.25) %>% setNames(c("Supel37", "xPUFAs", "fatAlcs", "BHT", "BPC"))) +
  theme_pubr() +
  # save it to a GIANT PDF
  ggsave(file="~/Downloads/20200815_unks_inspect_qcpass.pdf", width=10, height=2 * (nrow(files_unks) - nrow(files_bad)) + 2, units="in", limitsize=FALSE)
