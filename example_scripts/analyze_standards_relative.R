library(tidychrom)

## Step 1: Blanking
## STEP 2: CALCULATE MOLAR IONIZATION RATIOS

## Read in serial dilutions of standard mixes.
## Determine limit of linearity (LoL) for each standard's signal ion.
## Determine molar ionization ratios for the standards' signal ions. These may be used globally for relative quantification.

# USER PARAMETERS

# cores to use for parallel ops
n_cores <- detectCores()
# file selection criteria
# location of blanked data files
dir_data  <-  "/Users/jwinnikoff/Documents/MBARI/Lipids/GCMSData/cdf"
# filename patterns and dilutions
# suffix used to identify proper file type
suffix_in <- "blanked.mzxml"
# the analysis will *only* load files containing the following search patterns
# dilutions on the right are used to build standard curves
filename2dil <- c(
  # Supelco 37-FAME standards mix
  "Supel37_1-8_"     = 1/8,
  "Supel37_1-40_"    = 1/40,
  "Supel37_1-200_"   = 1/200,
  "Supel37_1-1k_"    = 1/1000,
  "Supel37_1-5k_"    = 1/5000,
  "Supel37_1-8_"     = 1/8,
  "Supel37_1-10_"    = 1/10,
  "Supel37_1-12_"    = 1/12,
  "Supel37_1-16_"    = 1/16,
  "Supel37_1-32_"    = 1/32,
  "Supel37_1-64_"    = 1/64,
  "Supel37_1-128_"   = 1/128,
  "Supel37_1-256_"   = 1/256,
  # mix of C20:1-OH (d11), C22:1-OH (d13)
  "fatAlcs_1-8_"     = 1/8,
  "fatAlcs_1-10_"    = 1/10,
  "fatAlcs_1-12_"    = 1/12,
  "fatAlcs_1-16_"    = 1/16,
  "fatAlcs_1-24_"    = 1/24,
  "fatAlcs_1-32_"    = 1/32,
  "fatAlcs_1-64_"    = 1/64,
  "fatAlcs_1-128_"   = 1/128,
  "fatAlcs_1-256_"   = 1/256,
  # mix of C18:4-ME and C22:5-ME
  "xPUFAs_1-8_"     = 1/8,
  "xPUFAs_1-10_"    = 1/10,
  "xPUFAs_1-12_"    = 1/12,
  "xPUFAs_1-16_"    = 1/16,
  "xPUFAs_1-24_"    = 1/24,
  "xPUFAs_1-32_"    = 1/32,
  "xPUFAs_1-64_"    = 1/64,
  "xPUFAs_1-128_"   = 1/128,
  "xPUFAs_1-256_"   = 1/256,
  # BHT-derived contaminants (placeholder)
  "BHT"             = 1/8
)

# datasheets (TSV format) for *each* standards mix
# Names in this list should be contained in the search patterns above, and should map 1:1.
# These files serve to specify (1) elution order and (2) quantities. Elution order varies by instrument and method.
# If you do not yet know elution order, you may need to use dummy CoA files with placeholder cpd IDs and quantities.
files_coa <- c(
  "Supel37" = "example_data/Supel37_DB23_EZ-oleic.tsv",
  "fatAlcs" = "example_data/fatAlcs_DB23.tsv",
  "xPUFAs"  = "example_data/xPUFAs_DB23.tsv",
  "BHT" = "example_data/BHT_contam.tsv"
)

# location of spectrum database
#dir_db <-     "~/Documents/MBARI/Lipids/GCMSData/db/supel37/MoNA" # downloaded spectra
dir_db <-     "/Users/jwinnikoff/Documents/MBARI/Lipids/GCMSData/cdf/20200608/spectra_std" # my own library

# mass spec parameters
# resolution of the mass analyzer
bin_width_mz <- 1
# local SNR threshold for peak detection
thres_snr_bpc <- 10
# local SNR threshold for signal ion integration
thres_snr_xic <- 100
# width in sec of regions of interest (ROIs) used for matching standard peaks
roi_width <- 4
# minimum number of shared ions to consider two spectra matching
thres_ions <- 7
# minimum acceptable cosine similarity to consider two spectra matching
cos_min <- 0.99
# minimum acceptable R^2 value for calibration curves
rsq_min <- 0.95
# hard intensity threshold for XICs
thres_intensity <- 0

# try to ID using spectrum library?
read_spec_db <- FALSE
# save your own best scans to the library?
save_spec_db <- FALSE

# where to save ROI summary data
file_scans_best <- file.path(dir_data, "20200810_scans_best.RData")

# END USER PARAMETERS

# search for standard data files
files_stds <- list.files(path = dir_data, recursive = T, pattern = suffix_in, full.names = T) %>%
  enframe(name = NULL) %>%
  dplyr::rename(file_data = value) %>%
  rowwise() %>%
  # filename is in names(filename2dil)
  filter(any(lapply(names(filename2dil), function(x){str_detect(file_data, x)}))) %>%
  #filter(str_count(filename, "20200608") == 1) %>% #TEST session filter
  # use filename pattern matching to identify the standards mix
  mutate(stds_mix = names(files_coa)[min(which(str_detect(file_data, names(files_coa)) == TRUE))]) %>%
  # and likewise to get the dilution
  mutate(dil = filename2dil[min(which(str_detect(file_data, names(filename2dil)) == TRUE))])
# ^beware, min() gets only the first match!

# import and collate the CoAs
coa_all <- lapply(
  names(files_coa),
  function(x){
    read_tsv(files_coa[x]) %>%
      # additional column identifying the standards mix
      mutate(stds_mix = x)
  }) %>%
  do.call(rbind, .) %>%
  group_by(stds_mix)

## DETECT peaks and extract spectra ##
# iterate over rows of the files table
# NOTE: If memory begins to overflow, there is probably something
# wrong with the files_stds tbl. Check it!
# it is ESSENTIAL to do garbage collection before running forked clusters
gc(full = TRUE)
peak_spectra_all <- files_stds %>%
  group_split() %>%
  pblapply(
    .,
    function(row){
      # load a run of a standards mix
      chromdata <- row %>%
        pull(file_data) %>%
        read_tidymass() %>%
        group_by(scan, rt)

      # get the appropriate CoA for that run
      coa_run <- coa_all %>%
        inner_join(row %>% select(stds_mix), by = "stds_mix")

      # how many peaks is it supposed to have?
      n_rois <- nrow(coa_run)

      # extract BPC for peak detection
      bpc <- chromdata %>%
        filter(intensity == max(intensity)) %>%
        ungroup() # important, or detect_peaks will not work!

      # detect peaks
      peaks_bpc <- bpc %>%
        detect_peaks(
          thres_snr = thres_snr_bpc, # user param
          thres_intensity = thres_intensity, # user param
          n = n_rois # derived from CoA
        ) %>%
        group_by(scan, rt)

      # extract peak spectra
      peak_spectra <- chromdata %>%
        inner_join(peaks_bpc %>% select(-mz, -intensity), by = c("scan", "rt")) %>%
        # label with peak filename
        mutate(file = row %>% pull(file_data))

      if(nrow(peak_spectra) > 0){
        return(peak_spectra)
      }else{
        # can't be a warning, because the warning gets returned
        message(paste(basename(row %>% pull(file_data)), "is empty!"))
      }
    },
    cl = 8L
  ) %>%
  do.call(rbind, .) %>%
  group_by(scan, rt, file)

## CLUSTER the features ##
# first by RT proximity (closer than roi_width?)
# then split up each of those peak clusters by cosine
peak_spectra_clustered <- peak_spectra_all %>%
  # by RT proximity
  cluster_scans(thres_drt = roi_width) %>%
  # by RT overlap
  #cluster_peaks(thres_lap = 0.9) %>%
  # then by re-cluster each of those clusters by cosine
  dplyr::rename(clust_id_lap = clust_id) %>%
  # remove singletons, save time
  filter(!is.na(clust_id_lap)) %>%
  group_by(clust_id_lap) %>%
  group_split() %>%
  pblapply(., function(tbl){
    tbl %>%
      group_by(scan, file, rt) %>%
      cluster_spectra(
        thres_cos = cos_min,
        thres_pks = thres_ions,
        cores = 1
      ) %>%
      dplyr::rename(clust_id_cos = clust_id)
    }, cl = 4L) %>% # link this up!
  do.call(rbind, .) %>%
  # again, remove singletons
  filter(!is.na(clust_id_cos)) %>%
  # map standard mix to each peak to split clusters
  left_join(files_stds %>% dplyr::rename(file = file_data) %>% select(file, stds_mix, dil), by = "file")

# In this experiment, no standard mixes share common analytes,
# so we can also split clusters that span multiple standard mixes.
# this creates a new clustering variable
# from the intersection of standard mix, RT overlap, and spectral similarity
peak_spectra_clustered$clust_id <- peak_spectra_clustered %>%
  group_by(stds_mix, clust_id_lap, clust_id_cos) %>%
  group_indices()

# pare down the clusters to just n compounds in each mix, choosing the most populous clusters
clusters_best <- peak_spectra_clustered %>%
  # summarize scans
  #group_by_at(setdiff(colnames(.), c("mz", "intensity"))) %>%
  group_by(file, scan, rt, rt_min, rt_max, clust_id, stds_mix) %>%
  summarize() %>%
  # map to number of standards in mix
  left_join(coa_all %>% summarize(n_stds = n()), by = "stds_mix") %>%
  # summarize clusters
  group_by(clust_id, stds_mix, n_stds) %>%
  summarize(
    # count scans per cluster
    n_scans = n(),
    # conservative summaries for assessing overlap
    rt_min = min(rt_min) - 2*roi_width, # added margins 20200811, increased from 6->8 s 20200813
    rt_max = max(rt_max) + 2*roi_width,
    rt = mean(rt)
  ) %>%
  # and get the n_stds most populous clusters
  group_by(stds_mix) %>%
  arrange(desc(n_scans)) %>%
  # first() not required, but suppresses warning
  slice(1:first(n_stds))

# join clusters to the CoA info
clusters_annotated <- clusters_best %>%
  group_by(stds_mix) %>%
  arrange(rt) %>%
  mutate(roi = seq(length(rt))) %>%
  ungroup() %>%
  left_join(coa_all, by = c("stds_mix", "roi"))

# extract spectra in the winning clusters
cluster_spectra <- peak_spectra_clustered %>%
  semi_join(clusters_annotated, by = "clust_id")

# get consensus spectra for the clusters
# this is used for deconvolution and library construction downstream
cluster_spectra_consensus <- cluster_spectra %>%
  normalize_spectra() %>%
  group_by(clust_id) %>%
  average_spectra()

# ascertain different clusters that coelute
clusters_conflicts <- clusters_annotated %>%
  # cluster_peaks() does a pairwise check of RT overlap
  dplyr::rename(clust_id_lap_cos = clust_id) %>%
  group_by(clust_id_lap_cos) %>%
  # for some reason, does not work with an ungrouped tbl -20200627
  cluster_scans(thres_drt = 6) %>% # maverick conflict detection: < 6 s between peaks
  #cluster_peaks(thres_lap = 0.2) %>% # conservative " " : 20% RT overlap
  # set colnames right again
  dplyr::rename(
    conflict = clust_id,
    clust_id = clust_id_lap_cos
  )

# find distinct signal ions for conflicting clusters
conflicts_resolved <- cluster_spectra_consensus %>%
  inner_join(
    clusters_conflicts %>%
      filter(!is.na(conflict)),
    by = "clust_id"
  ) %>%
  group_by(conflict) %>%
  group_split() %>%
  lapply(.,
    function(tbl){
      tbl %>%
        group_by(clust_id) %>%
        separate_signals(
          #NTS 20200627: pass these parameters from above!
          # if no ions meet these criteria, none are returned!
          thres_ortho = 0.9,
          thres_intensity = 0.1
        )
    }
    ) %>%
  do.call(rbind, .)

# determine appropriate signal ion for each cluster (analyte)
roi_signal_ions <- clusters_conflicts %>%
  left_join(
    # by default, use the most common basepeak m/z within the cluster
    peak_spectra_clustered %>%
      filter(intensity == max(intensity)) %>%
      group_by(clust_id) %>%
      summarize(mz = fmode(mz)),
    by = "clust_id"
  ) %>%
  # if a conflict has been resolved, insert the resolved m/z
  left_join(
    conflicts_resolved %>%
      select(clust_id, conflict, mz),
    by = c("clust_id", "conflict"),
    suffix = c("", ".resolved")
  ) %>%
  mutate(mz = ifelse(is.na(mz.resolved), mz, mz.resolved)) %>%
  select(-mz.resolved) %>%
  # assign ROIs from the CoAs!
  dplyr::rename(stds_mix_clust = stds_mix) %>%
  group_by(stds_mix_clust) %>%
  arrange(rt) %>%
  mutate(
    roi = coa_all %>%
      filter(stds_mix == dplyr::first(stds_mix_clust)) %>%
      pull(roi)
  ) %>%
  group_by(clust_id) %>%
  dplyr::rename(stds_mix = stds_mix_clust)

# set up palette for visualization of XIC integrals
# alternating palette accentuates adjacent ROIs
pal_roi <- RColorBrewer::brewer.pal(name="Dark2", 3)[1:2]

## INTEGRATE appropriate XICs ##
# it is ESSENTIAL to do garbage collection before running forked clusters
gc(full = TRUE)
areas_stds <- files_stds %>%
  #slice(-114) %>% # that one's problematic; pull it out of the pool!
  #slice(115:134) %>% #TEST
  # no need to quantify contaminants
  filter(stds_mix != "BHT") %>%
  group_split() %>%
  pblapply(
    .,
    function(row){
      # load a run of a standards mix
      data_filename <- row %>% pull(file_data)

      message(paste("loading", basename(data_filename)))
      chromdata <- data_filename %>%
        read_tidymass() %>%
        ungroup() %>%
        # make 0s explicit for peak integration
        complete(nesting(scan, rt), mz, fill = list(intensity = 0)) %>%
        #fill_spectra() #unfinished!
        group_by(scan, rt)

      # get all the BPC peaks and store their RTs, in order
      # this will be used to target the integration below
      scans_bpc <- peak_spectra_all %>%
        filter(file == data_filename) %>%
        group_by(scan, rt) %>%
        summarize() %>%
        arrange(rt)

      # extract all the ROIs determined above
      xics <- chromdata %>%
        inner_join(
          roi_signal_ions %>%
            # leave out analytes that definitely aren't there
            filter(stds_mix == row %>% pull(stds_mix)) %>%
            select(-rt),
          by = "mz"
        ) %>%
        filter((rt >= rt_min) & (rt <= rt_max)) %>%
        group_by(clust_id, roi, mz)

      # get bounds of the tallest peak in ROI
      # there's some code in here to tease apart pairs of ROIs that capture 2 peaks each.
      xic_peaks <- xics %>%
        # grouping very important here!
        group_by(clust_id, stds_mix, roi, mz) %>%
        # get bounds of top 3 peaks
        detect_peaks(thres_snr = thres_snr_xic, n=3) %>%
        # and then get the one nearest the BPC peak
        group_by(roi) %>%
        mutate(
          rt_target = scans_bpc %>% pull(rt) %>% .[[first(roi)]],
          drt = abs(rt - rt_target)
        ) %>%
        filter(drt == min(drt))

      #View(xic_peaks %>% select(roi, rt, mz, snr, intensity)) #TEST

      # rewindow each ROI to the largest peak therein
      xics_rewindowed <- xics %>%
        group_by(clust_id, stds_mix, roi, mz) %>%
        select(-rt_min, -rt_max) %>%
        inner_join(
          xic_peaks %>%
            group_by(clust_id, stds_mix, roi, mz) %>%
            select(clust_id, stds_mix, roi, mz, rt_min, rt_max),
          by = c("mz", "clust_id", "stds_mix", "roi")
        ) %>%
        filter((rt >= rt_min) & (rt <= rt_max))

      # There may be _no_ signal ions in the file!
      # Conditional prevents it from crashing the workflow.
      if(nrow(xics) > 0){
        # integrate the xics
        areas <- xics %>%
          auc() %>%
          mutate(x = 1) %>%
          # and bind the run info to the areas
          full_join(row %>% mutate(x = 1), by = "x") %>%
          select(-x) %>% # from here down added for viz
          ## Add visualized XICs to output
          # join ROI XICs
          left_join(
            xics %>%
              group_by(stds_mix, roi, mz) %>%
              summarize(xic_roi = geom_line(data = tibble(rt, intensity), aes(x = rt, y = intensity), color = pal_roi[[mod(first(roi),2)+1]]) %>% list() %>% list()),
            by = c("stds_mix", "roi", "mz")
          ) %>%
          # join peak XICs
          left_join(
            xics_rewindowed %>%
              group_by(stds_mix, roi, mz) %>%
              summarize(xic_peak = geom_polygon(data = tibble(rt, intensity), aes(x = rt, y = intensity), fill = pal_roi[[mod(first(roi),2)+1]], alpha = 0.4) %>% list() %>% list()),
            by = c("stds_mix", "roi", "mz")
          ) %>%
          # combine XICs
          mutate(xic = list(xic_peak, xic_roi) %>% list()) %>%
          select(-xic_roi, -xic_peak)

        return(areas)
      }else(
        # can't be a warning, because the warning gets returned
        message(paste(basename(row %>% pull(file_data)), "contains no signal ions!"))
      )
    },
    cl = 8L
  ) %>%
  do.call(rbind, .)

# visual inspection of integrated peaks by file
# view in a huge, long PDF
areas_stds %>%
  select(file_data, xic) %>%
  group_by(file_data) %>%
  group_split() %>%
  pblapply(., function(tb){
    ggplot() +
      tb %>% pull(xic) +
      theme_pubr() +
      ggtitle(tb %>% slice(1) %>% pull(file_data))
  },
  cl = NULL) %>%
  # takes awhile
  arrangeGrob(grobs=., ncol = 1) %>%
  ggsave(file="~/Downloads/20200810_stds_inspect.pdf", plot=., width=10, height=2 * areas_stds %>% pull(file_data) %>% unique() %>% length() + 2, units="in", limitsize=FALSE)

# visual inspection of each ROI on default graphics device
areas_stds %>%
  select(stds_mix, roi, xic) %>%
  group_by(stds_mix, roi) %>%
  group_split() %>%
  pblapply(., function(tb){
    ggplot() +
      tb %>% pull(xic) +
      theme_pubr() +
      ggtitle(paste(tb %>% slice(1) %>% pull(stds_mix), "ROI", tb %>% slice(1) %>% pull(roi)))
  },
  cl = NULL) %>%
  # takes awhile
  grid.arrange(grobs=., nrow=5, ncol=8)

## DETERMINE LDRs and MOLAR IONIZATION RATIOS ##

# begin by identifying the ROIs
areas_stds_annotated <- areas_stds %>%
  left_join(clusters_annotated %>% select(-clust_id, -rt_min, -rt_max), by = c("stds_mix", "roi")) %>%
  # and then the sessions (days) in which each run was performed
  group_by(file_data) %>%
  mutate(
    # data file's parent directory
    session = file_data %>% dirname() %>% basename(),
    # underscored suffix on the filename
    series = file_data %>%
      basename() %>%
      strsplit(split="_") %>%
      unlist() %>%
      .[[length(.)-1]] %>%
      as.integer()
    ) %>%
  group_by(stds_mix, roi, session, series) %>%
  # calculate max area for each compound in each series
  mutate(intb_max = lol(dil, intb, rsq = rsq_min, force_origin = FALSE, min_dils = 3))

# truncate each series at its saturation point
areas_stds_linear <- areas_stds_annotated %>%
  ungroup() %>%
  filter(intb <= intb_max) %>%
  # calculate ionization ratio in fraction (%) area per fraction molar concentration
  group_by(session, series, dil) %>%
  mutate(
    ion_coeff = (intb / sum(intb)) / (conc_molar / sum(conc_molar))
    )

# normalize the areas to the largest in that mix, ROI, and session
# enables superposition on a facet for assessment of rsq_min
areas_stds_norm <- areas_stds_linear %>%
  ungroup() %>%
  filter(intb >= 0) %>% # filter or force to 0?
  mutate(dil_norm = dil / max(dil)) %>%
  group_by(stds_mix, roi, session, series) %>%
  mutate(
    # the first part of the expression should account for truncated series
    into_norm = dil_norm * into / max(into),
    intb_norm = dil_norm * intb / max(intb)
  )

# plot normalized standard curves for each ROI
# these are subsaturating dilutions only, as determined by rsq_min
# combined R^2 in the facet title
areas_stds_norm %>%
  group_by(stds_mix, roi) %>%
  mutate(
    # calculate R-squareds
    #rsq = r_squared(dil, intb_norm, force_origin = TRUE),
    rsq = r_squared(dil, intb_norm, force_origin = FALSE),
    label = paste("ROI", roi, id, "rsq =", round(rsq, 3))
    ) %>%
  # all these acrobatics are apparently required to make it an ordered factor
  group_by(1) %>% mutate(label = as.factor(label)) %>% ungroup() %>% select(-`1`) %>%
  mutate(label = fct_reorder(label, rt)) %>%
  ggplot(aes(x = dil, y = intb_norm)) +
  facet_wrap(facets = ~label, nrow = 5, ncol = 8) +
  geom_point(aes(color = session), alpha = 0.4) +
  # note forcing thru origin
  geom_smooth(method = "lm", formula = y~x, color = "black") +
  theme_pubr() +
  #scale_y_continuous(labels = function(x) format(x, scientific = TRUE)) +
  ggtitle(paste("normalized peak area standard curves by ROI: series Rsq >=", rsq_min)) +
  xlab("dilution factor") +
  ylab("normalized baseline-adjusted area") +
  scale_color_brewer(palette = "Set3")

# plot molar ionization coefficients
# the coeff for each compound should show a tight central tendency
areas_stds_norm %>%
  mutate(label = paste(stds_mix, roi, id)) %>%
  # all these acrobatics are apparently required to make it an ordered factor
  group_by(1) %>% mutate(label = as.factor(label)) %>% ungroup() %>% select(-`1`) %>%
  mutate(label = fct_reorder(label, rt)) %>%
  ggplot(aes(x = label, y = ion_coeff)) +
    geom_point(aes(color = session), alpha = 0.4) +
    geom_boxplot(outlier.alpha = 0) +
    scale_color_brewer(palette = "Set3") +
    xlab("standard analyte") +
    ylab("relative molar ionization coefficient") +
    theme_pubr() +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))

# finally, save an quantification guide table for each analyte
quant_guide <- roi_signal_ions %>%
  left_join(
    areas_stds_norm %>%
      group_by(stds_mix, roi) %>%
      summarize(
        ion_coeff = median(ion_coeff), # there is some skew, so use median
        # does the saturation point vary between sessions?
        intb_max = max(intb_max)
        ),
    by = c("stds_mix", "roi")
  ) %>%
  # again, no reason to quant contaminants
  #filter(stds_mix != "BHT") %>%
  group_by(id)

## finally, save measured standard spectra to the database of authentic standards
#if(save_spec_db){
#  pbmapply(
#    scans_best_id$scan,
#    scans_best_id$file,
#    scans_best_id$id,
#    FUN = function(scan_pull, file_pull, id){
#      file_out <- file.path(dir_db, paste(id, ".mzXML", sep = ""))
#      spectra_best %>%
#        filter((scan == scan_pull) & (file == file_pull)) %>%
#        write_tidymass(file = file_out)
#      return(file_out)
#    }
#  )
#}
#
## Congrats, you determined the LoL for your standards and created a spectrum database!
## Next, on to Step 3: Relative Quantitation of Samples (analyze_samples.R)
#
## Generate back-to-back spectrum comparisons of standards vs. library.
##NTS 20200507: eliminate variable overwriting outside of loops!
#scans_best_id <- scans_best_id %>%
#  mutate(
#    b2b = spectra_best %>%
#      rename(
#        roi = "roi_pull",
#        scan = "scan_pull"
#      ) %>%
#      filter((roi_pull == roi) & (scan_pull == scan)) %>%
#      rename(
#        roi_pull = "roi",
#        scan_pull = "scan"
#      ) %>%
#      bind_rows(
#        file.path(dir_db, file_db) %>%
#          read_tidymass()
#      ) %>%
#      # the loaded file always has scan=1,
#      # so this makes sure it's on the bottom
#      group_by(as.factor(-1*scan)) %>%
#      plot_back2back() %>%
#      # overlay some explanatory elements
#      c(., list(
#        ggtitle(paste(
#          "ROI",
#          roi,
#          "\t\t\t\t",
#          file_db,
#          "\ncos =",
#          round(cos_db, 2),
#          "\t",
#          id
#        ))
#      )) %>%
#      list()
#  )
#
## to plot individual b2bs from this array:
##ggplot() + scans_best$b2b[4] # e.g. for ROI #4
#
## to plot out a grid
## extract and plotify
#b2b <- lapply(scans_best_id$b2b, function(x){ggplot() + x})
#pdf("20200513_cosineMatches_1_8.pdf", width = 50, height = 25)
#do.call("grid.arrange", c(b2b, nrow = 5, ncol = 7))
#dev.off()
#
