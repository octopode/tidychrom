## STEP 1: BLANKING

# Preprocessing script to read in blank chromatographic runs, average them scanwise,
# then subtract the average from a set of experimental runs

for(dir_data in c(
  "/Users/jwinnikoff/Documents/MBARI/Lipids/GCMSData/cdf/20200213",
  "/Users/jwinnikoff/Documents/MBARI/Lipids/GCMSData/cdf/20200214",
  "/Users/jwinnikoff/Documents/MBARI/Lipids/GCMSData/cdf/20200215"
)){
  # filter for file type
  cdfs <- list.files(path = dir_data, pattern = ".cdf", full.names = T)
  # blank runs
  files_blanks <- cdfs[which(str_detect(cdfs, "blank"))]
  # experimental runs
  files_exptal <- setdiff(cdfs, files_blanks)

  message(paste("preparing blank(s) in", basename(dir_data)))
  # load the blank runs
  chromdata_blanks <- files_blanks %>%
    # first into a list of tibbles
    lapply(., read_tidymass) %>%
    # with m/z binned to integer
    lapply(., function(x){bin_spectra(x, bin_width = 1)}) %>%
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
  #chromdata_blank_avg %>%
  #  mutate(file = "average") %>%
  #  bind_rows(chromdata_blanks) %>%
  #  group_by(file, scan) %>%
  #  # gets base peak intensity
  #  filter(intensity == max(intensity)) %>%
  #  ggplot() +
  #    geom_line(aes(x = rt, y = intensity, color = file)) +
  #    theme_pubr() +
  #    theme(legend.position = "right")

  # then, looping thru experimental files,
  for(file_in in files_exptal){
    message(paste("subtracting blank from", basename(file_in)))
    chromdata_exptal <- file_in %>%
      read_tidymass() %>%
      bin_spectra(bin_width = 1) %>%
      # subtract the blank scanwise
      left_join(chromdata_blank_avg, by = c("scan", "mz")) %>%
      replace_na(list(intensity.y = 0)) %>%
      mutate(
        rt = rt.x,
        intensity = intensity.x - intensity.y,
        intensity = ifelse(intensity < 0, 0, intensity)
      ) %>%
      select(scan, rt, mz, intensity) %>%
      filter(intensity > 0)

    # and save
    file_out <- str_replace(file_in, "\\.cdf", "_blanked.mzxml")
    chromdata_exptal %>% write_tidymass(file_out)
  }
}
