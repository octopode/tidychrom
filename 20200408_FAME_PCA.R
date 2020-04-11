library(ggpubr)
library(factoextra)
library(gridExtra)
library(phytools)
library(tidyverse)

palette.deepc <- colorRampPalette(c(
  "#E9EAC6",
  "#3877A4",
  "#8B7CB3",
  "#BF2E55"
))

colon2fa <- Vectorize(function(cpd){
  if(str_detect(cpd, ":0")) {return("SFA")}
  if(str_detect(cpd, ":1")) {return("MUFA")}
  if(str_detect(cpd, ":[2-9]")) {return("PUFA")}
  else {return("other")}
})

## TSV handling, borrowed from 20190201_FAMEs_facet.R but tidied somewhat
# load data with header from second row of the TSV and then drop unnamed columns
data_file <- "/Users/jwinnikoff/Documents/MBARI/Lipids/CtenoLipids2020/tidychrom/example_data/20200408_LipidExtracts.tsv"
num_headerrows = 3 # number of rows with header-type info
header_rownum = 3 # row to use for R headers

headers = read_tsv(data_file, skip = header_rownum-1, col_names = T, n_max = 1)
fame_data = read_tsv(data_file, skip = num_headerrows, col_names = F)
colnames(fame_data) = colnames(headers)

envi_data <- fame_data %>%
  select(
    eid,
    depth_med,
    depth_col,
    temp_med,
    temp_col
  ) %>%
  filter(!is.na(eid))

# taking areas_all_qc from analyze_samples.R in tidychrom
composition_wide <- areas_all_qc %>%
  ungroup() %>%
  select(samp, id, molar.per.cent) %>%
  spread(key = "id", value = "molar.per.cent") %>%
  replace_na(setNames(as.list(rep(0, ncol(.))), colnames(.)))

# run PCA
fame_pca <- composition_wide %>%
  select(-samp) %>%
  prcomp(center = TRUE, scale = TRUE)

envi_data <- composition_wide %>%
  left_join(envi_data, by = "samp")

palette.deepc <- colorRampPalette(c(
  "#E9EAC6",
  "#3877A4",
  "#8B7CB3",
  "#BF2E55"
))

# which PCs correlate best with depth and temperature?
cor(fame_pca$x, envi_data$depth_col, use="complete.obs") %>% tibble() %>% # depth corrs
  bind_cols(cor(fame_pca$x, envi_data$temp_col, use="complete.obs") %>% tibble()) %>% # temp corrs
  bind_cols(fame_pca$sdev^2 %>% tibble()) %>% # eigenvalues
  setNames(c("depth", "temp", "eigen")) %>%
  mutate(eigen = eigen / sum(eigen)) %>% # properly normalize the eigenvalues
  # add the PC#s as a column
  mutate(PC = seq(length(colnames(fame_pca$x)))) %>%
  gather("param", "weight", -PC) %>%
  mutate(param = factor(.$param, levels = c("eigen", "depth", "temp"))) %>%
  filter(PC <= 4) %>% # only show first n PCs
  ggplot(aes(x = PC, y = weight, fill = param)) +
  geom_bar(position="dodge", stat="identity") +
  scale_x_continuous(breaks = seq(1, length(colnames(fame_pca$x)), by = 1)) +
  discrete_scale("fill", "manual", palette.deepc) +
  ggtitle("eigenvalues and correlation strengths\nacross principal components") +
  theme_pubr() +
  theme(
    panel.background = element_rect(fill = "transparent",colour = NA),
    plot.background = element_rect(fill = "transparent",colour = NA)
  ) +
  theme(legend.position = "right")

# biplot with just vectors
fviz_pca_var(fame_pca, repel = T, col.var = "darkgrey") %>%
  fviz_add(., 1*cor(envi_data[c("depth_col")],
                    fame_pca$x, use="complete.obs"),
           color ="#6278AB", geom="arrow", labelsize = 5, linetype = "solid") %>%
  fviz_add(., 1*cor(envi_data[c("temp_col")],
                    fame_pca$x, use="complete.obs"),
           color ="#BE2D55", geom="arrow", labelsize = 5, linetype = "solid") +
  theme_pubr() +
  ggtitle("fatty acid PC loadings\ndepth and temp supplementary")
