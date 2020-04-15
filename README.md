# tidychrom

A Dead Simple Toolkit for Quantitative Chromatography

<flowchart>

## design principles

1. **ingredients you can pronounce**

No fancy algorithms <cough> ordered bijective interpolated time warping </cough>, 
though nothing explicitly prevents their use with this package.

**Caveat:** since `tidychrom` does not implement RT adjustment nor spectrum
deconvolution, it expects chromatographic data with (1) good separation and 
(2) fairly consistent RTs (within a few sec) across samples.

Something like this:

<img src="img/20200414_masterBPC%202.png" alt="base peak chromatogram" width="500px">

(master base peak chromatogram from [analyze_standards.R](analyze_standards.R))
	
Analyze highly complex mixtures at your own risk, and maybe with a dash of special sauce.

2. **data you can see and touch**

Existing R-based chromatography solutions rely on S3/4 objects with slots
that are not really standardized. This package attempts to keep all analysis
products in tibbles and facilitate downstream analysis with `dplyr`.

Visualization is implemented with `ggplot2`.

## project structure

This repo contains a couple pre-cooked workflows (see **workflows** below), but above all, tidychrom is meant to be modular. Take the handful of functions provided and use them in your own `dplyr`-based workflows, perhaps with inspiration from those provided here. To aid you in dissecting this repo, some tips on its organization:

1. There are no custom objects nor methods, only functions.

2. Functions are packaged 1 to a file.
	
3. Visualization is implemented in `ggplot2`. Custom plotting functions return a list of `gg` objects,
which can be

+ stored in a tibble column and _retrieved later_,

(from [analyze_samples.R](analyze_samples.R))
```
ggplot() + areas_all_qc %>%
  filter(
    samp == "JWL0012" &
      id == "C22:6"
    ) %>%
  pull(b2b)
```
<img src="img/20200414_JWL12_DHA_matchup%202.png" alt="spectrum matchup" width="350px">
	
+ arranged alongside other stored plots,

(as in [analyze_standards.R](analyze_standards.R))
```
b2b <- lapply(scans_best$b2b, function(x){ggplot() + x})
do.call("grid.arrange", c(b2b, nrow = 5, ncol = 7))
```
<img src="img/20200414_cosineMatches_1_40_newCoA%202.png" alt="ALL spectrum matchups" width="350px">
	
+ overlaid with other `gg` elements, like titles and other plots.

(from [analyze_samples.R](analyze_samples.R))
```
# to show all the ROIs integrated in a given sample
ggplot() + areas_all_qc %>%
  filter(
    samp == "JWL0138"
  ) %>%
  pull(xic) +
  ggtitle("all ROIs: sample JWL0138")
```
<img src="https://github.com/octopode/tidychrom/blob/master/img/20200414_JWL138_allROIs%202.png" alt="JWL0138 all ROIs" width="350px">

```
# or to show all the samples found in a given ROI
ggplot() + areas_all_qc %>%
  filter(
    id == "C22:6"
  ) %>%
  pull(xic) +
  ggtitle("C22:6 (DHA): 39 samples")
```
<img src="https://github.com/octopode/tidychrom/blob/master/img/20200414_DHA_allXICs%202.png" alt="DHA all XICs" width="350px">

## workflows

### targeted relative quantitation

This workflow was made to calculate the ratio (molar percentages) of fatty acid methyl esters (FAMEs) in a set of biological samples. It has 3 steps:

1. [Scanwise blanking](subtract_blanks_multidir.R)

2. [Standard ID and saturation point determination](analyze_standards.R)

3. [Sample ID and relative quantitation](analyze_samples.R)

Comments will walk you through each script.
All user-provided parameters (including data directories) are provided at the top.
Data files used in the scripts will be hosted at a later date.
