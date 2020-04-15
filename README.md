# tidychrom

A Dead Simple Toolkit for Quantitative Chromatography

![base peak chromatogram](img/20200414_masterBPC%202.png)

## design principles

1. **ingredients you can pronounce**

No fancy algorithms <cough> ordered bijective interpolated time warping </cough>, 
though nothing explicitly prevents their use with this package.

**Caveat:** since `tidychrom` does not implement RT adjustment nor spectrum
deconvolution, it expects

* chromatographic data with good separation and 
	
* fairly consistent RTs (within a few sec) across samples.
	
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
![spectrum matchup](img/20200414_JWL12_DHA_matchup%202.png)
	
+ arranged alongside other stored plots,

(as in [analyze_standards.R](analyze_standards.R))
```
b2b <- lapply(scans_best$b2b, function(x){ggplot() + x})
do.call("grid.arrange", c(b2b, nrow = 5, ncol = 7))
```
![ALL spectrum matchups](img/20200414_cosineMatches_1_40_newCoA%202.png)
	
+ overlaid with other `gg` elements, like titles and other plots.

## workflows

### targeted relative quantitation

This workflow was made to calculate the ratio (molar percentages) of fatty acid methyl esters (FAMEs) in a set of biological samples. It has 3 steps:

1. [Scanwise blanking](subtract_blanks_multidir.R)

2. [Standard ID and saturation point determination](analyze_standards.R)

3. [Sample ID and relative quantitation](analyze_samples.R)

