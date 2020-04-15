# tidychrom

A Dead Simple Toolkit for Quantitative Chromatography

## workflows

### targeted relative quantitation

## work-in-progress notes

### guiding design principles

1. **ingredients you can pronounce**

No fancy algorithms <cough> ordered bijective interpolated time warping </cough>, 
though nothing explicitly prevents their use with this package.

**Caveat:** since `tidychrom` does not implement RT adjustment nor spectrum
deconvolution, it expects 

	* chromatographic data with good separation and 
	
	* fairly consistent RTs (within a few sec) across samples.
	
Analyze highly complex mixtures at your own risk, and maybe with a dash of special sauce.

2. **data you can see and touch**

Existing R-based chromatography solutions rely on S3 and S4 objects with slots
that are not really standardized. This package attempts to keep all analysis
products in a single `tibble` and facilitate downstream analysis with `dplyr`.

Visualization is implemented with `ggplot2`.

### project structure

1. Functional programming style avoids custom objects.

2. Functions are packaged 1 to a file.
	
3. Visualization:

All implemented in ggplot. Custom plotting functions return a list of `gg` objects,
which can be:

	+ stored in a tibble and called up later,
	
	+ arranged alongside other stored plots,
	
	+ overlaid with other `gg` elements, like titles and other plots