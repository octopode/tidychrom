# tidychrom

A dead simple toolkit for quantitative chromatography

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

2. Functions are packaged 1 or 2 to a file:

	1. The first function may operate on any combination of arguments, possibly 
	including `.x` and `.y` for use in a `group_map()` call.

	2. The second function should take a `peaks` tibble as an argument, though
	it may require particular columns to be present in that tibble, as documented. 
	Absence of a required column will throw an `object not found` error.