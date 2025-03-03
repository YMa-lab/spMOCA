# spMOCA

## Spatially Informed Matrix Normal Model for Gene Co-expression Analysis in Spatial Transcriptomics

spMOCA is a spatially informed statistical model developed to infer gene co-expression networks across spatial locations in a biologically meaningful way. spMOCA builds upon the principle that gene expression in spatial transcriptomics data is influenced by both spatial dependency and gene-gene interaction. Using a count matrix as input, spMOCA integrates these elements to yield a more accurate and nuanced understanding of gene co-expression networks with spatial contexts. spMOCA employs an efficient optimization algorithm for estimating gene-gene correlation, which is scalable to datasets with tens of thousands of spatial locations and tens of thousands of genes, surpassing the capabilities of existing methods. 

Installation
------------
You can install the released version of spMOCA from Github with the following code, for more installation details or solutions that might solve related issues (specifically MacOS system) see the [link]().

## Dependencies 
* R version >= 4.2.2.
* R packages: Matrix, ggplot2, dplyr, sf, stats, reshape2, gtools, RcppArmadillo, Rcpp, SPARK, WGCNA, msigdbr, tidyverse

``` r
# install devtools if necessary
install.packages('devtools')

# install the spMOCA package
devtools::install_github('YMa-lab/spMOCA')

# load package
library(spMOCA)

```
The R package has been installed successfully on Operating systems: 
* MAC: Sequoia 15.0.1
* Linux: RedHat7

# Issues
All feedback, bug reports and suggestions are warmly welcomed! Please make sure to raise issues with a detailed and reproducible exmple and also please provide the output of your sessionInfo() in R! 


How to use `spMOCA`
-------------------
Details in [Tutorial](https://yma-lab.github.io/spMOCA/)

