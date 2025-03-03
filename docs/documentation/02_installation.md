---
layout: page
title: Installation
description: ~
---

`spMOCA` is implemented as an R package, which can be installed from GitHub by:

### Dependencies 
* R version >= 4.2.2.
* R packages: Matrix, ggplot2, dplyr, sf, stats, reshape2, gtools, RcppArmadillo, Rcpp, SPARK, WGCNA, msigdbr, tidyverse, shiny



#### 1. Install `devtools` if necessary
```r
install.packages('devtools')
```

#### 2. Install `spMOCA`
```r
devtools::install_github('YMa-lab/spMOCA')
```
#### 3. Load package
```r
library(spMOCA)
```

This package is supported for MAC and Linux. The package has been tested on the following systems:
- MAC: Sequoia 15.0.1
- Linux: RedHat7


